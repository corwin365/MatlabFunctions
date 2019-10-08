function FilteredData = weighted_modal_filter(Data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%"weighted modal filter"
%
%finds the mode in the area surrounding each gridbox, where the points in
%that region can be weighted to contribute more to the mode
%
%not sure if it's mathematically valid, but may be useful...
%
%Corwin Wright, c.wright@bath.ac.uk, 26/Sep/2019
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%inputs:
%%%%%%%%%%%
%
 %required
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
%    VarName     (type)   description
%    -----------------------------------------------------------------------------
%    Data        (array)  Data to be filtered, as an ND array
%
%
 %optional (all case-insensitive)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    VarName         (type,                    default)    description
%    -----------------------------------------------------------------------------
%    FilterSize      (array,                 [5,5,...])     size of volume to be smoothed, number of elements same as input data. If odd, will be rounded up
%    Weights         (array, ones of same size as data)     weights to apply
%    Coarseness      (numeric,                     100)     number of points to use for weighting at each input data point - essentially, the granularity of it
%
%
%%%%%%%%%%%
%outputs:
%%%%%%%%%%%
%
%  FilteredData:  array of same size as InputData, containing filtered data
%
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create input parser
%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

%inputs - required
%%%%%%%%%%%%%%%%%%%

%data must just exist and be numeric
addRequired(p,'Data', @isnumeric);

%inputs - optional
%%%%%%%%%%%%%%%%%%%

%filtersize must have same dimensionality as the input data
NDims = numel(size(Data));
CheckFilterSize = @(x) validateattributes(x,{'numeric'},{'size',[1,NDims],'positive'});
addParameter(p,'FilterSize',ones(1,NDims).*5,CheckFilterSize);

%coarseness needs to be a positive number
CheckCoarseness = @(x) validateattributes(x,{'numeric'},{'nonnegative'});
addParameter(p,'Coarseness',100,CheckCoarseness);

%weights must be the same size as the data
sz = size(Data);
CheckWeights = @(x) validateattributes(x,{'numeric'},{'size',sz,'>=',0});
addParameter(p,'Weights',ones(sz),CheckWeights); 

%parse the inputs, and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,Data,varargin{:})

%pull out the contents into struct "Inputs", used throughout rest of routine
Input = p.Results;

clearvars -except Input
Data       = Input.Data;
Weights    = Input.Weights;
FilterSize = Input.FilterSize;
Coarseness = Input.Coarseness;
clear Input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find the input array size. we'll use this a lot.
OriginalSize = size(Data);

%make an array of indices of the same size and shape
Indices = NaN.*Data;
Indices(1:end) = 1:1:numel(Data);

%create an array for the results
FilteredData = repmat(reshape(Data,prod(OriginalSize),1).*NaN,1,Coarseness);


%for each point, work out the weights and hence the filtered value
parfor iPoint = 1:1:numel(Data)
  
  %need to split out the column we're working on, in order for parfor to identify the variable type
  ThisChunk = FilteredData(iPoint,:);
  
  
  %work out where we are in the original array
  %this requires us to use ind2sub, but with a variable number of dimensions
  PointIndices = cell(1, numel(OriginalSize));
  [PointIndices{:}] = ind2sub(OriginalSize,iPoint);
  PointIndices = cell2mat(PointIndices);

  %ok. now we have the indices, work out which points fall into the domain
  Edges = NaN(numel(OriginalSize),2);
  for iDim=1:1:numel(OriginalSize)
    Edges(iDim,1) = PointIndices(iDim) - ceil(FilterSize(iDim)/2);
    Edges(iDim,2) = PointIndices(iDim) + ceil(FilterSize(iDim)/2);
  end
  
  %check they're actually in the domain, and truncate if not
  for iDim=1:1:numel(OriginalSize); Edges(iDim,:) = bound(Edges(iDim,:),1,OriginalSize(iDim)); end

  %ok. Now pull out the indices in this range
  %hardcode up to 3d for speed, then provide a generalised routine for
  %higher dimensionalities which is much slower
  if numel(OriginalSize) == 1
    IndicesInBox = Indices(Edges(1):Edges(2));
  elseif numel(OriginalSize) == 2;
    IndicesInBox = Indices(Edges(1,1):Edges(1,2),Edges(2,1):Edges(2,2));
  elseif numel(OriginalSize) == 3;
    IndicesInBox = Indices(Edges(1,1):Edges(1,2),Edges(2,1):Edges(2,2),Edges(3,1):Edges(3,2));
  else
    %general case
    IndicesInBox = Indices;
    for iDim=1:1:numel(OriginalSize)
      sz = size(IndicesInBox);
      NewOrder = unique([iDim,1:1:numel(sz)],'stable');
      D = permute(IndicesInBox,NewOrder);
      D = reshape(D,[sz(NewOrder(1)),prod(sz(NewOrder(2:end)))]);
      D = D(Edges(iDim,1):Edges(iDim,2),:);
      D = reshape(D,[size(D,1),sz(NewOrder(2:end))]);
      [~,idx] = sort(NewOrder);
      IndicesInBox = permute(D,idx);
    end
  end
  %So. We have the weights and we have the indices of the objects in the
  %window
  %from this, work out how many times each point should occur in the data
  WeightsInBox = Weights(IndicesInBox);
  Sigma = sum(WeightsInBox(:));
  if Sigma > 0;
    
    %what is the relative weight of each point?
    RelativeWeight = WeightsInBox ./ Sigma;
    
    %what is the relative weight of one element in the coarseness space?
    SingleWeight = Sigma./Coarseness;
    
    %loop over the weights and fill an array with the largest, subtracting
    %a single weight each time, until all values are less than a single
    %weight
    for iElement=1:1:Coarseness
      
      %find the heaviest, and remove one unit of weight from it
      [~,Biggest] = max(RelativeWeight(:));
      RelativeWeight(Biggest) = RelativeWeight(Biggest) - SingleWeight;
      
      %store the value of the heaviest
      ThisChunk(iElement) = Data(IndicesInBox(Biggest));
    end
  else
    %if the weights are invalid, then just find the mode of the region
    ThisChunk(:) = mode(Data(IndicesInBox(:))); 
  end
  
  %stick it where it belongs
  FilteredData(iPoint,:) = ThisChunk;
end

%take mode, and done!
FilteredData = mode(FilteredData,2);
FilteredData = reshape(FilteredData,OriginalSize);


%done!

end

