function Error = struct2nc(FileName,Data,Dimensions,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parse a Matlab data struct and produce a netCDF file containing
%the contents
%
%*********************************************************************
%IMPORTANT NOTES
%
%1. Currently the logic of the function REQUIRES that every dimension has
%   a unique number of elements, i.e. if two arrays have any dimension of
%   the same length, then it is the same dimension semantically.
%
%   If this cannot be assumed, then you'll need to use another approach 
%   to save these data.
%
%2. The function will also discard variable types I haven't configured my 
%   underlying functions to use - see 'ValidList' variable for a list of 
%   valid types
%
%3. Complex values can be saved, but will be automatically split into two
%   variables suffixed as '_imag' and '_real'.
%
%
%*********************************************************************
%
%required inputs:
%
%    [char] FileName - name of the file to write to
%  [struct] Data - struct containing the data
%    [cell] Dimensions - list of dimensions we want to declare
% 
% "Dimensions" must EXIST as an input, but can be the empty cell (i.e. "{}"). 
% Any unspecified dimensions will be automatically named and generated.
%
%optional inputs:
%
%  metadata:
%      [char] title      - file title
%      [char] long_title - file long title
%      [char] CreatedBy  - creator name
%    [struct] MetaData   - each field should be a string, which will be added to the metadata
%
%  dimension handling
%     [cell] DimUnits     - cell array of strings containing units of each specified dimension, in same order as Dimensions
%     [cell] DimNames     - cell array of strings containing names of each specified dimension, in same order as Dimensions
%
%  variable options
%     [struct] VarProperties - struct of structs, where each member:
%                              a. is named the same as the variable of interest
%                              b. contains fields which are strings, with possible names
%                                 Name - variable name to use. IF THIS IS THE SAME AS ANOTHER VARIABLE, BEHAVIOUR WILL BE UNPREDICTABLE. DO NOT DO THIS.
%                                 FullName - variable description
%                                 Units - units
%                                 Fill - fill value
%                              All of these fields are optional
%  
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/02/13
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assume failure
Error = 1;

%what variables types can this function save?
ValidList = {'char','string','logical','int8','int16','int32','single','double','logical'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create input parser and testing functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;


%inputs - required
%%%%%%%%%%%%%%%%%%%

addRequired(p,'FileName',  @ischar);  
addRequired(p,'Data',      @isstruct); 
addRequired(p,'Dimensions',@iscell); 

%inputs - optional
%%%%%%%%%%%%%%%%%%%

%file-level metadata
addParameter(p,'title',  '',  @ischar);  
addParameter(p,'long_title',  '',@ischar);  
addParameter(p,'CreatedBy',  '', @ischar);  
addParameter(p,'MetaData',  struct(), @isstruct);  

%dimension handling
addParameter(p,'DimUnits',  cell(numel(Dimensions),1), @iscell);  
addParameter(p,'DimNames',  cell(numel(Dimensions),1), @iscell);  

%variable properties
addParameter(p,'VarProperties',struct(), @isstruct);  

%parse inputs
%%%%%%%%%%%%%%%%%%%%%%%

parse(p,FileName,Data,Dimensions,varargin{:});
Inputs = p.Results; 
clear p varargin


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% specifify file-level metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create generic metadata
MetaData.title = Inputs.title;
MetaData.long_title = Inputs.long_title;
MetaData.CreatedBy = Inputs.CreatedBy;

%create additional metadata
Fields = fieldnames(Inputs.MetaData);
for iField=1:1:numel(Fields)
  Name = Fields{iField};
  Meta = Inputs.MetaData.(Fields{iField});
  if ~strcmp(class(Meta),'char')
    disp(['-->Metadata field ''',Name,''' is not a character string, skipping'])
  else
    MetaData.(Name) = Meta;
  end
end; clear iField Name Meta Fields


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check the type of each variable, and discard those we can't handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


VarList = fieldnames(Data);
for iVar=1:1:numel(VarList)

  %check the type
  Type = class(Data.(VarList{iVar}));
  
  %reject if not acceptable
  if ~sum(any(strcmp(Type,ValidList)))
    disp(['--> Variable ''',VarList{iVar},''' is not of a valid type and will not be saved. Valid Types are:'])
    for iType=1:1:numel(ValidList); disp(['----> ',ValidList{iType}]); end
    Data = rmfield(Data,VarList{iVar});
    continue
  end

  %if complex, split into a real and imag var
  if ~isreal(Data.(VarList{iVar})) & ~isstring(Data.(VarList{iVar}));
    Data.([VarList{iVar},'_real']) = real(Data.(VarList{iVar}));
    Data.([VarList{iVar},'_imag']) = imag(Data.(VarList{iVar}));
    Data = rmfield(Data,VarList{iVar});
    disp(['--> Complex variable ''',VarList{iVar},''' has been split into real and imaginary parts for storage'])
  end

  %if logical, convert to single
  if strcmp(Type,'logical'); Data.(VarList{iVar}) = single(Data.(VarList{iVar})); end;

end; clear iVar VarList ValidList iType Type


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% specify dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if we fed any explicit dims to the programme, and handle them if we did
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(Dimensions) ~= 0;

  %create structure for Dimensions
  DimData    = struct();
  DimL       = [];

  %populate input dimensions in the order specified
  for iDim=1:1:numel(Dimensions)

    DimData = cjw_nc_prepop_dim(DimData,Dimensions{iDim},Data.(Dimensions{iDim}));
    DimData(end).FullName = Inputs.DimNames{iDim};
    DimData(end).Units    = Inputs.DimUnits{iDim};
    DimData(end).Fill     = -9999;
    DimData(end).Axis     = 1:1:length(Data.(Dimensions{iDim}));
    DimData(end).Type     = 'double';


    %store the length, for identifying below
    DimL(iDim)      = length(Data.(Dimensions{iDim}));

  end; clear iDim

else

  %blank dims array as a baseline
  DimL = [];
  DimData = struct();
end

%now, create any others we might need, with arbitrary names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get a list of variables
VarList = fieldnames(Data);

%remove the dimensions we selected above from this list
for iDim=1:1:numel(Dimensions)
  cf = strcmp(Dimensions{iDim},VarList);
  if sum(cf) > 0; VarList(cf == 1) = []; end
end; clear iDim cf

%find the size of all the variables
VarSizes = NaN(numel(VarList),9999); %this will not work if any of the variables have 10000 or more dimensions. Sorry.
for iVar=1:1:numel(VarList)
  sz = size(Data.(VarList{iVar}));
  try;  VarSizes(iVar,1:numel(sz)) = sz;
  catch; disp('Routine cannot handle arrays with >9999 dimensions. Stopping.');return;
  end
end; clear iVar sz
VarSizes = VarSizes(:,1:max(sum(~isnan(VarSizes),2)));

%for any that we don't have yet, create an axis array
DimLengths = unique([DimL';unique(VarSizes(:))],'stable'); 
DimLengths = DimLengths(~isnan(DimLengths));
for iDL=numel(DimL)+1:1:numel(DimLengths)

  if ~any(DimL == DimLengths(iDL))
    %create a new dimension
    if numel(fieldnames(DimData)) == 0;ID = 1; else; ID = numel(DimData)+1; end

    DimData = cjw_nc_prepop_dim(DimData,['dim_',num2str(ID)],1:1:DimLengths(iDL));
    DimData(end).FullName = 'no name,  auto-generated';
    DimData(end).Units    = 'no units, auto-generated'; 
    DimData(end).Fill     = -9999;
    DimData(end).Axis     = 1:1:DimLengths(iDL);
    DimData(end).Type     = 'double';    
  end

end; clear iDL ID

clear Dimensions DimL VarSizes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%now, loop over the variables and prepare to stick them in the netCDF file
OutData = struct();

for iVar=1:1:numel(VarList)

  %get variable, variable name, class and size
  Var = Data.(VarList{iVar});
  Name = VarList{iVar};
  sz = size(Var);
  Class = class(Var);

  %define some default parameters, to be overwritten if needed
  FullName = '';
  Units = '';
  if any(strcmp(Class,{'single','double'})); Fill = -9999; 
  elseif isnumeric(Var);                     Fill = 0;
  else                                       Fill = '';
  end

  %get specified properties, if they exist
  if isfield(Inputs,'VarProperties');
    if isfield(Inputs.VarProperties,Name)
      if isfield(Inputs.VarProperties.(Name),'FullName'); FullName = Inputs.VarProperties.(Name).FullName; end
      if isfield(Inputs.VarProperties.(Name),   'Units'); Units    = Inputs.VarProperties.(Name).Units;    end
      if isfield(Inputs.VarProperties.(Name),    'Fill'); Fill     = Inputs.VarProperties.(Name).Fill;     end    
      if isfield(Inputs.VarProperties.(Name),    'Name'); Name     = Inputs.VarProperties.(Name).Name;     end         
    end
  end

  %using the size, identify the axes
  Axes = sz.*NaN;
  for iAxis=1:1:numel(sz)
    Axes(iAxis) = find(DimLengths == sz(iAxis));
  end

  %hence, create the data entry
  OutData = cjw_nc_prepop_data(OutData,Name,Class);
  OutData.(Name).Dims     = Axes;
  OutData.(Name).FullName = FullName;
  OutData.(Name).Units    = Units;
  OutData.(Name).Data     = Var;
  OutData.(Name).Fill     = Fill;

  clear Fill Units FillName

end; clear iVar Var Name sz Class Axes iAxis


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%create file
cjw_nc_create(FileName,MetaData,1);

%fire off dimensions to the netCDF file
cjw_nc_makedims(FileName,DimData);

%write the data
Success = cjw_nc_writedata(FileName,OutData);
if Success == 1; disp(['--> Data saved as ',FileName]); Error = 0;
else             disp(['--> Error writing file ',FileName,'; data not saved']); Error = 0;
end


%done!
Error = 0;



end
