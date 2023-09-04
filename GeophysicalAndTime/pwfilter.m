function [VarOut,PWStore] = pwfilter(NPWs,MinPC,Lon,Lat,Var,LonGrid,LatGrid,Alt,AltGrid)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simple planetary wave filter for application to scattered data as f(lon,lat,alt)
%altitude is optional - omit last two arguments to apply to 2D data only
%
%the first mode will always be the zonal mean, then the other modes
%will start from number 2 (i.e. modes 1,2,3 are in positions 2,3,4, etc)
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2022/09/23
%
%inputs (not in order - check above):
% NPWs    - number of planetary wave modes to fit. Zonal mean will always be included.
% MinPC   - minimum percentage of filled bins in a longitude band to compute PWs - will set to NaN otherwise
% Lon     - list of longitudes  at data points
% Lat     - list of latitude    at data points
% Alt     - list of altitudes   at data points (optional)
% Var     - list of data values at data points
% LonGrid - longitude grid to compute results on (PW outputs will also be on this grid)
% LatGrid - latitude  grid compute results on (PW outputs will also be on this grid)
% AltGrid - altitude  grid compute results on (PW outputs will also be on this grid)
%
%Lon, Lat, Alt and Var must all be the same size in the same order, and will be flattened
%LonGrid,LatGrid and AltGrid should be vectors, and will be meshgridded
%
%outputs:
% VarOut - the total PW contribution (inc zonal mean) at each point in the raw data
% PWStore - planetary wave estimated amplitude at each point on the output grid, for each mode (lon x lat x alt x PW)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check input validity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if we didn't specify altitude values and ranges, create dummy values
if nargin ~= 9; Alt = ones(size(Lon)); AltGrid = 1;end

%check input types, ranges and consistency.
p = inputParser;
addRequired(p,'NPWs',    @(x) validateattributes(x,{'numeric'},{'positive','integer','scalar'}));
addRequired(p,'MinPC',   @(x) validateattributes(x,{'numeric'},{'nonnegative','<=',100,'scalar'}));
addRequired(p,'Lon',     @isnumeric);
addRequired(p,'Lat',     @(x) validateattributes(x,{'numeric'},{'size',size(Lon)}));
addRequired(p,'Var',     @(x) validateattributes(x,{'numeric'},{'size',size(Lon)}));
addRequired(p,'LonGrid', @(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',360,'vector'}));
addRequired(p,'LatGrid', @(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90,'vector'}));
addRequired(p,'Alt',     @(x) validateattributes(x,{'numeric'},{'nonnegative','size',size(Lon)}));
addRequired(p,'AltGrid', @(x) validateattributes(x,{'numeric'},{'nonnegative','vector'}));
parse(p,NPWs,MinPC,Lon,Lat,Var,LonGrid,LatGrid,Alt,AltGrid);
clear p



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flatten input fields, retaining the input size
InSize = size(Lon);
Lon = Lon(:); Lat = Lat(:); Alt = Alt(:); Var = Var(:);

%create output meshgrids
[LonGrid,LatGrid,AltGrid] = meshgrid(LonGrid,LatGrid,AltGrid);

%bin the input data onto the grid
VarGrid = bin2matN(3,Lon,Lat,Alt,Var,LonGrid,LatGrid,AltGrid,'@nanmean');

%reshape to remove nested loop
sz = size(VarGrid); if numel(sz) < 3; sz = [sz,1,1]; end
VarGrid = reshape(permute(VarGrid,[2,1,3]),sz(2),sz(1)*sz(3));

%drop rows that don't meet MinFrac
Filled = sum(~isnan(VarGrid),1);
VarGrid(:,Filled./sz(2) <= MinPC./100) = NaN;

clear sz Filled Bad MinPC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PW analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract longitude axis
LonAxis = LonGrid(1,:,1);

%create storage arrays
PWStore = NaN([size(VarGrid),NPWs+1]); %map of PWs
VGrid   = zeros(size(VarGrid));        %sum of planetary waves at each point

%compute planetary wave value at each data point
for iPW=0:1:NPWs
  for iLine =1:1:size(VarGrid,2)

    if iPW==0;
      %take and repmat zonal mean
      yfit = repmat(nanmean(VarGrid(:,iLine)),[size(VarGrid,1),1]);
    else
      %fit PW
      yfit = nph_sinefit(LonAxis,VarGrid(:,iLine),360./iPW);
    end

    %remove from the data so we don't do so again
    VarGrid(:,iLine) = VarGrid(:,iLine) - yfit;

    %store results
    VGrid(:,iLine) = VGrid(:,iLine)+yfit;
    PWStore(:,iLine,iPW+1) = yfit;

    clear yfit F
  end; clear iLine
end; clear iPW
clear VarGrid LonAxis

%reshape back to original grid shape
sz = size(LatGrid); if numel(sz) < 3; sz = [sz,1,1]; end
VGrid   = permute(reshape(  VGrid,sz([2,1,3])),[2,1,3]);
PWStore = permute(reshape(PWStore,[sz([2,1,3]),NPWs+1]),[2,1,3,4]);
clear sz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% remove from raw profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create an interpolant object of the PW sum
%this requires permuting each input from meshgrid to ndgrid
%it also requires splitting out by the 2D and 3D case. Bah, so close.

nd = ndims(LonGrid);
if nd == 3;
  %3D case

  I = griddedInterpolant(permute(LonGrid,[2,1,3]),...
                         permute(LatGrid,[2,1,3]),...
                         permute(AltGrid,[2,1,3]),...
                         permute(VGrid,  [2,1,3]));


  %interpolate to the profile locations, difference from the raw data, 
  %and reshape to the original field shape
  VarOut = reshape(Var-I(Lon,Lat,Alt),InSize);

elseif nd == 2;

  I = griddedInterpolant(permute(LonGrid,[2,1]),...
                         permute(LatGrid,[2,1]),...
                         permute(VGrid,  [2,1]));


  %interpolate to the profile locations, difference from the raw data, 
  %and reshape to the original field shape
  VarOut = reshape(Var-I(Lon,Lat),InSize);

else
  disp('Error - invalid number of output dimensions')
end


clear I InSize Lon Lat Alt Var AltGrid LonGrid LatGrid VGrid nd

