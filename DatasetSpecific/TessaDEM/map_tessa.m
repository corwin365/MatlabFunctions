function [Alt,LonScale,LatScale,TileScript] = map_tessa(LonRange,LatRange,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to generate a map of a region using TESSA-DEM data
%
%requires TESSA tiles to be present on local system, but can 
%generate an SCP command to download them from eepc-0184 repository
%as an optional output
%
%Required inputs:
%
%  LonRange - range of longitudes required, in -180 to +180 format 
%  LatRange - range of latitudes required, in -90 to +90 format
%
%Optional Inputs:
%
%  ETFill - default true - fill polar regions without coverage with 0.1 degree easyTopo data (requires topo_V2.m)
%  DataDir - default [LocalDataDir,'/topography/tessaDEM/raw/'] - path to TESSA tiles
%  TileScript - default false - return an SCP command to get the ftiles rather than the tiles themselves
%  Resolution - default [1/112e3*30,1/112e3*30] - output resolution in [lon, lat] in degrees
%
%Outputs:
%  Alt - grid of altitudes on requested resolution
%  LonScale - longitude axis for this grid
%  LatScale - latitude axis for this grid
%  TileScript - SCP script to acquire desired tiles
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/04/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create parser object
%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

%validation of inputs
%%%%%%%%%%%%%%%%%%%%%%

%lat/lon, required
addRequired(p,'LonRange',@(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',180}))
addRequired(p,'LatRange',@(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90}))

%optional flags
addParameter(p,'ETFill',   true,@islogical); %use 0.1 easyTopo topography to fill polar regions
addParameter(p,'TileScript',false,@islogical); %if this is set to true, function just returns a lsit of required tiles for this map

%optional settings
addParameter(p,'DataDir',[LocalDataDir,'/topography/tessaDEM/raw/'],@ischar); %path to data
addParameter(p,'Resolution',[1/112e3*30,1/112e3*30],@isnumeric); %defaults to 30m resolution, i.e. dataset limit

%done - parse inputs
parse(p,LonRange,LatRange,varargin{:})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%work out the unique tiles we need
LonTiles = floor(min(LonRange)):1:ceil(max(LonRange));
LatTiles = floor(min(LatRange)):1:ceil(max(LatRange));
[LonTiles,LatTiles] = meshgrid(LonTiles,LatTiles);
TileList(1,:) = LonTiles(:);
TileList(2,:) = LatTiles(:);

%work out the final product grid
LonScale = min(LonRange):p.Results.Resolution(1):max(LonRange);
LatScale = min(LatRange):p.Results.Resolution(2):max(LatRange);
Alt = NaN(numel(LonScale),numel(LatScale));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if we just want a dwnload list of files, generate that
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.Results.TileScript == true;

  %generate the script
  TileScript = "scp -OT USERNAME@eepc-0184.bath.ac.uk:""" ;
  DataDir = "/data1/topography/tessaDEM/raw/";
  
  for iTile=1:1:size(TileList,2)
    FileName = num2str(TileList(2,iTile))+"_"+num2str(TileList(1,iTile));
    TileScript = TileScript+DataDir+FileName+" ";
  end

  TileScript = TileScript+""" ./";

  disp('===================================');
  disp('Variable ''TileScript'' contains an');
  disp('SCP command to download the required');
  disp('tiles. Please modify to include your');
  disp('username, then run in the target');
  disp('directory. You will need to set the');
  disp('''DataDir'' optional variable to this');
  disp('directory to access the downloaded data');
  disp('===================================');

  return; 
else
  TileScript = "";
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for each tile, load and put on the final grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iTile=1:1:numel(LonTiles)

  %load tile
  [Z,Lon,Lat] = load_tessaDEM_tile(LonTiles(iTile),LatTiles(iTile),'DataDir',p.Results.DataDir,'ETFill',p.Results.ETFill);

  %create an interpolant from the data
  I = griddedInterpolant({Lon,Lat},double(Z));

  %create an output grid, and interpolate
  idx = inrange(LonScale,minmax(Lon));
  jdx = inrange(LatScale,minmax(Lat));
  [idx,jdx] = ndgrid(idx,jdx);
  xi = LonScale(idx); yi = LatScale(jdx);
  zz = I(xi,yi);

  %store
  for iPoint=1:1:numel(idx); Alt(idx(iPoint),jdx(iPoint)) = zz(iPoint); end




end