function [Alt,LonPoints,LatPoints,TileScript] = map_tessa(LonPoints,LatPoints,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to generate a map of a region using TESSA-DEM data
%
%requires TESSA tiles to be present on local system, but can 
%generate an SCP command to download them from eepc-0184 repository
%as an optional output
%
%Required inputs:
%
%  LonRange - longitudes required, in -180 to +180 format 
%  LatRange - latitudes required, in -90 to +90 format
%
%Optional Inputs:
%
%  ETFill - default true - fill polar regions without coverage with 0.1 degree easyTopo data (requires topo_V2.m)
%  DataDir - default [LocalDataDir,'/topography/tessaDEM/raw/'] - path to TESSA tiles
%  TileScript - default false - return an SCP command to get the ftiles rather than the tiles themselves
%
%Outputs:
%  Alt - grid of altitudes on requested resolution
%  LonPoints - longitudes of the output
%  LatPoints - latitudes  of the output
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
addRequired(p,'LonPoints',@(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',180}))
addRequired(p,'LatPoints',@(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90,'size',size(LonPoints)}))

%optional flags
addParameter(p,'ETFill',   true,@islogical); %use 0.1 easyTopo topography to fill polar regions
addParameter(p,'ETPath',   [LocalDataDir,'/topography/easy_tenth_degree_topography/easy_topo.mat'],@ischar); %path to easyTopo
addParameter(p,'TileScript',false,@islogical); %if this is set to true, function just returns a lsit of required tiles for this map

%optional settings
addParameter(p,'DataDir',[LocalDataDir,'/topography/tessa/'],@ischar); %path to data
addParameter(p,'Resolution',[1/112e3*30,1/112e3*30],@isnumeric); %defaults to 30m resolution, i.e. dataset limit

%done - parse inputs
parse(p,LonPoints,LatPoints,varargin{:})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flatten the lats and lons. We'll reshape them back at the end



%work out the unique tiles we need
TileList(1,:) = floor(LonPoints(:));
TileList(2,:) = floor(LatPoints(:));
TileList = unique(TileList','rows')';

%work out the final product grid
Alt = NaN(size(LonPoints));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if we just want a dwnload list of files, generate that
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.Results.TileScript == true;

  %generate the script
  TileScript = "scp -T USERNAME@eepc-0184.bath.ac.uk:""" ;
  DataDir = "/data1/topography/tessa/";
  
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


for iTile=1:1:size(TileList,2)

  %load tile
  [Z,Lon,Lat,Error] = load_tessaDEM_tile(TileList(1,iTile),TileList(2,iTile),'DataDir',p.Results.DataDir);

  %make sure we actually got a tile 
  if Error ~= 0; continue; end

  %create an interpolant from the data
  %disable extrapolation disabled so it will only do values for this granule
  I = griddedInterpolant({Lon,Lat},double(Z),'linear','none');
  clear Z Lon Lat

  %interpolate onto the full set of points
  warning off %gives a warning about data format which is optimisation only and not relevant
  Vals = I(LonPoints,LatPoints);
  warning on

  %store the ones that work
  Alt(~isnan(Vals)) = Vals(~isnan(Vals));


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fill gaps with easyTopo lower-res data?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.Results.ETFill == true

  %find the empty regions
  Empty = find(isnan(Alt) | Alt == 0);
  if numel(Empty) == 0; return; end %we already have all the points!

  %what points do we need?
  ETLons = LonPoints(Empty);
  ETLats = LatPoints(Empty);

  %load the easyTopo data
  ET = load(p.Results.ETPath); ET = ET.topo;
  
  %generate interpolant, and fill those gaps
  ET.lons = ET.lons'; ET.lats = ET.lats'; ET.elev = ET.elev';
  I = griddedInterpolant(ET.lons,ET.lats,ET.elev);
  Alt(Empty) = I(ETLons,ETLats);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MinMax = minmax(Array)
MinMax = [nanmin(Array(:)),nanmax(Array(:))];
return

function InRange = inrange(Array,MinMax,NoEnds)
InRange = find(Array >  min(MinMax) & Array <  max(MinMax));
return



function [Alt,LonScale,LatScale,Error] = load_tessaDEM_tile(Lon,Lat,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to load and parse a tile from the TESSA DEM
%see further details below
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/04/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%data obtained from https://tessadem.com/, see there for rules on usage
%
%Format:
% % Raw binary data
% % Total size: 405 GB
% % Data splitting: 22,551 files of 1° latitude x 1° longitude.
% % File name: [latitude]_[longitude] of the lower left corner.
% % File rows: 3,600
% % File columns: 3,600 below 50° latitude; 2,400 between 50° and 60° latitude; 1,800 between 60° and 70° latitude; 1,200 between 70° and 80° latitude; 720 above 80° latitude.
% % Data type: 16-bit signed short integer (little endian byte order).
% % Data unit: meter
% % No data value: -9999
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create parser object
%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

%validation of inputs
%%%%%%%%%%%%%%%%%%%%%%

%lat/lon, required
addRequired(p,'Lon',@(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',180}))
addRequired(p,'Lat',@(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90}))

%other optional
addParameter(p,'DataDir',[LocalDataDir,'/topography/tessa/'],@ischar);
% % % addParameter(p,'ETFill',true,@islogical);  %option removed and moved upstream

%done - parse inputs
parse(p,Lon,Lat,varargin{:})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% work out spacing of points from latitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if     abs(Lat) < 50; Cols = 3600;
elseif abs(Lat) < 60; Cols = 2400;
elseif abs(Lat) < 70; Cols = 1800;
elseif abs(Lat) < 80; Cols = 1200;
else                  Cols = 720;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify file needed for this tile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LatString = num2str(floor(Lat));
LonString = num2str(floor(Lon));
FileName = [p.Results.DataDir,'/',LatString,'_',LonString];
clear LatString LonString 

if exist(FileName,'file')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load tile
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fid = fopen(FileName);
  Alt = fliplr(fread(fid,[Cols,3600], '*int16', 0, 'ieee-le'));
  fclose(fid);

  LatScale = linspace(Lat,Lat+1,3600);
  LonScale = linspace(Lon,Lon+1,Cols);

  Error = 0;
  return

else

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% no tile: do nothing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Alt      = [];
  LonScale = [];
  LatScale = [];
  Error = 2;
  return

end


return