function [Topography,Map,CloudMap,ShapeData] = topo_etc(LonBox,LatBox,GreyScale,Modis,NatEarth,HRNatEarth)


%logic has some trouble with ending at exactly 180 or 90 degrees:
if LonBox(1) == -180; LonBox(1) = -179.999; end
if LatBox(1) ==  -90; LatBox(1) =  -89.999; end
if LonBox(2) ==  180; LonBox(2) =  179.999; end
if LatBox(2) ==   90; LatBox(2) =   89.999; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get topography, coastline etc layers for a given region's 3D plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%replace map with ...
if exist('GreyScale')  ==0; GreyScale = 0;   end;  %greyscale map
if exist('Modis'    )  ==0; Modis = 0;       end;  %modis map
if exist('NatEarth')   ==0; NatEarth = 0;    end;  %natural earth map
if exist('HRNatEarth') ==0; HRNatEarth = 0;  end;  %high-res natural earth map
%earliest-specified will be chosen, so set any not wanted before choice to zero

%define region to work in
minlon = LonBox(1); maxlon = LonBox(2);
minlat = LatBox(1); maxlat = LatBox(2);

%load topo
%%%%%%%%%%%%%%%%%%%%%%%%%%%

tp = load([LocalDataDir,'/topography/easy_tenth_degree_topography/easy_topo.mat']);
% regionally subset
tp = tp.topo;
% tp.elev = tp.topo; %unsmoothed unfiltered topography
tp.lons = -180:0.1:180;
tp.lats = -90:0.1:90;
% tp = rmfield(tp,'topo');
InLonRange = find(tp.lons > LonBox(1) & tp.lons <= LonBox(2));
InLatRange = find(tp.lats > LatBox(1) & tp.lats <= LatBox(2));
tp.elev = tp.elev(InLatRange,InLonRange);
tp.lons = tp.lons(InLonRange);
tp.lats = tp.lats(InLatRange);
[tp.lons,tp.lats] = meshgrid(tp.lons,tp.lats);
clear Lon Lat InLOnRange InLatRange
%smooth away the biggest spikes
tp.elev = smoothn(tp.elev,[1,1].*1);
tp.elev(tp.elev <= -1) = -5; %shallower sea

tp.elev = tp.elev ./ 1000; % CONVERT TO KM!!!

Topography = tp; clear tp;


%% LOAD MAP ===============================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GreyScale == 1;
  Map.Map = flipud(imread([LocalDataDir,'/topography/imagery/greyscale.png']));  
  Map.Map = repmat(Map.Map,1,1,3);
elseif Modis == 1;
  Map.Map = flipud(imread([LocalDataDir,'/topography/imagery/MODIS_Map.jpg']));
elseif NatEarth == 1;
  Map.Map = flipud(imread([LocalDataDir,'/topography/ne/rasterI/NE1_50M_SR_W.tif']));
elseif HRNatEarth == 1;
  Map.Map = flipud(imread([LocalDataDir,'/topography/ne/rasterI/HYP_HR_SR_OB_DR.tif']));
else
  Map.Map = flipud(imread([LocalDataDir,'/topography/imagery/land_ocean_ice_8192.png']));
end

%work out which part we need
MapLonScale = -180:360/size(Map.Map,2):180;
MapLatScale = -90:180/size(Map.Map,1):90;
InLonRange = find(MapLonScale > LonBox(1) & MapLonScale <= LonBox(2));
InLatRange = find(MapLatScale > LatBox(1) & MapLatScale <= LatBox(2));
Map.Map = Map.Map(InLatRange,InLonRange,:);
Map.LatScale = MapLatScale(InLatRange);
Map.LonScale = MapLonScale(InLonRange);
clear InLonRange InLatRange MapLonScale MapLatScale

%% LOAD CLOUDS ============================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%
CloudMap.Map = flipud(imread([LocalDataDir,'/topography/imagery/MODIS_map.jpg']));
%work out which part we need
MapLonScale = -180:360/size(CloudMap.Map,2):180;
MapLatScale = -90:180/size(CloudMap.Map,1):90;
InLonRange = find(MapLonScale > LonBox(1) & MapLonScale <= LonBox(2));
InLatRange = find(MapLatScale > LatBox(1) & MapLatScale <= LatBox(2));
CloudMap.Map = CloudMap.Map(InLatRange,InLonRange,:);
CloudMap.LatScale = MapLatScale(InLatRange);
CloudMap.LonScale = MapLonScale(InLonRange);
clear InLonRange InLatRange MapLonScale MapLatScale

%% LOAD COASTLINE =========================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ShapeData = shaperead([LocalDataDir,'/topography/ne/10mcoasts/ne_10m_coastline.shp'],'BoundingBox',[minlon minlat;maxlon maxlat],'UseGeoCoords', true);


end

