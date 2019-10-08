function [ShapeData] = coasts_only(LonBox,LatBox)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get coastline layer for a given region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define region to work in
minlon = LonBox(1); maxlon = LonBox(2);
minlat = LatBox(1); maxlat = LatBox(2);


%% LOAD COASTLINE =========================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ShapeData = shaperead([LocalDataDir,'/topography/ne/10mcoasts/ne_10m_coastline.shp'],'BoundingBox',[minlon minlat;maxlon maxlat],'UseGeoCoords', true);


end

