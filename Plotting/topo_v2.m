function [Topo,Coasts,Image] = topo_v2(LonBox,LatBox,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load topography, coastline and image data for use in generating plots
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/09/16
%
%inputs:
%  required:
%    lonbox: two-element array [min lon, max lon], in range -360 to +360
%    latbox: two-element array [min lat, max lat], in range -90 to +90
%  optional:
%    'Image': image to load, one of 'GreyScale', 'Modis','NatEarth','HRNatEarth', 'land_ocean_ice', 'pale','land_ocean_ice_cloud'
%
%outputs:
%  Topo:   struct containing grid of lons, lat and elevation values for the chosen region
%  Coasts: shapefile data of coastlines in the selected region
%  Image:  Image data of the surface in the selected region
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%set the paths we will use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DataDir  = [LocalDataDir,'/topography/'];
Settings.ImageDir = [LocalDataDir,'/topography/'];
Settings.ShapeDir = [LocalDataDir,'/topography/'];

%load topography data, and duplicate into the full range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load
Topo = load([Settings.DataDir,'/easy_tenth_degree_topography/easy_topo.mat']);
Topo = Topo.topo;

%modifications
Topo.elev(Topo.elev < -5) = -5; %oceans are flat
Topo.elev = Topo.elev./1000; %scale to km


%duplicate out so we can use -180 to +180 or 0-360 scales
Topo.lons = Topo.lons(:,1:end-1);
Topo.lats = Topo.lats(:,1:end-1);
Topo.elev = Topo.elev(:,1:end-1);
Topo.elev = [Topo.elev,    Topo.elev,Topo.elev    ];
Topo.lats = [Topo.lats,    Topo.lats,Topo.lats    ];
Topo.lons = [Topo.lons-360,Topo.lons,Topo.lons+360];

%load imagery data, and duplicate into full range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = find(contains(varargin,'Image'));

if numel(idx) > 0;
  
  %find the path to the image
  ImageToLoad = varargin(idx+1);
  Skip = 0;
  switch ImageToLoad{1}
    case 'GreyScale';            Path = 'imagery/greyscale.png';
    case 'Modis';                Path = 'imagery/MODIS_Map.jpg';
    case 'NatEarth';             Path = 'ne/rasterI/NE1_50M_SR_W.tif';
    case 'HRNatEarth';           Path = 'ne/rasterI/HYP_HR_SR_OB_DR.tif';
    case 'land_ocean_ice';       Path = 'imagery/land_ocean_ice_8192.png';
    case 'land_ocean_ice_cloud'; Path = 'imagery/land_ocean_ice_cloud_8192.png';
    case 'faded';                Path = 'imagery/faded.jpg';
    case 'pale';                 Path = 'imagery/pale.png';
      
    otherwise;
      disp('Image specified not location, skipping load')
      Skip = 1;
  end
  
  Path = [Settings.ImageDir,Path];
  
  %check the file exists
  if ~exist(Path,'file')
    disp('Image file not in specified location, skipping load');
    Skip = 1;
  end
  
  %load the file
  Image.Map = imread(Path);
  
  %greyscale needs duplicating out to have three colours, even though they're all the same
  if strcmp(ImageToLoad{1},'GreyScale'); Image.Map = repmat(Image.Map,1,1,3); end
    
  %create corresponding lat and lon arrays (assumes Mercator projection)
  Image.Lon = linspace(-180,180,1+size(Image.Map,2)); Image.Lon = Image.Lon(1:end-1);
  Image.Lat = linspace(  90,-90,1+size(Image.Map,1)); Image.Lat = Image.Lat(1:end-1);

  %duplicate out
  Image.Map = [Image.Map,Image.Map,Image.Map];
  Image.Lon = [Image.Lon-360,Image.Lon,Image.Lon+360]; 

end



%load coastline data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Coasts = shaperead([Settings.ShapeDir,'/ne/10mcoasts/ne_10m_coastline.shp'], ...
                      'BoundingBox',[[-180,180];LatBox]',                           ...
                      'UseGeoCoords', true);

%now, duplicate out outside the range (if needed)
if max(LonBox) > 180 | min(LonBox) < -180;
  for iC=1:1:numel(Coasts)
    
    C = Coasts(iC);
    C.Lon = C.Lon + 360;
    C.BoundingBox(:,1) = C.BoundingBox(:,1)+360;
    Coasts(end+1) = C;
    
    C = Coasts(iC);
    C.Lon = C.Lon - 360;
    C.BoundingBox(:,1) = C.BoundingBox(:,1)+360;
    Coasts(end+1) = C;  
    
  end
end
                   

%now, subset all the datasets to the chosen region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%topography
x = inrange(Topo.lons(1,:),LonBox);
y = inrange(Topo.lats(:,1),LatBox);
Topo.lons = Topo.lons(y,:);  Topo.lons = Topo.lons(:,x);
Topo.lats = Topo.lats(y,:);  Topo.lats = Topo.lats(:,x);
Topo.elev = Topo.elev(y,:);  Topo.elev = Topo.elev(:,x);


%imagery
if exist('Image','var');
  
  %being exactly at the limits causes a problem - this is the last step, so
  %adjust if this happens
  if LatBox(2) ==  90; LatBox(2) =  89.99999; end
  if LatBox(1) == -90; LatBox(1) = -89.99999; end

  %ok, now subset the map
  x = inrange(Image.Lon,LonBox);
  y = inrange(Image.Lat,LatBox);
  Image.Lon = Image.Lon(x);
  Image.Lat = Image.Lat(y); Image.Lat = fliplr(Image.Lat);
  Image.Map = Image.Map(y,:,:);
  Image.Map = Image.Map(:,x,:);
  Image.Map = flipud(Image.Map);

else
  Image = 'no image loaded';
end

%coasts
Coasts2 = Coasts(1:2); %this will be ignored, is just to initialise the structure correctly
k = 0;
for iC=1:1:numel(Coasts)
  if max(Coasts(iC).Lon) < min(LonBox); continue; end
  if min(Coasts(iC).Lon) > max(LonBox); continue; end
  k = k+1;
  Coasts2(k) = Coasts(iC);
end
Coasts = Coasts2;



