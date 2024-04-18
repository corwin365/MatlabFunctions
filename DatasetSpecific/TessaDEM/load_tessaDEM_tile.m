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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %following option removed and moved upstream
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % elseif p.Results.ETFill == true
% % 
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %% no tile: use easyTopo to fill the gap 
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% %   if Lat > 84 | Lat < -80;
% %     %we are outside TESSA coverage - use easyTopo 0.1 degree
% %     ET = topo_v2(Lon+[0,1],Lat+[0,1]);
% %     Alt = fliplr(ET.elev);
% %     LonScale = squeeze(ET.lons(1,:));
% %     LatScale = squeeze(ET.lats(:,1));
% %     clear ET
% %   else
% %     %open ocean - all-zero
% %     Alt = [0,0;0,0];
% %     LonScale = Lon+[0,1];
% %     LatScale = Lat+[0,1];
% % 
% %   end
% % 
% %   Error = 1;
% %   return

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


end
