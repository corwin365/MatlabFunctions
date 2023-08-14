function [Az,El,ID] = find_star_azel(Lat,Lon,Time,MagLimit,FilePath)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%find starfield at a given (lat,lon) as a function of (az,el)
%%cut off stars below a specified magnitude
%%Corwin Wright, c.wrigh@bath.ac.uk
%%06/AUG/2017
%%
%%inputs:
%%Lat - latitude of observation site (deg)
%%Lon - longitude of observation site (deg)
%%Time - time of observation (UTC)
%%MagLimit - magnitude limit of stars to include
%%FilePath - path to data csv (optional)
%%
%%outputs:
%%Az - azimuth of stars (deg)
%%El - elevation of stars (deg)
%%ID - ID number of star in catalogue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% % % % % % %%test settings
% % % % % % clear all
% % % % % % FilePath = 'yale_stars.csv';
% % % % % % Lat = 51.61;
% % % % % % Lon = 1.25;
% % % % % % Time = datenum(2017,8,6,14,47,0);
% % % % % % MagLimit = 3.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get star database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin <5; FilePath = [LocalDataDir,'/Miscellany/yale_stars.csv']; end
%load the star database
Stars = csvread(FilePath,2,0);
clear FilePath

Mag = Stars(:,7);
RA  = Stars(:,8);
Dec = Stars(:,9);
clear Stars


%discard dim stars (brightest = lowest mag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Good = find(Mag <= MagLimit);
if numel(Good) == 0; Az = NaN; El = NaN; Good = NaN; return; end %no stars found
RA  = RA(Good);
Dec = Dec(Good);
ID = Good;
clear Mag MagLimit


%convert what's left to az and el from our location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Az,El] = RaDec2AzEl(RA,Dec,Lat,Lon,datestr(Time,'yyyy/mm/dd HH:MM:SS'));

%done!