function [Moon,Stars] = moon_and_stars(Lat,Lon,Time,MagCutOff,Alt)


% % % %testing parameters
% % % clear all
% % % Lat = 42.6;
% % % Lon = -71.5;
% % % Alt = 0;
% % % Time = [0,1e-3]+736545.187152778+145;
% % % MagCutOff = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the relative location of the moon and major stars to a given place
%
%Corwin Wright, c.wright@bath.ac.uk, 09/AUG/2017
%
%inputs:
% Lat,Lon,Alt - lat,lon,altitude [m] of site (alt optional, 0 if not specified)
% Time        - time (UTC) of interest
% MagCutOff   - maximum magnitude of stars(optional, 9e99 if not specified)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5; Alt = 0; end
if nargin < 4; MagCutOff  = 9e99; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% moon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MoonPhase = load([LocalDataDir,'/miscellany/Moon/lunar_phases.mat']);
MoonPhase = MoonPhase.MoonPhase;


Moon = NaN(numel(Time),3);

for iTime=1:1:numel(Time);
  
  %compute the properties
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [MoonAz,MoonEl] = LunarAzEl(Time(iTime),Lat,Lon,Alt);
  Moon(iTime,[1,2]) = [MoonAz,MoonEl];
  clear MoonAz MoonEl
  
  [~,tidx] = min(abs(MoonPhase(1,:) - Time(iTime)));
  Moon(iTime,3) = MoonPhase(2,tidx);

end; clear iTime tidx Alt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iTime=1:1:numel(Time);
  
  %get star field
  [StarAz,StarEl,StarID] = find_star_azel(Lat,Lon,Time(iTime),MagCutOff);
  
  %got number of stars on first iteration - create array
  if iTime == 1; Stars = NaN(numel(Time),3,numel(StarAz)); end
  
  %fill array
  for iStar=1:1:numel(StarAz);
    Stars(iTime,:,iStar) = [StarID(iStar),StarAz(iStar),StarEl(iStar)];
  end
  
  clear StarAz StarEl StarID
  
end; clear iTime

%%done!
return
  