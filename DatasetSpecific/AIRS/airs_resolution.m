function [R,Z,N,O] = airs_resolution(Day,DoY,Lat,ZScale,v1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to compute a weighted resolution profile
%for the Hoffmann and Alexander (2009), based on 
%linear dependence on latitude, day-of-year, and
%whether it's day or night
%
%requires: 
%     'airs3d_resolution.mat' - file containing the coefficients and baseline profiles
%
%in:
%  Day: 1 for daytime profile, 0 for nighttime.
%  DoY: day of year. Can be fractional
%  Lat: latitude
%  ZScale: height scale to interpolate results onto, in km
%  v1: optional - set to 1 to use the old 8-profile version, otherwise uses new ERA5-based version
%
%out:
%  R - resolution at each height
%  Z - duplicate of requested Z scale
%  N - noise estimate at each point
%  O - estimated signal fraction at each point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('v1'); v1 = 0; elseif v1 ~= 1; v1 = 0; end


if v1 == 1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% old version
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %load data
  load([LocalDataDir,'/AIRS/airs3d_resolution.mat'])

  %stick very large values where the resoution is poorly-defined
  Profiles(Profiles == 0) = 30;

  %get resolution profile
  [~,idx1] = min(abs(DoY-DayOfYear));
  [~,idx2] = min(abs(Lat-Latitude));
  if Day == 1; idx3 = 2; else idx3 = 1; end %1 is nighttime, 2 is daytime

  Contribs = squeeze(RetrievalMix(idx1,idx2,idx3,:));

  Profile = Profiles(:,1)*Contribs(1);
  for iProf=2:1:8; Profile = sum([Profile,Profiles(:,iProf)*Contribs(iProf)],2,'omitnan'); end
  Profile = Profile./sum(Contribs);

  %for forwards compatability, add outputs for noise and sig frac set to NaN
  Noise   = NaN(size(Profile));
  ObsFrac = NaN(size(Profile));


else

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% new version
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %load data
  Grid = load([LocalDataDir,'/AIRS/airs3d_resolution_e5.mat']);

  %find gridbox of interest
  [~,idx1] = min(abs(Grid.DoY-DoY));
  [~,idx2] = min(abs(Grid.Lat-Lat));
  [~,idx3] = min(abs(Grid.Day-Day));

  %get data in format for interpolation
  Z = Grid.Z;
  Profile = squeeze(Grid.Res(   idx1,idx2,idx3,:));
  Noise   = squeeze(Grid.Noise( idx1,idx2,idx3,:));
  ObsFrac = squeeze(Grid.Signal(idx1,idx2,idx3,:));
  

end


%interpolate to chosen scale
R = interp1(Z,Profile,ZScale);
N = interp1(Z,Noise,ZScale);
O = interp1(Z,ObsFrac,ZScale);
Z = ZScale;

%done!
return
end
