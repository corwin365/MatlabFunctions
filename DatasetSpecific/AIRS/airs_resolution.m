function [R,Z] = airs_resolution(Day,DoY,Lat,ZScale)

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
%
%out:
%  R - resolution at each height
%  Z - duplicate of requested Z scale
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%interpolate to chosen scale
R = interp1(Z,Profile,ZScale);
Z = ZScale;

%done!
return
end
