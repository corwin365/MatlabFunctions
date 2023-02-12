function [R,Z] = airs_resolution(Night,DoY,Lat,ZScale)

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
%  Night: 0 for daytime profile, 1 for nighttime
%  DoY: day of year. Can be fractional
%  Lat: latitude
%  ZScale: height scale to interpolate results onto, in km
%
%out:
%  R - resolution at each height
%  Z - duplicate of requested Z scale
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load([LocalDataDir,'/AIRS/airs3d_resolution.mat'])

%get resolution profile
[~,idx1] = min(abs(DoY-DayOfYear));
[~,idx2] = min(abs(Lat-Latitude));
idx3  = Night+1; %first element is day, second is night

Contribs = squeeze(RetrievalMix(idx1,idx2,idx3,:));

Profile = Profiles(:,1)*Contribs(1);
for iProf=2:1:8; Profile = Profile +  Profiles(:,iProf)*Contribs(iProf); end
Profile = Profile./sum(Contribs);

%interpolate to chosen scale
R = interp1(Z,Profile,ZScale);
Z = ZScale;

%done!
return
end
