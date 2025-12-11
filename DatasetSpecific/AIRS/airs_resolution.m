function [R,Z,N,O] = airs_resolution(Day,DoY,Lat,ZScale,Format,Channel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to compute a weighted resolution profile
%for the Hoffmann and Alexander (2009), based on 
%linear dependence on latitude, day-of-year, and
%whether it's day or night
%
%updated 2025/12/02 to:
%  (a) provide 2D AIRS channels or 3D retrieved as options
%  (b) remove legacy flag for original version
%
%requires file containing the coefficients and baseline profiles:
%     'airs3d_resolution_e5.mat' and 'airs2d_resolution.mat'
%
%in:
%  Day: 1 for daytime profile, 0 for nighttime. Ignored for 2D data.
%  DoY: day of year. Can be fractional. Ignored for 2D case
%  Lat: latitude
%  ZScale: height scale to interpolate results onto, in km
%  Format: optional - set to '2' to use 2D brightness temperature kernels at fixed locations, 4 for 2D kernels interpolated in space and time, otherwise 3D used
%  Channel: '4mu'	'15mu_high' or	'15mu_low'. Only used in 2D case
%
%out:
%  R - resolution at each height
%  Z - duplicate of requested Z scale
%  N - noise estimate at each point (only for 3D, NaN in 2D case)
%  O - estimated signal fraction at each point (only for 3D, NaN in 2D case)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%handle inputs
if ~exist('Format','var');
  Format = 3; %i.e. 3D
elseif Format == 1;
  error('Error: this is no longer a valid option for airs_resolution()') %just in case
elseif Format == 2 | Format == 4;
  %we're in the 2D regime. What channel?
  if nargin < 6
    error('Error: 2D weighting functions requested but no channel requested: valid options are ''4mu'', ''15mu_high'' or	''15mu_low'' ')
  elseif  strcmpi(Channel,      '4mu'); idx2 = 1;
  elseif  strcmpi(Channel,'15mu_high'); idx2 = 2;
  elseif  strcmpi(Channel, '15mu_low'); idx2 = 3;
  else
    error('Error: 2D weighting functions requested but invalid channel requested: valid options are ''4mu'', ''15mu_high'' or	''15mu_low'' ')
  end
else
  Format = 3; %i.e. 3D
end




if Format == 3;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% 3D retrieval kernels
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
  

elseif Format == 2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% 2D brightness temperature kernels
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %load data
  Grid = load([LocalDataDir,'/AIRS/airs2d_resolution.mat']);

  %work out which geographic regime we're in
  if     abs(Lat) < 25; idx1 = find(strcmp(Grid.Region, 'Equatorial')); 
  elseif abs(Lat) < 60; idx1 = find(strcmp(Grid.Region, 'Midlatitudes'));
  else 
    %polar regime. This is more complex as it depends on the time of year
    if DoY > 80 && DoY < 265
      if   Lat < 0; idx1 = find(strcmp(Grid.Region, 'PolarWinter'));
      else          idx1 = find(strcmp(Grid.Region, 'PolarSummer'));  
      end
    else
      if   Lat > 0; idx1 = find(strcmp(Grid.Region, 'PolarWinter'));
      else          idx1 = find(strcmp(Grid.Region, 'PolarSummer'));  
      end 
    end
  end

  %extract the data
  Z       = Grid.Z;  
  Profile = Grid.Resolution(:,idx1,idx2);
  Noise   = NaN(size(Z));
  ObsFrac = NaN(size(Z));

  %normalise  
  Profile = Profile./nansum(Profile);

elseif Format == 4

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% 2D brightness temperature kernels, interpolated
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %load data
  Grid = load([LocalDataDir,'/AIRS/airs2d_resolution.mat']);

  %first, interpolate in time in the polar regions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %produce a time weighting taper that is structured as:
    % DJF - 0
    % MAM - linear from 0 to 1
    % JJA - 1
    % SON - linear from 1 to 0
  Taper = NaN(365,1);
  Taper(1:date2doy(datenum(1,2,28))) = 0; % JF
  Taper(date2doy(datenum(1,12,1)):end) = 0; %D
  Taper(date2doy(datenum(1,6,1)):date2doy(datenum(1,8,31))) = 1; %JJA
  %MAM:
  mam = date2doy(datenum(1,3,1)):date2doy(datenum(1,5,31));
  Taper(mam) = linspace(0,1,numel(mam));
  %SON:
  son = date2doy(datenum(1,9,1)):date2doy(datenum(1,11,30));
  Taper(son) = linspace(1,0,numel(son));
  clear mam son

  %work out which season is summer and which is winter
  if Lat > 0; 
    jja_idx = find(strcmp(Grid.Region, 'PolarSummer'));
    djf_idx = find(strcmp(Grid.Region, 'PolarWinter'));
  else
    djf_idx = find(strcmp(Grid.Region, 'PolarSummer'));
    jja_idx = find(strcmp(Grid.Region, 'PolarWinter'));   
  end

  %hence, produce a time-weighted profile and renormalise
  Polar = ((1-Taper(closest(DoY,1:1:365))).*Grid.Resolution(:,djf_idx,idx2)) ...
        + ((  Taper(closest(DoY,1:1:365))).*Grid.Resolution(:,jja_idx,idx2));
  Polar = Polar./sum(Polar);

  %now, weight the regions in latitude
  %remember, we care about RELATIVE weight, not absolute, at each latitude
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  % % % %define a taper in a similar way that gives absolute and mixed regions
  lats = -90:0.1:90;
  tropics = cosd(2.*lats);          tropics(tropics < 0) = 0;
  midlats = sind(abs(2.*lats));     midlats(midlats < 0) = 0;
  polar = sind(abs(2.*lats) - 90);  polar(    polar < 0) = 0;
  weights(:,1) = polar;
  weights(:,2) = midlats;
  weights(:,3) = tropics;



  %hence, find the relative weights we want
  idx = closest(Lat,lats); weights = weights(idx,:); weights = weights./sum(weights);

  %and finaly produce a weighted resolution profile, including renormalisation
  Profile = weights(1)*Polar                                                            + ...
            weights(2)*Grid.Resolution(:,find(strcmp(Grid.Region,'Midlatitudes')),idx2) + ...
            weights(3)*Grid.Resolution(:,find(strcmp(Grid.Region,  'Equatorial')),idx2);

  %done!
  Z       = Grid.Z;  
  Noise   = NaN(size(Z));
  ObsFrac = NaN(size(Z));
 


end

%interpolate to chosen scale
R = interp1(Z,Profile,ZScale);
N = interp1(Z,Noise,ZScale);
O = interp1(Z,ObsFrac,ZScale);
Z = ZScale;

%done!
return
end
