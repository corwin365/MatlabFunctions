function OUT = pw_hindley_v2(Data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reimplementation of Hindley23 PW filter as a function
%
%Corwin Wright, c.wright@bath.ac.uk, 2024/12/20
%note that this a direct reimplementation of Neil's older version, no logical changes have been made

%OUTPUTS:
%
%   OUT: struct containing the output fields, on a profile x height grid
%     FailReason - reason why no data present (0: data present; 1: MaxdX or Maxdt exceeded; 2: MinFracInProf not met ; 3: MindPhi exceeded; 4. input data was NaN; 5. no peak found )
%     Lat - latitude (deg)
%     Lon - longitude (deg)
%     Alt - altitude (km)
%     Time - time (Matlab units)
%     Trust - number of standard deviations from the edge of the data a point was fitted. Higher is better.
%     Temp_Residual - perturbaiton after PW removal (Kelvin)
%     Temp_PW - pws fitted, with an additional third dimension of PW mode (K)
%     Note - a written note.
%
%INPUTS:
%
%  required:
%    InstrumentData [struct] - data to use for the analysis.  This must contain at least the following:
%         Lat:          latitude of each point, in degrees
%         Lon:          longitude of each point, in degrees
%         Alt:          altitude of each point, units arbitrary but must be internally consistent
%         Time:         time, in Matlab datenum() units
%   The function get_limbsounders() will return data in this format, and is probably the easiest way to get it.
%
%  optional:
%
%    VarName         (type,                  default)  description
%    -----------------------------------------------------------------------------
%    TimeScale       (real,  all unique dates in data)  time scale for PW fit, i.e. unique periods to compute over
%    ZScale          (real, mean heightscale of input)  height scale for PW fit
%    LatScale        (real,                  -90:2:90)  lat scale for PW fit
%    FWHM            (real,                   [4,2,1])  FWHM of fit, in [deg lat, km alt, days time]
%    MinPoints       (integer,                     20)  Minimum points needed to fit a box
%    PWModes         (integer,                  1:1:9)  Specific PW modes to fit. Last will also contain residual (zonal mean)
%    Verbose         (logical,                   true)  Print progress to screen
%    NFWHMs          (integer,                      2)  Number of FWHMs to use data within in time, to control runtime
%    MinContrib      (real,                      0.05)  Minimum relative point weight to consider, to control runtime
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create input parser and testing functions
p = inputParser;
IsPositive        = @(x) validateattributes(x,{'numeric'},{'positive'});
IsPositiveInteger = @(x) validateattributes(x,{'numeric'},{'positive','integer'});

%inputs - required
addRequired(p,'Data',  @isstruct);  %just a type check

%inputs - optional but meaningful
addParameter(p,  'TimeScale',floor(min(Data.Time(:),[],'omitnan')):1:ceil(max(Data.Time(:),[],'omitnan')), @isnumeric);  %time scale for PW fit
addParameter(p,     'ZScale',mean(Data.Alt,1,'omitnan'), IsPositive);  %height scale for PW fit
addParameter(p,   'LatScale', -90:2:90, @isnumeric       );  %latitude scale for PW fit
addParameter(p,       'FWHM',  [4,2,1], IsPositive       );  %FWHM in lat [deg], alt [km],time [days]
addParameter(p,  'MinPoints',       20, IsPositiveInteger);  %number of points required for a sine fit
addParameter(p,    'PWModes',    1:1:9, IsPositiveInteger);  %PW modes to fit
addParameter(p,    'Verbose',     true,        @islogical);  %print out status updates

%inputs - optional timesavers, use at a cost in accuracy
addParameter(p,    'NFWHMs',    2, IsPositive       );  %number of FWHMs in time to include in fit data. Infinite is perfect, but is expensive.
addParameter(p,'MinContrib', 0.05, IsPositive       );  %minimum relative weight to contribute to a fit. Again, zero is best but expensive.

%parse inputs and tidy up
parse(p,Data,varargin{:});
Settings = p.Results; Settings = rmfield(Settings,'Data');
clear IsPositive p varargin IsPositiveInteger

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% other setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scale FWHMs to stdevs
Settings.Std = Settings.FWHM./2.355;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up fitting parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create storage arrays
FitData = struct();
Fields  = {'walt','wlat','wtim','walt_diff','wlat_diff','wlat_tim'};
for iF=1:1:numel(Fields); FitData.(Fields{iF}) = NaN(numel(Settings.ZScale), numel(Settings.LatScale),numel(Settings.TimeScale)); end
Fields  = {'A','B','C'};
for iF=1:1:numel(Fields); FitData.(Fields{iF}) = NaN(numel(Settings.ZScale), numel(Settings.LatScale),numel(Settings.TimeScale),numel(Settings.PWModes)); end
clear Fields iF

%check if we have a uniform z scale, and if so precompute weightings
dz = unique(diff(Data.Alt,1,2)); dz = dz(~isnan(dz));
if range(dz) < 1e-10; %i.e. the same to within numerical precision, given we're working with geophysical data
  ZWeights = cell(numel(Settings.ZScale),1);
  for iZ=1:1:numel(Settings.ZScale); ZWeights{iZ} = exp(-((Settings.ZScale - Settings.ZScale(iZ)).^2) ./ (2 * Settings.Std(2).^2)); end;
end
clear dz iZ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loop over time
%%%%%%%%%%%%%%%%
for iTime=1:1:numel(Settings.TimeScale)

  if Settings.Verbose == true
    disp(['----> Fitting PWs for ',datestr(Settings.TimeScale(iTime))]);
  end

  %trim data to within a narrow range of this particular date
  InTimeRange = inrange(mean(Data.Time,2,'omitnan'),Settings.TimeScale(iTime)+[-1,1].*Settings.FWHM(3).*Settings.NFWHMs);
  Working = reduce_struct(Data,InTimeRange,{},1);

  %determine relative weights of each time point
  w_tim = exp(- ((Working.Time - Settings.TimeScale(iTime)).^2) ./ (2 * Settings.Std(3)^2));

  %loop over latitude
  %%%%%%%%%%%%%%%%%%%%
  for iLat=1:1:numel(Settings.LatScale)

    %determine relative weights of each lat point
    w_lat = exp(- ((Working.Lat - Settings.LatScale(iLat)).^2) ./ (2 * Settings.Std(1)^2));

    %loop over altitude
    %%%%%%%%%%%%%%%%%%%%
    for iZ=1:1:numel(Settings.ZScale)

      %check if we have predetermined z weightings, and use these if they exist
      %otherwise, compute new weights
      if exist('ZWeights','var'); w_alt = repmat(ZWeights{iZ},size(Working.Alt,1),1);
      else;                       w_alt = exp(- ((Working.Alt - Settings.ZScale(iZ)).^2) ./ (2 * Settings.Std(2)^2));
      end

      %combine the weightings
      w  = w_tim .* w_lat .*w_alt;
      idx = find(w > Settings.MinContrib);

      %time to do the weighted fit
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %do we have enough points to do the fit?
      T = Working.Temp(idx); T = T(~isnan(T));
      if numel(T) < Settings.MinPoints; clear w idx T w_alt; continue; end

      %we do! Do it. Turn warnigns off, they're going to happen often but we're going to handle them here.
      warning off
      [~,F] = nph_sinefit(Working.Lon(idx),Working.Temp(idx),(360./Settings.PWModes),'weights',w(idx));
      %if the fit had too few points, then throw it away
      if any(regexp(lastwarn,'Matrix is close to singular or badly scaled.')); F = NaN(length(Settings.PWModes),3); end
      lastwarn('warning_reset')
      warning on

      %store the results
      FitData.A(iZ,iLat,iTime,:) = F(:,1); 
      FitData.B(iZ,iLat,iTime,:) = F(:,2); 
      FitData.C(iZ,iLat,iTime,:) = F(:,3); 

      %store the weighted mean locations...
      FitData.wlat(iZ,iLat,iTime) = wmean(Working.Lat( idx),w(idx));
      FitData.walt(iZ,iLat,iTime) = wmean(Working.Alt( idx),w(idx));
      FitData.wtim(iZ,iLat,iTime) = wmean(Working.Time(idx),w(idx));

      %... and their deviations from bin centre
      FitData.wlat_diff(iZ,iLat,iTime) = Settings.LatScale(  iLat) - FitData.wlat(iZ,iLat,iTime);
      FitData.walt_diff(iZ,iLat,iTime) = Settings.ZScale(      iZ) - FitData.walt(iZ,iLat,iTime);
      FitData.wtim_diff(iZ,iLat,iTime) = Settings.TimeScale(iTime) - FitData.wtim(iZ,iLat,iTime);    


      clear w idx T F w_alt
    end; clear iZ
  end; clear iLat w_lat   
end; clear iTime InTimeRange Working w_tim ZWeights



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now we need to interpolate back to the original grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%there are lots of caveats on the logic here, which I reproduce from Neil's original
%code at the very bottom of this file, below the functions.
%you should probably read them.


%find the average of the true weighted centres along each scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avg_walt = squeeze(mean(FitData.walt,[2,3],'omitnan'));
avg_wlat = squeeze(mean(FitData.wlat,[1,3],'omitnan'));
avg_wtim = squeeze(mean(FitData.wtim,[1,2],'omitnan'));

%replace any NaNs with their scale values (these have no data anyway, this is just to stabilise the numerics)
avg_walt(isnan(avg_walt)) = Settings.ZScale(   isnan(avg_walt));
avg_wlat(isnan(avg_wlat)) = Settings.LatScale( isnan(avg_wlat));
avg_wtim(isnan(avg_wtim)) = Settings.TimeScale(isnan(avg_wtim));

%replace extrema with the true data limits, to remove any risk of extrapolation
avg_walt([1,end]) = Settings.ZScale(   [1,end]);
avg_wlat([1,end]) = Settings.LatScale( [1,end]);
avg_wtim([1,end]) = Settings.TimeScale([1,end]);

%create interpolants using dummy vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create fitting interpolants using dummy vectors
sz = size(FitData.A);
F = struct();
for iPW=1:1:numel(Settings.PWModes)
  F(iPW).A = griddedInterpolant({1:sz(1),1:sz(2),1:sz(3)},FitData.A(:,:,:,iPW),'linear','none');
  F(iPW).B = griddedInterpolant({1:sz(1),1:sz(2),1:sz(3)},FitData.B(:,:,:,iPW),'linear','none');
  F(iPW).C = griddedInterpolant({1:sz(1),1:sz(2),1:sz(3)},FitData.C(:,:,:,iPW),'linear','none');
end; clear iPW sz

%now estimate PWs, treating each unit of time (i.e. each unit of diff(Settings.TimeScale)) individually 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create output struct, duplicating the format Neil used previously
OUT = spawn_uniform_struct({'Lat','Lon','Alt','Time','Temp_Residual'},size(Data.Lat));
OUT.Temp_PW = NaN([size(Data.Lat),numel(Settings.PWModes)]);
OUT.Note = 'Note that the final PW mode (:,:,end) in Temp_PW also contains the zonal mean. Generated using 2024/12/20 code version.';

%process
for iTime=1:1:numel(Settings.TimeScale)-1

  if Settings.Verbose == true
    disp(['----> Evaluating PWs for ',datestr(Settings.TimeScale(iTime))]);
  end

  %get the data for this time unit
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  idx = inrange(mean(Data.Time,2,'omitnan'),Settings.TimeScale([iTime,iTime+1]));
  if numel(idx) == 0; continue; end

  %interpolate the sample locations in units of the dummy vectors
  %this is vile.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Samples.walt = interp1(linspace(avg_walt(1),avg_walt(end),numel(avg_walt)), ...
                                 1:1:numel(avg_walt), ...
                                 Data.Alt(idx,:));

  Samples.wlat = interp1(linspace(avg_wlat(1),avg_wlat(end),numel(avg_wlat)), ...
                                 1:1:numel(avg_wlat), ...
                                 Data.Lat(idx,:));

  Samples.wtim = interp1(linspace(avg_wtim(1),avg_wtim(end),numel(avg_wtim)), ...
                                 1:1:numel(avg_wtim), ...
                                 Data.Time(idx,:));  

  %and hence reconstruct the PWs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for iPW=1:1:numel(Settings.PWModes)
    OUT.Temp_PW(idx,:,iPW) = F(iPW).A(Samples.walt,Samples.wlat,Samples.wtim) .* cos(2.*pi.*(Data.Lon(idx,:))./(360./Settings.PWModes(iPW))) + ...
                             F(iPW).B(Samples.walt,Samples.wlat,Samples.wtim) .* sin(2.*pi.*(Data.Lon(idx,:))./(360./Settings.PWModes(iPW))) + ...
                             F(iPW).C(Samples.walt,Samples.wlat,Samples.wtim);
  end; clear iPW

  %write to output arrays
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  OUT.Lat( idx,:) = Data.Lat(idx,:);
  OUT.Lon( idx,:) = Data.Lon(idx,:);
  OUT.Alt( idx,:) = Data.Alt(idx,:);
  OUT.Time(idx,:) = Data.Time(idx,:);

  OUT.Temp_Residual(idx,:) = Data.Temp(idx,:) - sum(OUT.Temp_PW(idx,:,:),3);

  clear idx

end; clear iTime Samples avg_walt avg_wlat avg_wtim


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finally: how much do we trust the results?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create an array telling us how much to trust the results
%this is just how many FWHMs the point is from the edge in its WORST dimension
%higher values are better, since they are further inside the data field

TimeMax =  (max(Settings.TimeScale) - Data.Time) ./ Settings.FWHM(3);
TimeMin = -(min(Settings.TimeScale) - Data.Time) ./ Settings.FWHM(3);
Trust = cat(3,TimeMax,TimeMin);
AltMax  =  (max(Settings.ZScale) - Data.Alt) ./ Settings.FWHM(2);
Trust = cat(3,Trust,AltMax);
AltMin  = -(min(Settings.ZScale) - Data.Alt) ./ Settings.FWHM(2);
Trust = cat(3,Trust,AltMin);
LatMax  =  (max(Settings.LatScale) - Data.Lat) ./ Settings.FWHM(1);
Trust = cat(3,Trust,LatMax);
LatMin  = -(min(Settings.LatScale) - Data.Lat) ./ Settings.FWHM(1);
Trust = cat(3,Trust,LatMin);
Trust = min(Trust,[],3);
Trust(isnan(Trust)) = -999; %really bad trust if it's a NaN!
clear TimeMax TimeMin AltMax AltMin LatMax LatMin


OUT.Trust = Trust; clear Trust


return
end













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BELOW BE FUNCTIONS. I WOULDN'T VENTURE DOWN THERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










function y = wmean(x,w,dim)
%WMEAN   Weighted Average or mean value.
%   For vectors, WMEAN(X,W) is the weighted mean value of the elements in X
%   using non-negative weights W. For matrices, WMEAN(X,W) is a row vector 
%   containing the weighted mean value of each column.  For N-D arrays, 
%   WMEAN(X,W) is the weighted mean value of the elements along the first 
%   non-singleton dimension of X.
%
%   Each element of X requires a corresponding weight, and hence the size 
%   of W must match that of X.
%
%   WMEAN(X,W,DIM) takes the weighted mean along the dimension DIM of X. 
%
%   Class support for inputs X and W:
%      float: double, single
%
%   Example:
%       x = rand(5,2);
%       w = rand(5,2);
%       wmean(x,w)

% EDIT: NPH 2017: Added nansum() rather than sum() to make it nan safe, also
% removed that requirement that all weights must be non-zero, it just spits
% out a NaN now if that's the case rather than an error.

if nargin<2
    error('Not enough input arguments.');
end

% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    error('Inputs x and w must be the same size.');
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
% if (all(w(:)==0))
%     error('At least one weight must be non-zero.');
% end

if nargin==2
  % Determine which dimension SUM will use
  dim = find(size(x)~=1, 1 );
  if isempty(dim), dim = 1; end
end

y = sum(w.*x,dim,'omitnan')./sum(w,dim,'omitnan');

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%caveats referred to above relating to interpolation back to the obs locations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to remove the PW perturbations from our original profiles, we need to
% interpolate the fitted PWs back on to the lat/lon/time/height of the
% original data.

% to do this, we create interpolants of our fit coefficients A, B and C and
% then evaluate these at every lat/time/height. Then, using these values
% and longitude, we evalute the temperature perturbation at every location.

% why am I doing it this way? because you avoid having to build and evaluate
% a 4-D interpolant, which could get huuuge.

% CAVEAT: note that the fits of PWs weren't always exactly fitted where you
% centred the Gaussian, i.e. the weighted average of all the times of the
% data might not actually be the central time of the weighting. This means
% that your fit was for a slightly different time than you thought. I
% solved this fully for the meteor radar work, but I'm relucatant here
% because it will mean a very large scatteredInterpolant object rather than
% a quicker griddedInterpolant. We'll see how we get on.

% EDIT: A quick spot check with GNSS-RO suggests:
% - For time, sometimes we can be up to an hour off. This could get
% better/worse for some sunsynchronous satellites, but is worth considering.
% - For latitude bins, looks like we're typically less than half a degree out,
% but there's some interesting zonal patterns in this with the GNSS-RO data.
% - For time (revised), looks like on some days this can be as much as 2-3
% hours out, presumably depending on the distribution of occultations that
% day versus the prev/next.

% % % % % Let's try the scatteredInterpolant with the true weighted locations and
% % % % % see how we get on. - EDIT: This was unbearably slow!!!!
% % % % %
% % % % % F = scatteredInterpolant(FitData.walt(:),FitData.wlat(:),FitData.wtim(:),...
% % % % %     ones(size(linearise(sq(FitData.A(:,:,:,1))))),'linear');
% % % % % % F.Values = linearise(sq(FitData.A(:,:,:,1)));
% % % % % Ai = F(altrep,latrep,timrep);

% New plan: let's use the AVERAGE weighted fit locations and a gridded
% interpolant. The latitude offsets seems to be consistent each day, and
% the time offsets seem to be consistent each latitude. Roughly.
% Note that thereis a trick here - we're going to make matlab think it's
% got a regularly gridded interpolant but actually we're gonna then sample
% it irregularly, which is faster.