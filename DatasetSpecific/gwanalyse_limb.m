function [OutData,PW] = gwanalyse_limb(Data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%general-purpose function to measure GWs using vertical profile data
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/10/29
%
%
%===================================
%OUTPUTS:
%
%   OutData: struct containing the output fields, on a profile x height grid
%     FailReason - reason why no data present (0: data present; 1: MaxdX or Maxdt exceeded; 2: MinFracInProf not met ; 3: MindPhi exceeded; 4. input data was NaN )
%     A - amplitude (input units)
%     Lz - vertical wavelength (km)
%     Lh - horizontal wavelength (km)
%     Lat - latitude (deg)
%     Lon - longitude (deg)
%     Alt - altitude (input units)
%     Tp - temperature perturbation (input units)
%     MF  - momentum flux (mPa)
%     MF and Lh will be present but all NaN if 'filter' is set to anything other than 2
%
%   PWs:     array containing PWs. Format varies depending on filter type used.
%
%===================================
%
%
%INPUTS:
%
%  required:
%    InstrumentData [struct] - data to use for the analysis.  This must contain at least the following:
%         Lat:          latitude of each point, in degrees
%         Lon:          longitude of each point, in degrees
%         Alt:          altitude of each point, units arbitrary but must be internally consistent
%         Time:         time, in Matlab datenum() units
%
%   For all filters except 'Hindley23', we also need:
%         Temp:         temperature, units arbitrary but must be internally consistent
%
%   If using the 'Hindley23' planetary wave filter, which is computed externally, then we instead need:
%         Temp_Residual: residual temperatures after the PW filtering, also arbitrary-but-consistent
%         Temp_PW: planetary wave temperatures for each mode at each location, same units
%      You probably don't want to do this, but we can also use these inputs for the other filters instead of 'Temp'- the
%      routine will merge them to produce 'Temp' if so and then proceed normally.
%
%
%  optional:
%
%    VarName         (type,           default)  description
%    -----------------------------------------------------------------------------
%    Analysis        (integer,              2)  type of analysis to use, from (1) 1DST, (2) Alexander et al 2008 
%    Filter          (char,          'SGolay')  type of detrending filter to use (see below)
%    STScales        (vector,   1:1:NLevels/2)  number of scales to use in 1D ST
%    STc             (positive real,        1)  value of 'c' to use in ST
%    STPadSize       (positive,            20)  levels of zero-padding to put at each end of the data before S-Transforming
%    RegulariseZ     (logical            true)  interpolate the data to a regular height grid
%    Verbose         (logical,           true)  report to the user what's happening
%
%-----------------------------------
%if 'Analysis' is set to 2, then the following options can be used:
%    MindPhi         (positive real,    0.025)  minimum fractional phase change permitted for Alexander et al 2008 analysis 
%    MaxdX           (positive real,      400)  maximum interprofile distance [km] permitted for Alexander et al 2008 analysis 
%    Maxdt           (positive real,      900)  maximum interprofile time [s] permitted for Alexander et al 2008 analysis 
%    MinFracInProf   (positive real,      0.5)  minimum fraction useful levels remaining in a profile after above two filters to proceed
%-----------------------------------
%if 'Filter' is set to 'Hindley23', then the following options can be used:
%    H23_OutRem (logical,                true)   remove outliers from H23 fitting, as in get_limbsounders()
%-----------------------------------
%if 'Filter' is set to 'SGolay', then the following options can be used:
%    SGOrder         (integer,              2)  order of Savitzky-Golay filter to use
%    SGLength        (real,                18)  frame length of SG filter, in km
%When using this filter singularities like the tropopause must be removed in advance
%-----------------------------------
%if 'Filter' is set to 'PWgrid', then the following options can be used:
%    NPWs            (integer,              3)  number of PWs to fit (plus zonal mean)
%    PWWindow        (positive real,        1)  number of days to use in each PW fit
%    PWTimeRes       (positive real,        1)  time resolution to export PW fits on
%    PWLonGrid       (vector      -180:20:180)  longitude grid to compute PWs on
%    PWMinPC         (positive real,       66)  fraction of longitude bins which must be filled for each lat band
%    PWLatGrid       (vector         -90:5:90)  longitude grid to compute PWs on
%    PWAltGrid       (vector, input grid mean)  altitude grid to compute PWs on
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create parser object
%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

%validation of inputs
%%%%%%%%%%%%%%%%%%%%%%

%data structure
addRequired(p,'Data',@isstruct); %input data must be a struct. Will be hand-parsed later.

%analysis to use
addParameter(p,'Analysis',2,@(x) validateattributes(x,{'numeric'},{'integer','<=',2})) %type of analysis to use (see above)

%filter to use
addParameter(p,'Filter','SGolay',@ischar) %type of filter to use

%ST properties
addParameter(p,'STScales',1:1:size(Data.Alt,2)/2,@isvector  ) %scales to compute on ST
addParameter(p,'STc',                          1,@ispositive) %'c' parameter for ST
addParameter(p,'STPadSize',                   20,@ispositive) %levels of zero-padding to put at each end of the data before S-Transforming

%Alex08 horizontal wavelength properties
addParameter(p,'MaxdX',         300,@ispositive) %maximum distance between profiles
addParameter(p,'Maxdt',         900,@ispositive) %maximum time between profiles
addParameter(p,'MinFracInProf', 0.5,@ispositive) %maximum fraction of profile remaining after above filters
addParameter(p,'MindPhi',       0.025,@ispositive) %minimum fractional phase change to compute Lh
 
%verbosity
addParameter(p,'Verbose',  true,@islogical) %print progress updates?


%interpolate to a regular height grid? (you probably want to keep this set to true - it checks if it's already true,
%and only does the interpolation if needed)
addParameter(p,'RegulariseZ',true,@islogical)

%additional variables - planetary wave grid fit filter ('PWGrid')
addParameter(p,      'NPWs',                  3,@isreal) %number of planetary waves to fit
addParameter(p,   'PWMinPC',                 66,@(x) validateattributes(x,{'numeric'},{'>',0})) %fraction of longitude bins that must be filled
addParameter(p,  'PWWindow',                  1,@(x) validateattributes(x,{'numeric'},{'>',0})) %window width of period used, in days
addParameter(p, 'PWTimeRes',                  1,@(x) validateattributes(x,{'numeric'},{'>',0})) %time resolution
addParameter(p, 'PWLonGrid',        -180:20:180,@(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',180,'vector'})) %longitude bins
addParameter(p, 'PWLatGrid',           -90:5:90,@(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90,'vector'})) %latitude bins
addParameter(p, 'PWAltGrid',nanmean(Data.Alt,1),@(x) validateattributes(x,{'numeric'},{'>=',   0,         'vector'})) %altitude bins

%additional parameters - SGolay filter ('SGolay')
addParameter(p,'SGOrder',  2,@(x) validateattributes(x,{'numeric'},{'>',0})) %SGolay filter order
addParameter(p,'SGLength',18,@(x) validateattributes(x,{'numeric'},{'>',0})) %SGolay filter length (km)

%additional parameters - Hindley 2023 filter ('Hindley23')
addParameter(p,'H23_OutRem',  true,@islogical) %remove outliers?

%parse inputs and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,Data,varargin{:})

%pull out the remaining arguments into struct "Settings", used throughout rest of routine
Settings = p.Results;
Data = Settings.Data; Settings = rmfield(Settings,'Data');
Settings.TimeRange = [floor(nanmin(Data.Time,[],'all')),ceil(nanmax(Data.Time,[],'all'))];
clearvars -except InstInfo Settings Data

%check contents of input data struct:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%if we have loaded one of the Hindley23 data files but are not using the Hindley 23 filter (not sure why you'd do this 
% normally, but might be useful for cross-testing) then merge the temperature and PW fields to produce an integrated product
if ~strcmpi(Settings.Filter,'Hindley23') && ~isfield(Data,'Temp') && isfield(Data,'Temp_PW') && isfield(Data,'Temp_Residual');
  Data.Temp = sum(Data.Temp_PW,3) + Data.Temp_Residual;
  Data = rmfield(Data,{'Temp_PW','Temp_Residual','Note'});
end

%%do we have at least Lat, Lon, Alt? These are used by all filters
if ~isfield(Data,'Lat') ||  ~isfield(Data,'Lon') || ~isfield(Data,'Alt')
  error('Missing field - struct must contain Lat, Lon, and Alt fields.')
  return
end

%for all except Hindley23, we also need a 'Temp' field containing temperature
if ~strcmpi(Settings.Filter,'Hindley23') &  ~isfield(Data,'Temp');
  error('Missing field - struct must contain Temp field.')
else Data.Temp = NaN(size(Data.Lon)); end %this will be ignored anyway in the analysis

%%are these fields all the same size?
if ~isequal(size(Data.Lat),size(Data.Lon)) ||  ~isequal(size(Data.Lat),size(Data.Alt)) || ~isequal(size(Data.Lat),size(Data.Temp));
  error('Lat, Lon, Alt and Temp fields must all be the same size.')
  return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ensure the data is on a regular height grid
% for Hindley23 removal, this must be done AFTER detrending, as new data are loaded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.RegulariseZ == true && ~strcmpi(Settings.Filter,'Hindley23'); Data = regularise_data_z(Data,Settings); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% detrend the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tg

%start by duplicating temperatures. All filter functions will act on this,
%in order to let us combine filters if we want
Data.Tp = Data.Temp;

%now, apply the filter. Currently only one at a time, but simple to rewrite in the future to apply multiple filters if needed

if     strcmpi(Settings.Filter,   'PWgrid'); [Data,PW] = filter_pwgrid(   Data,Settings);
elseif strcmpi(Settings.Filter,   'SGolay'); Data      = filter_sgolay(   Data,Settings);
elseif strcmpi(Settings.Filter,'Hindley23'); [Data,PW] = filter_hindley23(Data,Settings); 
else
  error('Filter type not included in programme, stopping')
  return
end


%create empty PW output if we didn't generate one. Note that this output will be 
%fairly meaningless anyway if the PW filter wasn't the first-applied
if ~exist('PW','var'); PW.Comment = 'Planetary wave filter not used, no PW data computed'; end


%regularisation if using Hindley23 detrending
if Settings.RegulariseZ == true && strcmpi(Settings.Filter,'Hindley23'); Data = regularise_data_z(Data,Settings); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S-Transform profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%produce storage arrays
NProfiles = size(Data.Tp,1);
NLevs     = size(Data.Tp,2);
OutData   = spawn_uniform_struct({'A','Lz','Lh','Lat','Lon','Alt','Tp','MF','Time','FailReason'},[NProfiles,NLevs]);
Mask      = ones([NProfiles,NLevs]); %this is used to mask out bad data later

%some approaches require two adjacent profiles to be computed. To avoid duplicate computation in this case,
%it is marginally more efficient if we work backwards and store the 'previous' (i.e. next) profile to permit this.

if Settings.Verbose == 1; textprogressbar('--> Computing gravity waves '); end
for iProf=NProfiles:-1:1

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %S-Transform the profile
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %fill any NaNs with a zero for purposes of the ST
  %we'll replace with NaNs at the end of the programme
  NoData = find(isnan(Data.Tp(iProf,:)));
  Tp = Data.Tp(iProf,:); Tp(NoData) = 0;
  Mask(iProf,NoData) = 0;

  %zero-pad the data to prevent FFT wraparound
  Tp = [zeros(1,Settings.STPadSize),Tp,zeros(1,Settings.STPadSize)];

  %compute ST
  ThisST = nph_ndst(Tp,                                ...
                    Settings.STScales,                 ...
                    nanmean(diff(Data.Alt(iProf,:))),  ...
                    Settings.STc);
  clear Tp 

  %remove the zero-padding region from all output variables
  Fields = {'IN','F1','A','R','HA','HR','allgws','BoostFactor','C'};
  for iF=1:1:numel(Fields); 
    F = ThisST.(Fields{iF});
    F = F(Settings.STPadSize+1:end-Settings.STPadSize);
    ThisST.(Fields{iF}) = F;
  end; clear Fields iF F
  ThisST.ST = ThisST.ST(:,Settings.STPadSize+1:end-Settings.STPadSize);


  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %compute and store results depending on analysis
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if Settings.Analysis == 1

    %simple 1DST approach - just store
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OutData.A(   iProf,:) = ThisST.A;
    OutData.Lz(  iProf,:) = 1./ThisST.F1;
    OutData.Lat( iProf,:) = Data.Lat( iProf,:); OutData.Lon(iProf,:) = Data.Lon(iProf,:);
    OutData.Alt( iProf,:) = Data.Alt( iProf,:); OutData.Tp( iProf,:) = Data.Tp( iProf,:);
    OutData.Time(iProf,:) = Data.Time(iProf,:); OutData.FailReason(iProf,:) = 0;
    

  elseif Settings.Analysis == 2

    %Alexander et al (JGR, 2008) approach
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %if this is the first profile to be processed, retain the data then move to the next
    if iProf == NProfiles; NextST = ThisST; continue; end
    
    %assume success until proven otherwise
    OutData.FailReason(iProf,:) = 0;

    %compute cospectrum
    CoSpectrum = ThisST.ST .* conj(NextST.ST);

    %find time and space separation of the two profiles, 
    dx = nph_haversine([Data.Lat(iProf,  :);Data.Lon(iProf,  :)]', ...
                       [Data.Lat(iProf+1,:);Data.Lon(iProf+1,:)]')';
    dt = abs(Data.Time(iProf,:) - Data.Time(iProf+1,:)).*60.*60.*24;
    Bad = find(dt > Settings.Maxdt | dx > Settings.MaxdX);

    %discard the profile completely if we fall below the minimum acceptable fraction of safe data
    if numel(Bad)/numel(dx) > 1-Settings.MinFracInProf; 
      OutData.FailReason(iProf,:) = 2; 
      clear Bad dt dx CoSpectrum; 
      continue; 
    end

    %if we pass the above, discard any heights where we failed either individually
    if numel(Bad) > 0; CoSpectrum(:,Bad) = NaN; OutData.FailReason(iProf,Bad) = 1; end

    %drop modal frequency
    CoSpectrum(1,:) = NaN;

    %locate maximum at each height
    [A,idx] = nanmax(CoSpectrum,[],1,'omitnan');
    A = sqrt(abs(A));
    
    %find vertical waveLENGTHS
    Lz = 1./ThisST.freqs(idx);

    %find phase change, and discard values where it's too small to be meaningful
    dx(dx > Settings.MaxdX) = NaN;
    dPhi = angle(CoSpectrum(idx))./(2*pi);

    OutData.FailReason(iProf,dPhi < Settings.MindPhi) = 3; 
    dPhi(dPhi < Settings.MindPhi) = NaN; 

    %compute Lh and MF
    Lh = abs(dx./dPhi);
    MF = cjw_airdensity(Data.Pres(iProf,:),Data.Temp(iProf,:))  ...
      .* (Lz ./ Lh)                                            ...
      .* (9.81./0.02).^2                                        ...
      .*A./(Data.Temp(iProf,:)-Data.Tp(iProf,:)).^2;


    %adjust lat/lon to the midpoint of the profile-pair
    [latmean,lonmean] = meanm(Data.Lat(iProf+[0,1],:), ...
                              Data.Lon(iProf+[0,1],:));

    %store results
    OutData.A(   iProf,:) = A;
    OutData.Lz(  iProf,:) = Lz;
    OutData.Lh(  iProf,:) = Lh;
    OutData.MF(  iProf,:) = MF;
    OutData.Lat( iProf,:) = latmean;
    OutData.Lon( iProf,:) = lonmean;
    OutData.Alt( iProf,:) = Data.Alt(iProf,:);
    OutData.Time(iProf,:) = Data.Time(iProf,:);
    OutData.Tp(  iProf,:) = Data.Tp(iProf,:);


    %store the new ST for the next pass
    NextST = ThisST; 

    clear CoSpectrum A idx Lz dx dPhi Lh latmean lonmean MF dt Bad

  end



 if Settings.Verbose == 1; if mod(iProf,200); textprogressbar(100.*(NProfiles-iProf)./NProfiles); end; end
end; clear iProf NextST ThisST NextST
if Settings.Verbose == 1; textprogressbar(100); textprogressbar('!'); end


%apply mask
f= fieldnames(OutData);
for iF=1:1:numel(f); F = OutData.(f{iF}); F(Mask == 0) = NaN; OutData.(f{iF}) = F; end; 
OutData.FailReason(Mask == 0) = 4;
clear iF f F Mask
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% planetary wave filter, grid-based approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data,PW] = filter_pwgrid(Data,Settings)


%compute a time grid to work on, and also store metadata about the PWs
PW.Time = min(Settings.TimeRange):Settings.PWTimeRes:max(Settings.TimeRange);
PW.Lon = Settings.PWLonGrid; PW.Lat = Settings.PWLatGrid; PW.Alt = Settings.PWAltGrid; PW.PWs = 0:1:Settings.NPWs;
PW.WindowSize = Settings.PWWindow; PW.MinPercent = Settings.PWMinPC;

%generate a store array for the PWs, and for the output T'
PW.PW = NaN(numel(Settings.PWLonGrid),numel(Settings.PWLatGrid),numel(Settings.PWAltGrid),Settings.NPWs+1,numel(PW.Time));
A = Data.Tp.*NaN; %working variable used internally to simplify logic

%fill it, stepping over day-by-day using a time window as specified
if Settings.Verbose == 1; textprogressbar('--> Computing planetary waves '); end

for iStep=1:1:numel(PW.Time)

  %select the data we need by finding find the indices of the points in the UseWindow (points to compute from) and
  %the Output window (points to store, higher resolution)
  %***logic assumes idxO is a subset of idxU***
  idxU = inrange(Data.Time,PW.Time(iStep)+[-1,1].*(Settings.PWWindow)./2);
  idxO = inrange(Data.Time,PW.Time(iStep)+[-1,1].*(mean(diff(PW.Time)))./2);
  if numel(idxU) == 0 || numel(idxO) == 0; continue; end

  %reduce data down to just what we want to use
  PWCalcData = reduce_struct(Data,idxU,[],0);

  %compute the PWs in the use window, and store it in placeholder A
  %A will overwrite most loops - this is fine as long as idxO is a subset of idxU
  [A(idxU),b] = pwfilter(Settings.NPWs,Settings.PWMinPC,                      ...
                         PWCalcData.Lon,PWCalcData.Lat,PWCalcData.Tp,   ...
                         Settings.PWLonGrid,Settings.PWLatGrid,               ...
                         PWCalcData.Alt,Settings.PWAltGrid);
  %store the Tp and PW data
  Data.Tp(idxO) = A(idxO);
  PW.PW(:,:,:,:,iStep) = permute(b,[2,1,3,4]);

  if Settings.Verbose == 1; textprogressbar(100.*iStep./numel(PW.Time)); end

end; clear iDay OutWindow UseWindow idxU PWCalcData a b idxO A iStep
if Settings.Verbose == 1; textprogressbar(100); textprogressbar('!');end

return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D Savitzky-Golay filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Data = filter_sgolay(Data,Settings)

%is the data all the same resolution?
if nanstd(flatten(diff(Data.Alt,1,2)))./nanmean(flatten(diff(Data.Alt,1,2))) < 0.01;

  %if so, we can do this in a single pass
  FrameLen = abs(make_odd(round(Settings.SGLength./nanmean(flatten(diff(Data.Alt,1,2))))));
  Data.Tp = Data.Temp-sgolayfilt(Data.Tp',Settings.SGOrder,FrameLen)';
else

  %if not, we need a loop
  for iProf=1:1:numel(Data.Alt,1)
    FrameLen = abs(make_odd(round(Settings.SGLength./nanmean(Data.Alt(iProf,:)))));
    Data.Tp(iProf,:) = sgolayfilt(Data.Tp(iProf,:),Settings.SGOrder,FrameLen);
  end
end

if Settings.Verbose == 1; disp(['--> Savitzky-Golay filter applied vertically']); end

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hindley23 PW filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data,PW] = filter_hindley23(Data,Settings)


 %copy the T' and T fields generated by Neil's code
 Data.Tp   = Data.Temp_Residual;
 Data.Temp = nansum(Data.Temp_PW,3);

 %store PWs
 PW = Data.Temp_PW;

 %remove outliers
 if Settings.H23_OutRem == true

   %latitude and longitude have physical limits
   Bad = [];
   Bad = [Bad;find(Data.Lat <  -90 | Data.Lat >  90 | Data.Lon < -180 | Data.Lon > 180)];

   %temperature should be >100K always, and  <400K at altitudes below the mesopause
   Bad = [Bad;find(Data.Temp < 100)];
   Bad = [Bad;find(Data.Temp > 400 & Data.Alt < 100)];

   %do it
   Bad = unique(Bad);
   Fields = fieldnames(Data);
   for iF=1:1:numel(Fields)
     F = Data.(Fields{iF});
     F(Bad) = NaN;
     Data.(Fields{iF}) = F;
   end
 end


return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ensure data is on a regular Z grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Data = regularise_data_z(Data,Settings)


%first, check if the data is ALREADY regular AND ascending in z. Counts if the full distribution is within 5% of the mean.
%if it's fine, we don't need to proceed
dZ_distrib = unique(diff(Data.Alt,1,2));


if min(diff(Data.Alt,1,2),[],'all') > 0 & range(dZ_distrib) < 0.1.* nanmean(dZ_distrib); return; end


%ok, we need to regularise. First, work out a scale
NewZ = nanmin(Data.Alt(:)):nanmean(abs(dZ_distrib)):nanmax(Data.Alt(:));
clear dZ_distrib


%create struct to store new data
NewData = struct();
Fields  = fieldnames(Data); Fields(ismember(Fields,'Alt')) = [];
for iF=1:1:numel(Fields); NewData.(Fields{iF}) = NaN([size(Data.Alt,1),numel(NewZ)]); end

%now interpolate EVERYTHING onto the new scale
Fields  = fieldnames(Data); Fields(ismember(Fields,'Alt')) = [];
if Settings.Verbose == 1; disp('--> Data not on a regular and common height grid, interpolating to [ min Z: mean dZ : max Z ] '); end

for iF=1:1:numel(Fields)

  %get data fields
  Fin  =    Data.(Fields{iF});
  Fout = NewData.(Fields{iF}); 

  %check we're working on a 2D profiles x heights array, and just pass it straight through if not
  %do this by chekcing if it's the same size as the 'Alt' array (arbitrary choice)
  if ~isequal(size(Fin),size(Data.Alt)); NewData.(Fields{iF}) = Fin; continue; end


  %interpolate profile-by-profile
  for iProfile=1:1:size(Data.Alt,1);
    Good = find(~isnan(Data.Alt(iProfile,:)));
    if numel(Good) < 2; continue; end
    Fout(iProfile,:) = interp1(Data.Alt(iProfile,Good),Fin(iProfile,Good),NewZ);
  end

  %store data
  NewData.(Fields{iF}) = Fout;
end

%store new altitudes
NewData.Alt = repmat(NewZ,size(Data.Alt,1),1);

%copy over, and return
Data = NewData;





return
end

