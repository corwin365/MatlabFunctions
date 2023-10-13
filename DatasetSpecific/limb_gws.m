function [OutData,PW] = limb_gws(Data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%general-purpose function to find GWs from limb sounders
%
%the routine assumes when using a vertical filter that singularities 
%like the tropopause have been removed - make sure this is the case.
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/09/04
%
%inputs:
%  required:
%    InstrumentData [struct] - data to use, in the format produced by get_limbsounders.m
%
%  optional:
%
%    VarName         (type,           default)  description
%    -----------------------------------------------------------------------------
%    Filter          (char,          'PWgrid')  type of detrending filter to use (see below)
%    STScales        (vector,   1:1:NLevels/2)  number of scales to use in 1D ST
%    STc             (positive real,     0.25)  value of 'c' to use in ST
%    Analysis        (integer,              1)  type of analysis to use, from (1) 1DST, (2) Alexander et al 2008. 
%    RegulariseZ     (logical            true)  ensure the data is on a regular height grid
%
%
%if Filter is set to 'Hindley23_ISSI', then the input 'Data' struct must contain
%the following fields:
%
%    H23_ISSI_Source (char,             'Obs')  subdataset to use, i.e. 'Model' or Obs'
%    H23_ISSI_OutRem (logical           true)   remove outliers from H23 fitting, as in get_limbsounders()
%
%if Filter is set to 'PWgrid', then the following options can be used:
%
%    NPWs            (integer,              3)  number of PWs to fit (plus zonal mean)
%    PWWindow        (positive real,        1)  number of days to use in each PW fit
%    PWTimeRes       (positive real,        1)  time resolution to export PW fits on
%    PWLonGrid       (vector      -180:20:180)  longitude grid to compute PWs on
%    PWMinPC         (positive real,       66)  fraction of longitude bins which must be filled for each lat band
%    PWLatGrid       (vector         -90:5:90)  longitude grid to compute PWs on
%    PWAltGrid       (vector, input grid mean)  altitude grid to compute PWs on
%
%if Filter is set to 'SGolay', then the following options can be used:
%
%    SGOrder         (integer,              2)  order of Savitzky-Golay filter to use
%    SGLength        (real,                18)  frame length of SG filter, in km
%
%
%outputs:
%   Data: struct containing all variables, on a [profiles x height] grid
%   PWs:  array containing PWs, on a lon x lat x height x wave x days grid
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
addParameter(p,'Analysis',1,@(x) validateattributes(x,{'numeric'},{'integer','<=',2})) %type of analysis to use (see above)

%filter to use
addParameter(p,'Filter','PWgrid',@ischar) %type of filter to use

%ST properties
addParameter(p,'STScales',1:1:size(Data.Alt,2)/2,@isvector  ) %scales to compute on ST
addParameter(p,'STc',                       0.25,@ispositive) %'c' parameter for ST

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
addParameter(p,'H23_ISSI_Source', 'Obs',@ischar) %obs or model?
addParameter(p,'H23_ISSI_OutRem',  true,@islogical) %remove outliers?

%parse inputs and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,Data,varargin{:})

%check contents of input data struct:. 
  %%do we have at least Lat, Lon, Temp, Alt? These are what we use
if ~isfield(Data,'Lat') ||  ~isfield(Data,'Lon') || ~isfield(Data,'Alt') || ~isfield(Data,'Temp');
  error('Missing field - struct must contain Lat, Lon, Alt and Temp fields.')
  return
end
  %%are these fields all the same size?
if ~isequal(size(Data.Lat),size(Data.Lon)) ||  ~isequal(size(Data.Lat),size(Data.Alt)) || ~isequal(size(Data.Lat),size(Data.Temp));
  error('LAt, Lon, Alt and Temp fields must all be the same size.')
  return
end


%pull out the remaining arguments into struct "Inputs", used throughout rest of routine
Input = p.Results;
Data = Input.Data; Input = rmfield(Input,'Data');
Input.TimeRange = [floor(nanmin(Data.Time,[],'all')),ceil(nanmax(Data.Time,[],'all'))];
clearvars -except InstInfo Input Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ensure the data is on a regular height grid
% for Hindley23 removal, this must be done AFTER detrending, as new data are loaded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Input.RegulariseZ == true && ~strcmpi(Input.Filter,'Hindley23_ISSI'); Data = regularise_data_z(Data,Input); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% detrend the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start by duplicating temperatures. All filter functions will act on this,
%in order to let us combine filters if we want
Data.Tp = Data.Temp;

%now, apply the filter. Currently only one at a time, but simple to rewrite
%this to apply multiple filters if needed

if     strcmpi(Input.Filter,        'PWgrid'); [Data,PW] = filter_pwgrid(       Data,Input);
elseif strcmpi(Input.Filter,        'SGolay'); Data      = filter_sgolay(       Data,Input);
elseif strcmpi(Input.Filter,'Hindley23_ISSI'); [Data,PW] = filter_hindley23issi(Data,Input); 
else
  disp('Filter type not included in programme, stopping')
  stop
end


%create empty PW output if we didn't generate one. Note that this output will be 
%fairly meaningless anyway if the PW filter wasn't the first-applied
if ~exist('PW','var'); PW.Comment = 'Planetary wave filter not used, no PW data computed'; end


%regularisation if using Hindley23 detrending
if Input.RegulariseZ == true && strcmpi(Input.Filter,'Hindley23_ISSI'); Data = regularise_data_z(Data,Input); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S-Transform profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%produce storage arrays
NProfiles = size(Data.Tp,1);
NLevs     = size(Data.Tp,2);
OutData   = spawn_uniform_struct({'A','Lz','Lh','Lat','Lon','Alt','Tp','MF'},[NProfiles,NLevs]);
Mask      = ones([NProfiles,NLevs]); %this is used to mask out bad data later

%some approaches require two adjacent profiles to be computed. To avoid duplicate computation,
%this is marginally more efficient if we work backwards and store the next profile for this case.


textprogressbar('--> Computing gravity waves ')
for iProf=NProfiles:-1:1


  %fill any NaNs with a zero for purposes of the ST
  %replace with NaNs afterwards
  NoData = find(isnan(Data.Tp(iProf,:)));
  Mask(iProf,NoData) = NaN;
  Tp = Data.Tp(iProf,:); Tp(NoData) = 0;

  %compute ST
  ThisST = nph_ndst(Tp,                  ...
                    Input.STScales,                    ...
                    nanmean(diff(Data.Alt(iProf,:))),  ...
                    Input.STc);
  clear Tp NoData

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %compute and store results depending on analysis
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if Input.Analysis == 1

    %simple 1DST approach - just store
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OutData.A(  iProf,:) = ThisST.A;
    OutData.Lz( iProf,:) = 1./ThisST.F1;
    OutData.Lat(iProf,:) = Data.Lat(iProf,:); OutData.Lon(iProf,:) = Data.Lon(iProf,:);
    OutData.Alt(iProf,:) = Data.Alt(iProf,:); OutData.Tp( iProf,:) = Data.Tp( iProf,:);

  elseif Input.Analysis == 2

    %Alexander et al (JGR, 2008) approach
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %if this is the first profile to be processed, retain the data then move to the next
    if iProf == NProfiles; NextST = ThisST; continue; end

    %compute cospectrum
    CoSpectrum = ThisST.ST .* conj(NextST.ST);

    %find maximum amplitude and wavelength at each level
    [A,idx] = max(abs(CoSpectrum),[],1);
    A = sqrt(A);
    Lz= 1./ThisST.freqs(idx);

    %find along-track distance between profiles, phase difference, and hence Lh
    dx   = nph_haversine([Data.Lat(iProf,  :);Data.Lon(iProf,  :)]',[Data.Lat(iProf+1,:);Data.Lon(iProf+1,:)]');
    dPhi = NaN(size(Lz)); for iLev=1:1:size(Data.Alt,2); dPhi(iLev) = CoSpectrum(idx(iLev),iLev); end
    dPhi = atan(imag(dPhi)./real(dPhi))./(2.*pi);
    Lh = abs(dx./dPhi');

    %compute MF
    MF = cjw_airdensity(Data.Pres(iProf,:),Data.Temp(iProf,:))  ...
      .* (Lz ./ Lh')                                            ...
      .* (9.81./0.02).^2                                        ...
      .*A./(Data.Temp(iProf,:)-Data.Tp(iProf,:)).^2;


    %adjust lat/lon to the midpoint of the profile-pair
    [latmean,lonmean] = meanm(Data.Lat(iProf+[0,1],:), ...
                              Data.Lon(iProf+[0,1],:));

    %store results
    OutData.Lat(iProf,:) = latmean;
    OutData.Lon(iProf,:) = lonmean;
    OutData.A(  iProf,:) = A;
    OutData.Lz( iProf,:) = Lz;
    OutData.Lh( iProf,:) = Lh;
    OutData.Alt(iProf,:) = Data.Alt(iProf,:);
    OutData.Tp( iProf,:) = Data.Tp(iProf,:);
    OutData.MF( iProf,:) = MF;

    %store the new ST for the next pass
    NextST = ThisST; 

    clear CoSpectrum A idx Lz dx dPhi Lh latmean lonmean MF

  end



  % % % if numel(NoData) ~= 0;
  % % %   F = fieldnames(OutData);
  % % %   for iF=1:1:numel(F)
  % % %     FieldData = OutData.(F{iF});
  % % %     FieldData(iProf,NoData) = NaN;
  % % %     OutData.(F{iF}) = FieldData;
  % % %   end
  % % %   clear F iField
  % % % end

 if mod(iProf,100); textprogressbar(100.*(NProfiles-iProf)./NProfiles); end
end; clear iProf NextST ThisST NextST
textprogressbar(100); textprogressbar('!')

%remove the data where we put zeros
if sum(Mask,'all') ~= numel(Mask);
  f = fieldnames(OutData);
  for iF=1:1:numel(f)
    OutData.(f{iF}) = OutData.(f{iF}).*Mask;
  end; clear iF f
end

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% planetary wave filter, grid-based approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data,PW] = filter_pwgrid(Data,Input)


%compute a time grid to work on, and also store metadata about the PWs
PW.Time = min(Input.TimeRange):Input.PWTimeRes:max(Input.TimeRange);
PW.Lon = Input.PWLonGrid; PW.Lat = Input.PWLatGrid; PW.Alt = Input.PWAltGrid; PW.PWs = 0:1:Input.NPWs;
PW.WindowSize = Input.PWWindow; PW.MinPercent = Input.PWMinPC;

%generate a store array for the PWs, and for the output T'
PW.PW = NaN(numel(Input.PWLonGrid),numel(Input.PWLatGrid),numel(Input.PWAltGrid),Input.NPWs+1,numel(PW.Time));
A = Data.Tp.*NaN; %working variable used internally to simplify logic

%fill it, stepping over day-by-day using a time window as specified
textprogressbar('--> Computing planetary waves ')

for iStep=1:1:numel(PW.Time)

  %select the data we need by finding find the indices of the points in the UseWindow (points to compute from) and
  %the Output window (points to store, higher resolution)
  %***logic assumes idxO is a subset of idxU***
  idxU = inrange(Data.Time,PW.Time(iStep)+[-1,1].*(Input.PWWindow)./2);
  idxO = inrange(Data.Time,PW.Time(iStep)+[-1,1].*(mean(diff(PW.Time)))./2);
  if numel(idxU) == 0 || numel(idxO) == 0; continue; end

  %reduce data down to just what we want to use
  PWCalcData = reduce_struct(Data,idxU,[],0);

  %compute the PWs in the use window, and store it in placeholder A
  %A will overwrite most loops - this is fine as long as idxO is a subset of idxU
  [A(idxU),b] = pwfilter(Input.NPWs,Input.PWMinPC,                      ...
                         PWCalcData.Lon,PWCalcData.Lat,PWCalcData.Tp,   ...
                         Input.PWLonGrid,Input.PWLatGrid,               ...
                         PWCalcData.Alt,Input.PWAltGrid);
  %store the Tp and PW data
  Data.Tp(idxO) = A(idxO);
  PW.PW(:,:,:,:,iStep) = permute(b,[2,1,3,4]);

  textprogressbar(100.*iStep./numel(PW.Time))

end; clear iDay OutWindow UseWindow idxU PWCalcData a b idxO A iStep
textprogressbar(100); textprogressbar('!')

return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D Savitzky-Golay filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Data = filter_sgolay(Data,Input)

%is the data all the same resolution?
if nanstd(flatten(diff(Data.Alt,1,2)))./nanmean(flatten(diff(Data.Alt,1,2))) < 0.01;

  %if so, we can do this in a single pass
  FrameLen = abs(make_odd(round(Input.SGLength./nanmean(flatten(diff(Data.Alt,1,2))))));
  Data.Tp = Data.Temp-sgolayfilt(Data.Tp',Input.SGOrder,FrameLen)';
else

  %if not, we need a loop
  for iProf=1:1:numel(Data.Alt,1)
    FrameLen = abs(make_odd(round(Input.SGLength./nanmean(Data.Alt(iProf,:)))));
    Data.Tp(iProf,:) = sgolayfilt(Data.Tp(iProf,:),Input.SGOrder,FrameLen);
  end
end

disp(['--> Savitzky-Golay filter applied vertically'])

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hindley23 PW filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data,PW] = filter_hindley23issi(Data,Input)

 %copy the T' and T fields Neil generated over
 if strcmpi(Input.H23_ISSI_Source,'Obs') | strcmpi(Input.H23_ISSI_Source,'Model');
   Data.Tp   = Data.(Input.H23_ISSI_Source).Temp_Residual;
   Data.Temp = Data.(Input.H23_ISSI_Source).Temp;
 else
   error("This data source is not available for the Hindley23 filter method, must be one of 'Obs' or 'Model");
 end

 %store PWs
 PW = Data.(Input.H23_ISSI_Source).Temp_PW;
 Data = rmfield(Data,{'Obs','Model'}); %they cause trouble in the logic, and we don't need them any mores

 %remove outliers
 if Input.H23_ISSI_OutRem == true

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

function Data = regularise_data_z(Data,Input)


%first, check if the data is ALREADY regular in z. Counts if the full distribution is within 5% of the mean.
%if it's fine, we don't need to proceed
dZ_distrib = unique(diff(Data.Alt,1,2));
if range(dZ_distrib) < 0.1.* nanmean(dZ_distrib); return; end


%ok, we need to regularise. First, work out a scale
NewZ = nanmin(Data.Alt(:)):nanmean(dZ_distrib):nanmax(Data.Alt(:));
clear dZ_distrib


%create struct to store new data
NewData = struct();
Fields  = fieldnames(Data); Fields(ismember(Fields,'Alt')) = [];
for iF=1:1:numel(Fields); NewData.(Fields{iF}) = NaN([size(Data.Alt,1),numel(NewZ)]); end

%now interpolate EVERYTHING onto the new scale
Fields  = fieldnames(Data); Fields(ismember(Fields,'Alt')) = [];
disp('--> Data not on a regular and common height grid, interpolating to [ min Z: mean dZ : max Z ] ')

for iF=1:1:numel(Fields)

  %get data fields
  Fin  =    Data.(Fields{iF});
  Fout = NewData.(Fields{iF}); 

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

