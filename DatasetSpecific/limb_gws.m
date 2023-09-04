function [OutData,PWs] = limb_gws(Data,varargin)

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
%  optional:
%
%    VarName         (type,           default)  description
%    -----------------------------------------------------------------------------
%    Filter          (char,              'PW')  type of detrending filter to use
%    STScales        (vector,   1:1:NLevels/2)  number of scales to use in 1D ST
%    STc             (positive real,     0.25)  value of 'c' to use in ST
%    Analysis        (integer,              1)  type of analysis to use, from (1) 1DST, (2) Alexander et al 2008. 
%
%if Filter is set to 'PW', then the following options can be used:
%
%    NPWs            (integer,              6)  number of PWs to fit (plus zonal mean)
%    PWWindow        (positive real,        1)  number of days to use in each PW fit
%    PWTimeRes       (positive real,        1)  time resolution to export PW fits on
%    PWLonGrid       (vector      -180:20:180)  longitude grid to compute PWs on
%    PWMinPC         (positive real,     0.66)  fraction of longitude bins which must be filled for each lat band
%    PWLatGrid       (vector         -90:5:90)  longitude grid to compute PWs on
%    PWAltGrid       (vector, input grid mean)  altitude grid to compute PWs on
%
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
CheckAnalysis = @(x) validateattributes(x,{'numeric'},{'integer','<=',2}); 
addParameter(p,'Analysis', 1,CheckAnalysis) %type of analysis to use (see above)

%filter to use
addParameter(p,'Filter','PW',@ischar) %type of filter to use

%ST properties
addParameter(p,'STScales',1:1:size(Data.Alt,2)/2,@isvector  ) %scales to compute on ST
addParameter(p,'STc',                       0.25,@ispositive) %'c' parameter for ST


%additional variables - planetary waves
addParameter(p,      'NPWs',                  6,@isinteger ) %number of planetary waves to fit
addParameter(p,   'PWMinPC',               0.66,@ispositive) %fraction of longitude bins that must be filled
addParameter(p,  'PWWindow',                  1,@ispositive) %window width of period used
addParameter(p, 'PWTimeRes',                  1,@ispositive) %time resolution
addParameter(p, 'PWLonGrid',        -180:20:180,@isvector  ) %longitude bins
addParameter(p, 'PWLatGrid',           -90:5:90,@isvector  ) %latitude bins
addParameter(p, 'PWAltGrid',nanmean(Data.Alt,1),@isvector  ) %altitude bins



%parse inputs and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,Data,varargin{:})

%check contents of input data struct. 
disp('Currently no checks on contents of Data struct, fix this!')
disp('need to check that altitudes monotonically increase and are non-NaN')


%pull out the remaining arguments into struct "Inputs", used throughout rest of routine
Input = p.Results;
Data = Input.Data; Input = rmfield(Input,'Data');
Input.TimeRange = [floor(nanmin(Data.Time,[],'all')),ceil(nanmax(Data.Time,[],'all'))];
clearvars -except InstInfo Input Data








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute temperature perturbations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(Input.Filter,'PW')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% planetary wave filter
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %compute a time grid to work on, and also store metadata about the PWs
  PW.Time = min(Input.TimeRange):Input.PWTimeRes:max(Input.TimeRange);
  PW.Lon = Input.PWLonGrid; PW.Lat = Input.PWLatGrid; PW.Alt = Input.PWAltGrid; PW.PWs = 0:1:Input.NPWs;
  PW.WindowSize = Input.PWWindow; PW.MinPercent = Input.PWMinPC;

  %generate a store array for the PWs, and for the output T'
  PW.PW = NaN(numel(Input.PWLonGrid),numel(Input.PWLatGrid),numel(Input.PWAltGrid),Input.NPWs+1,numel(PW.Time));
  Data.Tp = Data.Temp.*NaN;
  A = Data.Tp.*NaN; %working directory used internally to simplify logic

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
    [A(idxU),b] = pwfilter(Input.NPWs,Input.PWMinPC,                        ...
                           PWCalcData.Lon,PWCalcData.Lat,PWCalcData.Temp, ...
                           Input.PWLonGrid,Input.PWLatGrid,                   ...
                           PWCalcData.Alt,Input.PWAltGrid);
    %store the Tp and PW data
    Data.Tp(idxO) = A(idxO);
    PW.PW(:,:,:,:,iStep) = permute(b,[2,1,3,4]);

    textprogressbar(100.*iStep./numel(PW.Time))

  end; clear iDay OutWindow UseWindow idxU PWCalcData a b idxO A iStep
  textprogressbar(100); textprogressbar('!')

else

  disp('Filter type not included in programme, stopping')
  stop

end

%create empty PW output if we use a filter other than a PW
if ~strcmp(Input.Filter,'PW'); PW.Comment = 'Planetary wave filter not used, no PW data computed'; end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S-Transform profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%produce storage arrays
NProfiles = size(Data.Tp,1);
NLevs     = size(Data.Tp,2);

OutData.A      = NaN(NProfiles,NLevs);
OutData.Lz     = OutData.A;
OutData.Lh     = OutData.A;
OutData.Lat    = OutData.A;
OutData.Lon    = OutData.A;
OutData.Alt    = OutData.A;
OutData.Tp     = OutData.A;


%some approaches require two adjacent profiles to be computed. To avoid duplicate computation,
%this is marginally more efficient if we work backwards and store the next profile for this case.
textprogressbar('--> Computing   gravity waves ')
for iProf=NProfiles:-1:1

  %compute ST
  ThisST = nph_ndst(Data.Tp(iProf,:),                  ...
                    Input.STScales,                    ...
                    nanmean(diff(Data.Alt(iProf,:))),  ...
                    Input.STc);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %compute and store results depending on analysis
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if Input.Analysis == 1;

    %simple 1DST approach - just store
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OutData.A(  iProf,:) = ThisST.A;
    OutData.Lz( iProf,:) = 1./ThisST.F1;
    OutData.Lat(iProf,:) = Data.Lat(iProf,:);
    OutData.Lon(iProf,:) = Data.Lon(iProf,:);
    OutData.Alt(iProf,:) = Data.Alt(iProf,:);
    OutData.Tp( iProf,:) = Data.Tp(iProf,:);

  elseif Input.Analysis == 2;

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

    %store the new ST for the next pass
    NextST = ThisST; 

    clear CoSpectrum A idx Lz dx dPhi Lh latmean lonmean

  end





 textprogressbar(100.*(NProfiles-iProf)./NProfiles)
end; clear iProf NextST ThisST 
textprogressbar(100); textprogressbar('!')



end
