function [ST,Airs,Error,ErrorInfo] = gwanalyse_airs_3d(Airs,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%generalised function to ST AIRS data prepared by prep_airs_3d
%
%Corwin Wright, c.wright@bath.ac.uk, 10/AUG/2019
%extensively modified 07/JUL/2020 to add "2D+1" option, using phase differences to compute Lz
%
%IMPORTANT: this routine does not allow most of the optional flags of
%nph_ndst (e.g. Scales, Full Mode) - only the NUMBER of scales and c can be
%controlled. To use the more advanced features, use the NDST routine
%directly.
%
%ASSUMES DATA ARE ALREADY REGULARLY SPACED, OR AT LEAST WILL BE AFTER
%HEIGHT SUBSETTING. 
%
%REQUIRES the following external functions:
% nph_ndst
% gwanalyse_airs_2d (only if we request 2D+1 output)
%
%INCLUDES the following functions written by others:
% Neil Hindley: nph_haversine,scale_height (n.hindley@bath.ac.uk)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%
%outputs:
%%%%%%%%%%%
%
%  ST:  STed AIRS data, provided no critical errors hit. Empty struct otherwise.
%  Airs: input AIRS data with small modifcations (e.g. same height range as ST output)
%
%  Error: error state; 0: no error; 1: critical error; 2: minor error
%  ErrorInfo: critical error details 
%
%%%%%%%%%%%
%inputs:
%%%%%%%%%%%
%
 %required
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
%    VarName (type)   description
%    -----------------------------------------------------------------------------
%    Airs    (struct) Airs data, usually as output from prep_airs_3d
%
%
 %optional (all case-insensitive)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    VarName         (type,           default)  description
%    -----------------------------------------------------------------------------
%    NScales         (numeric,           1000)  number of frequency scales to use
%    ZRange          (array,          [20,60])  height range of data to use
%    HeightScaling   (logical,           true)  scale data for height before ST (undone afterwards)
%    Spacing         (array,    [NaN,NaN,NaN])  element spacing computed externally
%    c               (array, [0.25,0.25,0.25])  s-transform c parameter
%    MinWaveLength   (array,   [1,1,1].*99e99)  minimum output wavelength
%    MaxWaveLength   (array,          [0,0,0])  maximum output wavelength
%    NotAirsData     (logical,          false)  overrides some of the sanity checks on AIRS data formatting
%    TwoDPlusOne     (logical,          false)  *ADDITIONALLY* compute vertical wavelengths with 2D+1 method
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set default outputs
ST = struct();
Error = 0;
ErrorInfo = '';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this function is written for AIRS data, but can be used for non-AIRS data
%in this case, some sanity checks below need to be overridden
%result: Status = 0 to apply the checks, otherwise do not
Status = find(strcmp(varargin,'NotAirsData') ~= 0);
if numel(Status) ~= 0; Status = varargin{Status+1}; else Status = 0; end



%create input parser
%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;


%inputs - required
%%%%%%%%%%%%%%%%%%%
addRequired(p,'Airs',  @isstruct);


%inputs - fully optional
%%%%%%%%%%%%%%%%%%%%%%%%%

%NScales must be a positive integer
CheckNScales = @(x) validateattributes(x,{'numeric'},{'nonnegative','integer'});  
addParameter(p,'NScales',1000,CheckNScales);  %assume 1000 scales unless specified  

%ZRange must be an array with two positive values between 0 and 90
if Status == 0; CheckZRange = @(x) validateattributes(x,{'numeric'},{'size',[1,2],'>=',   0,'<=', 90});
else            CheckZRange = @(x) validateattributes(x,{'numeric'},{'size',[1,2],'>=',-Inf,'<=',Inf});
end
addParameter(p,'ZRange',[20,60],CheckZRange);  %assumes only the 3km regular step range

%Spacing and c must be arrays with three positive values 
CheckSpacing = @(x) validateattributes(x,{'numeric'},{'size',[1,3],'positive'});
addParameter(p,'Spacing',[NaN,NaN,NaN],CheckSpacing);  %assumes we want data regularly spaced
addParameter(p,'c',[0.25,0.25,0.25],CheckSpacing);  %default values

%HeightScaling is logical
addParameter(p,'HeightScaling',true,@islogical);  %assumes we want to scale the data with height for STing (put back afterwards)

%NotAirsData is logical
addParameter(p,'NotAirsData',false,@islogical);  %assumes we are feeding the routine AIRS data

%TwoDPlusOne is logical, and its settings are a struct (which should be blank if not set)
addParameter(p,        'TwoDPlusOne',   false,@islogical); 
addParameter(p,'TwoDPlusOneSettings',struct(),@isstruct); 

%MaxWaveLength must be an positive real number
CheckLambda = @(x) validateattributes(x,{'numeric'},{'>=',0});
addParameter(p,'MaxWaveLength',[1,1,1].*99e99,CheckLambda);  %crazy-large default

%MinWaveLength must be an positive real number
addParameter(p,'MinWaveLength',[1,1,1].*0,CheckLambda);  %zero default





%inputs - extra settings for 2D+1 ST
%%%%%%%%%%%%%%%%%%%%%%%%%

%these will come in as a struct, which must exist and is parsed here
%note that the individual flags WILL NOT BE SANITY-CHECKED - caveat emptor...

%first, check if we're doing a 2D call. Then parse out the fields
a = find(strcmp(varargin,'TwoDPlusOne'));
if numel(a) > 0;
  if varargin{a+1} == 1;
    
    %defaults - override if they don't exist in the call
    TwoDSettings.c1          = [1,1].*0.5;               %c for  first pass (coarse in space, fine in freq)
    TwoDSettings.c2          = [1,1].*0.25;              %c for second pass (coarse in freq, fine in space)
    TwoDSettings.NPeaks      = 1;%3;                     %number of spectral peaks to identify for each spatial point
    TwoDSettings.Threshold   = 0;                        %threshold in spectral space to be identified as a peak. This is an amplitude for a lone 2DST (as opposed to a cospectrum)
    TwoDSettings.Filt        = fspecial('gaussian',5,1); %characteristic size of point in spectral space. This is a little fatter than the default, to avoid very close peaks.
    TwoDSettings.Thin        = 1;                        %thin out the number of scales (large runtime reduction, but changes the results)
    TwoDSettings.Steps       = 3;%[1,2,3,4,5];           %number of steps to take phase difference over. '0' takes it from a basis level, defined above, while nonzero values use the phase shift with that many levels *above*
    TwoDSettings.Weight      = 1;                        %height-weight the vertical layers
    
    b = find(strcmp(varargin,'TwoDPlusOneSettings'));
    if numel(b) > 0;
      InStruct = varargin{b+1};  Flags = fieldnames(TwoDSettings);
      for iFlag=1:1:numel(Flags)
        if isfield(InStruct,Flags{iFlag}); TwoDSettings.(Flags{iFlag}) = InStruct.(Flags{iFlag}); end
      end
    end
  end  
end; clear a b iFlag Flags InStruct






%parse the inputs, and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ok, do the parsing!
parse(p,Airs,varargin{:});

%pull out the contents into struct "Inputs", used throughout rest of routine
Input = p.Results;

%tidy up
clearvars -except Input ST Airs Error ErrorInfo TwoDSettings
clear CheckNScales CheckZRange CheckSpacing CheckLambda



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reduce to desired height range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InZRange = find(Airs.ret_z >= Input.ZRange(1) ...
              & Airs.ret_z <= Input.ZRange(2));

            
%if we're using the 2D+1 method, we'l need some extra levels for phase fitting   
%retain these in a separate array - they'll be thrown away after use
if Input.TwoDPlusOne | Input.TwoDPlusOne_ind; 
  if min(InZRange) > 1;
    Extra.Bottom.ret_z    = squeeze(Airs.ret_z(       min(InZRange)-1));
    Extra.Bottom.Tp       = squeeze(Airs.Tp(      :,:,min(InZRange)-1));
    Extra.Bottom.BG       = squeeze(Airs.BG(      :,:,min(InZRange)-1));
    Extra.Bottom.ret_temp = squeeze(Airs.ret_temp(:,:,min(InZRange)-1));
  else
    Extra.Bottom.ret_z = [];
  end
  if max(InZRange) < numel(Airs.ret_z);    
    Extra.Top.ret_z    = squeeze(Airs.ret_z(       min(InZRange)+1));
    Extra.Top.Tp       = squeeze(Airs.Tp(      :,:,min(InZRange)+1));
    Extra.Top.BG       = squeeze(Airs.BG(      :,:,min(InZRange)+1));
    Extra.Top.ret_temp = squeeze(Airs.ret_temp(:,:,min(InZRange)+1));
  else
    Extra.Top.ret_z = [];
  end
end

%Safely stockpiled. Ok, trim down the variables for normal use
Airs.ret_z    = Airs.ret_z(       InZRange);
Airs.Tp       = Airs.Tp(      :,:,InZRange);
Airs.BG       = Airs.BG(      :,:,InZRange);
Airs.ret_temp = Airs.ret_temp(:,:,InZRange);

clear InZRange



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate point spacing, and check it's regular
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isnan(sum(Input.Spacing))
  
  %compute and check XT spacing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Loc1 = [Airs.l1_lat(:,1),Airs.l1_lon(:,1)];
  Loc2 = circshift(Loc1,1,1); Loc2(1,:) = NaN;
  dXT = nph_haversine(Loc1,Loc2);  
  dXT = nanmean(dXT(:));
  
  PointSpacing(1) = nanmean(dXT(:));
  
  %compute and check AT spacing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Loc1 = [Airs.l1_lat(1,:);Airs.l1_lon(1,:)]';
  Loc2 = circshift(Loc1,1,1); Loc2(1,:) = NaN;
  dAT = nph_haversine(Loc1,Loc2);  
  PointSpacing(2) = nanmean(dAT(:));
  
  %compute and check Z spacing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  dZ = diff(Airs.ret_z);
  PointSpacing(3) = nanmean(dZ(:));
  
  clear dAT dXT dZ Loc1 Loc2 Wiggle

else
  PointSpacing = Input.Spacing;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do ST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scale in height, to compensate for wave amplitude growth
if Input.HeightScaling;
  %find scale height at reference level
  
  NormHeight = 42; %shouldn't matter at all as the caling is uniform, so let's leave this hardcoded for once
  H = scale_height(NormHeight);
  
  %and hence compute correction factor for each height
  CFac = exp(-(Airs.ret_z-NormHeight) ./ (2*H));  
  
  %repmat out to full size of data
  sz = size(Airs.Tp);
  CFac = repmat(permute(CFac,[2,3,1]),sz(1),sz(2),1);

  %apply it to the TPrime data
  Airs.Tp = Airs.Tp .* CFac;

  clear NormHeight H sz
end

% do the ST
ST = nph_ndst(Airs.Tp,                                ...
              Input.NScales,                          ...
              PointSpacing,                           ...
              Input.c,                                ...
              'maxwavelengths', Input.MaxWaveLength , ...
              'minwavelengths', Input.MinWaveLength);
clear PointSpacing


if Input.HeightScaling;
  Airs.Tp = Airs.Tp ./ CFac;
  ST.A    = ST.A    ./ CFac;
  ST.HA   = ST.HA   ./ CFac;
  ST.R    = ST.R    ./ CFac;
  ST.HR   = ST.HR   ./ CFac;
  
  
  clear CFac
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if requested, also use the 2D+1 method to compute vertical wavelength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if     Input.TwoDPlusOne;    
  
  %standard version - uses fitted horizontal wavelengths from 3DST
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %settings for 2D+1. consider moving these to primary inputs after testing.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  NewFields = get_2dp1_lambdaz(Airs,Extra,TwoDSettings);
  clear TwoDSettings
  
  
  %add those fields to the main ST array
  Fields = fieldnames(NewFields);
  for iF=1:1:numel(Fields);
    ST.([Fields{iF},'_2dp1']) = NewFields.(Fields{iF});
  end; clear iF Fields NewFields

end
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ok. geometry time. convert from F1, F2 (,F3)  to k, l (,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% COMPUTE ANGLES AND PROJECT:
ST.kh = quadadd(ST.F1,ST.F2);
sz = size(ST.A);

% find angle CLOCKWISE from along track direction...
ang_at = atan2d(ST.F1,ST.F2);

% find azimuth of AT direction CLOCKWISE from north...
xt_mid = floor(sz(1)./2);
[~,az_at] = distance(Airs.l1_lat(xt_mid,1:end-1,1),Airs.l1_lon(xt_mid,1:end-1,1),Airs.l1_lat(xt_mid,2:end,1),Airs.l1_lon(xt_mid,2:end,1));
az_at(end+1) = az_at(end);
az_at = repmat(az_at,size(ST.kh,1),1,size(ST.kh,3));


% add angles... (since they both should be clockwise from north)
ang_north = wrapTo180(az_at + ang_at);

% finally, find k and l:
ST.k = ST.kh .* sind(ang_north);
ST.l = ST.kh .* cosd(ang_north);

%finally, m 
ST.m = ST.F3;


clear sz ang_at az_at ang_north xt_mid


%done!
return
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nph_haversine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [km] = nph_haversine(loc1,loc2)

loc1 = deg2rad(loc1); loc2 = deg2rad(loc2);

R = 6371;                                 % Earth's radius in km
delta_lat = loc2(:,1) - loc1(:,1);        % difference in latitude
delta_lon = loc2(:,2) - loc1(:,2);        % difference in longitude
a = sin(delta_lat./2).^2 + cos(loc1(:,1)) .* cos(loc2(:,1)) .* ...
    sin(delta_lon./2).^2;
c = 2 .* atan2(sqrt(a), sqrt(1-a));
km = R .* c;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scale_height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function H = scale_height(alt)

% define pressure "base" levels:

p0 = 1013.25; % (hPa) (sea level)
p1 = 226.3210; 
p2 = 54.7489;
p3 = 8.6802;
p4 = 1.1091;
p5 = 0.6694;
p6 = 0.0396;

% and their corresponding altitude levels (km)

z0 = 0;
z1 = 11;
z2 = 20;
z3 = 32;
z4 = 47;
z5 = 51;
z6 = 71;

% now calculate scale heights...

H0 = - (z1-z0) / (ln(p1/p0));
H1 = - (z2-z1) / (ln(p2/p1));
H2 = - (z3-z2) / (ln(p3/p2));
H3 = - (z4-z3) / (ln(p4/p3));
H4 = - (z5-z4) / (ln(p5/p4));
H5 = - (z6-z5) / (ln(p6/p5));

H6 = H5;



if alt <= z1, H = H0; end
if alt > z1 && alt <= z2, H = H1; end
if alt > z2 && alt <= z3, H = H2; end
if alt > z3 && alt <= z4, H = H3; end
if alt > z4 && alt <= z5, H = H4; end
if alt > z5 && alt <= z6, H = H5; end
if alt > z6, H = H6; end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ln
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Out = ln(In)

Out = log(In);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% quadadd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c] = quadadd(a,b)
%Pythagorian addition of the two terms

c = sqrt(a.^2 + b.^2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D+1 st, version 1. This computes independent horizontal wavenumber fits.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NewFields = get_2dp1_lambdaz(Airs,Extra,Settings)

    
  %glue the extra levels we retained to the top and bottom
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Vars = {'Tp','ret_z','BG','ret_temp'};
  Dims = [3,1,3,3];
  for iVar=1:1:numel(Vars)
    Working.(Vars{iVar}) = Airs.(Vars{iVar});
    if isfield(Extra,'Bottom');
      if numel(Extra.Bottom.ret_z) > 0;
        Working.(Vars{iVar}) = cat(Dims(iVar),Extra.Bottom.(Vars{iVar}),Working.(Vars{iVar}));
      end;
    end
    if isfield(Extra,'Top');
      if numel(Extra.Top.ret_z) > 0;
        Working.(Vars{iVar}) = cat(Dims(iVar),Working.(Vars{iVar}),Extra.Top.(Vars{iVar}));
      end;
    end;
    Airs.(Vars{iVar}) = Working.(Vars{iVar});
  end; clear iVar Working Vars Dims

  
  %% 2DST every level, and store as a height-weighted sum
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %define scales
  Scales{1} = 1:1:size(Airs.l1_lat,1)/2;
  Scales{2} = 1:1:size(Airs.l1_lat,2)/2; Scales{2} = [reverse(-Scales{2}),Scales{2}];
  
  %thin scales?
  if Settings.Thin
    
    %some coarser scales that should still resolve wave structure well
    %these numbers assume a standard-sized AIRS granule
    lx = unique(floor((20.* 90)./logspace(log10(20),log10( 20.*90),30))); %UP TO 30 freqs xt
    ly = unique(floor((18.*135)./logspace(log10(18),log10(18.*135),30))); %UP TO 30 freqs at    
    
    Scales{1} = lx; Scales{2} = ly;
    clear lx ly
    
% % %     s1 = Scales{1}; Scales{1} = s1([1:1:5,6:2:end]);
% % %     s2 = Scales{2}; Scales{2} = s2([1:1:5,6:2:end]);
% % %     clear s1 s2
  end

  textprogressbar('Initial ST pass ')
  for iLevel=1:1:numel(Airs.ret_z);
    
    %st level
    ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel), ...
                           'FullST',true,           ...
                           'c',Settings.c1,         ...
                           'Scales',Scales);
    
    %compute weight
    if Settings.Weight == 1; CFac = exp(-(Airs.ret_z(iLevel)-42) ./ (2*scale_height(42))); %42 is arbitrary, but must be consistent across levels
    else                     CFac = 1;
    end
    ST.ST  = ST.ST .* CFac;
    
    %store
    if iLevel == 1;
      %first level: just store
      Store.ST   = ST.ST;
      freqs = ST.freqs; %needed later
    else
      %subsequent level: integrate into existing average
      Store.ST = ((iLevel-1).*Store.ST + ST.ST)./iLevel;
    end
    
    %and loop
    textprogressbar(iLevel./numel(Airs.ret_z).*100)
  end
  textprogressbar('!')
  clear iLevel CFac ST
 
  
  %% find spectral maxima in the two fields
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %find spectral peaks
  PeakStore = NaN([size(Airs.l1_lat),Settings.NPeaks,3]);
  
  textprogressbar('Finding peaks ')
  for iX=1:1:size(Airs.l1_lat,1);
    for iY=1:1:size(Airs.l1_lat,2);
      
      Frame = abs(Store.ST(  :,:,iX,iY));
      
      if Settings.NPeaks > 1; cent = FastPeakFind(Frame,0,Settings.Filt); %find n peaks
      else                    cent = [];                                  %used in logic below
      end 
      
      if numel(cent) == 0;
        %if we either found no peaks or only wanted a single peak,
        %take the overall-maximum as the largest peak
        [~,idx] = max(Frame(:));
        [x,y] = ind2sub(size(Frame),idx);
        cent(2) = x; cent(1) = y;
        clear idx x y
      end
      

      %convert from list to [x1,y1;x2,y2;...];
      Cent(:,1) = cent(1:2:end); %position on k_(along-track)  axis
      Cent(:,2) = cent(2:2:end); %position on k_(across-track) axis
      
      %find values of these indices
      for iZ=1:1:size(Cent,1); Cent(iZ,3) = Frame(Cent(iZ,2),Cent(iZ,1)); end
      
      %remove small values
      Good = find(Cent(:,3) > Settings.Threshold);
      Cent = Cent(Good,:);
      clear Good
      
      %sort the data, select the number of peaks desired, and store
      [~,idx] = sort(Cent(:,3),'desc');
      if size(Cent,1) > Settings.NPeaks; Cent = Cent(idx(1:Settings.NPeaks),:);
      else;                              Cent = Cent(idx,:);
      end
      clear idx
      PeakStore(iX,iY,1:size(Cent,1),:) = Cent;
      clear Cent
      
    end
    textprogressbar(iX./size(Airs.l1_lat,1).*100)
  end
  clear iX iY iZ Frame iMode Working cent
  textprogressbar('!')
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% we now have the N strongest spectral signals at each geographic point,
  % found using a coarse c to provide reasonable spectral fit.
  %now put back and do a finer fit to these frequencies
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %find all the unique scales of importance in the data
  clear Scales
  s1 = unique(PeakStore(:,:,:,1,:)); s1(isnan(s1)) = []; Scales{1} = s1; clear s1;
  s2 = unique(PeakStore(:,:,:,2,:)); s2(isnan(s2)) = []; Scales{2} = s2; clear s2;
  
  %frequencies from original ST, to fix horizontal wavelengths
  f1 = freqs{1};
  f2 = freqs{2};
  
  %create storage arrays for the fine STs
  Store.STs  = complex(NaN([numel(Scales{1}),numel(Scales{2}), ...
                       size(Airs.l1_lat),numel(Airs.ret_z)], ...
                       'single'));
  
  textprogressbar('Computing fine STs ')
  for iLevel=1:1:numel(Airs.ret_z);
    
    %do ST
    ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel), ...
                           'FullST',true,          ...
                           'c',Settings.c2(1:2),   ...
                           'Scales',Scales);
    %store ST
    Store.STs(:,:,:,:,iLevel) = ST.ST;
    
    textprogressbar(iLevel ./ size(Airs.ret_z,1) .* 100);
    
  end
  clear ST iLevel
  textprogressbar('!')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% find the requested cospectra, and store the resulting GW properties
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %results arrays
  Store.A  = NaN([size(Airs.l1_lat),size(Airs.ret_z,1),numel(Settings.Steps),Settings.NPeaks]);
  Store.L  = Store.A;
  Store.F1 = Store.A;
  Store.F2 = Store.A;
  
  
  textprogressbar('Computing cospectra ')
  for iLevel=1:1:size(Airs.ret_z,1);
    for iStep = 1:1:numel(Settings.Steps)
      
      %identify the level-pair to compute a cospectrum over
      ST1 = Store.STs(:,:,:,:,iLevel); %this level is one of them...
      if iLevel+Settings.Steps(iStep) <= size(Airs.ret_z,1);      
        ST2 = Store.STs(:,:,:,:,iLevel+Settings.Steps(iStep));
      else; continue; end
      dZ = Airs.ret_z(iLevel) - Airs.ret_z(iLevel+Settings.Steps(iStep));
      
      %take the cospectrum.
      CoSpectrum = ST1 .* conj(ST2);
      
      %then pull out the indices we need
      for iX=1:1:size(Airs.l1_lat,1);
        for iY=1:1:size(Airs.l1_lat,2);
          for iZ=1:1:Settings.NPeaks;
              if ~isnan(PeakStore(iX,iY,iZ,3));
                
                %find new indices
                idx1 = closest(PeakStore(iX,iY,iZ,1),Scales{1});
                idx2 = closest(PeakStore(iX,iY,iZ,2),Scales{2});
                
                %cospectral value
                CS = CoSpectrum(idx1,idx2,iX,iY);
                
                %amplitude
                Store.A(iX,iY,iLevel,iStep,iZ) = sqrt(abs(CS));
                
                %vertical wavelength
                dPhi = angle(CS)./(2*pi);
                
                Store.L(iX,iY,iLevel,iStep,iZ) = dZ./dPhi;
                
                %horizontal wavelengths
                Store.F1(iX,iY,iLevel,iStep,iZ) = f1(PeakStore(iX,iY,iZ,2));
                Store.F2(iX,iY,iLevel,iStep,iZ) = f2(PeakStore(iX,iY,iZ,1));
                
            end
          end
        end
      end
    end
    textprogressbar(iLevel ./ size(Airs.ret_z,1) .* 100);
  end;
  clear iLevel iX iY iZ iStep CoSpectrum CS dPhi
  clear idx1 idx2 dZ ST1 ST2 f1 f2 freqs PeakStore
  Store = rmfield(Store,{'STs','ST'});
  textprogressbar('!')


  %combine and return
  NewFields.A  =    squeeze(Store.A( :,:,:,1,1));
  NewFields.m  = 1./squeeze(Store.L( :,:,:,1,1));
  NewFields.F1 =    squeeze(Store.F1(:,:,:,1,1));
  NewFields.F2 =    squeeze(Store.F2(:,:,:,1,1));

return
