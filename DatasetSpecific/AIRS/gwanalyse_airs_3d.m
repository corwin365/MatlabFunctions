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
%note that the individual flags WILL NOT BE SANITY-CHECKED - cave usor...

%first, check if we're doing a 2D call. Then parse out the fields
a = find(strcmp(varargin,'TwoDPlusOne'));
if numel(a) > 0;
  if varargin{a+1} == 1;
    
    %defaults - override if they don't exist in the call
    TwoDSettings.c1          = [1,1].*1;               %c for  first pass (coarse in space, fine in freq)
    TwoDSettings.c2          = [1,1].*0.25;              %c for second pass (coarse in freq, fine in space)
    TwoDSettings.NPeaks      = 1;                        %number of spectral peaks to identify for each spatial point
    TwoDSettings.Threshold   = 0;                        %threshold in spectral space to be identified as a peak. This is an amplitude for a lone 2DST (as opposed to a cospectrum)
    TwoDSettings.Filt        = fspecial('gaussian',5,1); %characteristic size of point in spectral space. This is a little fatter than the default, to avoid very close peaks.
    TwoDSettings.Thin        = 1;                        %thin out the number of scales (large runtime reduction, but changes the results)
    TwoDSettings.Steps       = 2;%[1,2,3,4,5];           %number of steps to take phase difference over. '0' takes it from a basis level, defined above, while nonzero values use the phase shift with that many levels *above*
    TwoDSettings.Weight      = 0;                        %height-weight the vertical layers

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
if Input.TwoDPlusOne; 
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
    
  NewFields = get_airs_2dp1_lambdaz(Airs,Extra,TwoDSettings);
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

