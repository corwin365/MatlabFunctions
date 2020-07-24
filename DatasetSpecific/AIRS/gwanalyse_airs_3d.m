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
%    TwoDPlusOne     (logical,          false)  *ADDITIONALLY* compute vertical wavelengths with 2D+1 method.
%    TwoDPlusOne_ind (logical,          false)  *ADDITIONALLY* compute vertical wavelengths with 2D+1 method, using independent fits rather than the 3DST fits
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
if Status ==0; CheckZRange = @(x) validateattributes(x,{'numeric'},{'size',[1,2],'>=',   0,'<=', 90});
else           CheckZRange = @(x) validateattributes(x,{'numeric'},{'size',[1,2],'>=',-Inf,'<=',Inf});
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

%TwoDPlusOne is logical
addParameter(p,'TwoDPlusOne',    true,@islogical);  %assumes we *don't* want to compute wavelengths via 2D+1 (this is much slower, but often more accurate)
addParameter(p,'TwoDPlusOne_ind',true,@islogical);  %as above, but computes fits independently rather than using 3DST fits

%MaxWaveLength must be an positive real number
CheckLambda = @(x) validateattributes(x,{'numeric'},{'>=',0});
addParameter(p,'MaxWaveLength',[1,1,1].*99e99,CheckLambda);  %crazy-large default

%MinWaveLength must be an positive real number
addParameter(p,'MinWaveLength',[1,1,1].*0,CheckLambda);  %zero default


%parse the inputs, and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ok, do the parsing!
parse(p,Airs,varargin{:});

%pull out the contents into struct "Inputs", used throughout rest of routine
Input = p.Results;

%tidy up
clearvars -except Input ST Airs Error ErrorInfo
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if requested, also use the 2D+1 method to compute vertical wavelength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if     Input.TwoDPlusOne;    
  
  %standard version - uses fitted horizontal wavelengths from 3DST
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  NewFields = get_2dp1_lambdaz_v2(Airs,ST,Extra,Input);
  %add those fields to the main ST array
  Fields = fieldnames(NewFields);
  for iF=1:1:numel(Fields);
    ST.([Fields{iF},'_2dp1']) = NewFields.(Fields{iF});
  end; clear iF Fields NewFields

end
  
if Input.TwoDPlusOne_ind;
  
  %alternate version - fits horizontal wavelengths independently and uses these
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  NewFields = get_2dp1_lambdaz(Airs,Extra,Input); 
  %add those fields to the main ST array
  Fields = fieldnames(NewFields);
  for iF=1:1:numel(Fields);
    ST.([Fields{iF},'_2dp1_ind']) = NewFields.(Fields{iF});
  end; clear iF Fields NewFields
    
end



%undo height scaling
if Input.HeightScaling;
  Airs.Tp = Airs.Tp ./ CFac;
  ST.A    = ST.A    ./ CFac;
  ST.HA   = ST.HA   ./ CFac;
  ST.R    = ST.R    ./ CFac;
  ST.HR   = ST.HR   ./ CFac;
  
  if Input.TwoDPlusOne;     ST.A_2dp1     = ST.A_2dp1     ./ CFac; end
  if Input.TwoDPlusOne_ind; ST.A_2dp1_ind = ST.A_2dp1_ind ./ CFac; end
  
  clear CFac
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

function NewFields = get_2dp1_lambdaz(Airs,Extra,Input)


  %glue the extra levels we retained to the top and bottom
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

  %take the 2D S-Transform for each level, and retain the information we need later
  Vars = {'F1','F2','A','ST','idx_F1','idx_F2','k','l'};
  for iLevel=1:1:numel(Airs.ret_z)
    
    %do the 2DST
    ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel),'FullST',true,'c',Input.c(1:2), ...
                           'MaxWavelength', Input.MaxWaveLength,'MinWavelength', Input.MinWaveLength);
    
    %find the **indices** of the fitted horizontal waves
    F1s = unique(ST.F1); idx_F1 = ST.F1; for iF=1:1:numel(F1s); idx_F1(ST.F1 == F1s(iF)) = closest(ST.freqs{1},F1s(iF)); end
    F2s = unique(ST.F2); idx_F2 = ST.F2; for iF=1:1:numel(F2s); idx_F2(ST.F2 == F2s(iF)) = closest(ST.freqs{2},F2s(iF)); end
    ST.idx_F1 = idx_F1; ST.idx_F2 = idx_F2; 
    clear F1s F2s iF idx_F1 idx_F2
    
    %store for later. Append in dimension 5 as that's safe for everything ane makes the code cleaner, we can fix it later/
    for iVar=1:1:numel(Vars);
      if iLevel == 1; 
        Store.(Vars{iVar}) = ST.(Vars{iVar}); 
        freqs = ST.freqs; %we need freqs later, but only one copy as it's the same each loop
      else
        Store.(Vars{iVar}) = cat(5,Store.(Vars{iVar}),ST.(Vars{iVar}));
      end
    end
  end; clear iLevel ST iVar
  
  %trim the dimensionality of all except Store.ST
  for iVar=1:1:numel(Vars)
    if strcmp(Vars{iVar},'ST'); continue; end
    Store.(Vars{iVar}) = permute(Store.(Vars{iVar}),[1,2,5,3,4]);
  end
  clear iLevel iVar

  %produce complex cospectra
  CC = Store.ST.*NaN;
  for iLevel=2:1:numel(Airs.ret_z)-1
    CC(:,:,:,:,iLevel) = Store.ST(:,:,:,:,iLevel-1) .* conj(Store.ST(:,:,:,:,iLevel+1));
  end; clear iLevel

  %pull the peak-amplitude wavelength signals we want out of the array
  %there's definitely a way to vectorise this, but it's so fast it's not worth the hassle
  sz = size(CC);
  for iX=1:1:sz(3); for iY=1:1:sz(4); for iZ=1:1:sz(5)
    CC(1,1,iX,iY,iZ) = CC(Store.idx_F1(iX,iY,iZ),Store.idx_F2(iX,iY,iZ),iX,iY,iZ);
  end; end; end
  clear iX iY iZ
  CC = squeeze(CC(1,1,:,:,:));
  
  %convert each levels CC values to covarying phase between the two voxels at peak wavelength
  AlldP = angle(CC);  clear CC
  
  %convert phase change to wavelength
  sz = size(Airs.ret_z); if sz(2) > sz(1); Airs.ret_z = Airs.ret_z'; end
  dZ = [diff(Airs.ret_z)+ circshift(diff(Airs.ret_z),1);0]; %the levels that this makes wonky will be removed later
  sz = size(AlldP);
  Lambda = permute(repmat(dZ,1,sz(1),sz(2)),[2,3,1])./AlldP.*2*pi;
  clear sz AlldP dZ
  
  %force sign convention 
  Negative = find(Lambda < 0);
  Lambda(  Negative) = -1.*Lambda(  Negative);
  Store.F1(Negative) = -1.*Store.F1(Negative);
  Store.F2(Negative) = -1.*Store.F2(Negative);
  Store.k( Negative) = -1.*Store.k( Negative);
  Store.l( Negative) = -1.*Store.l( Negative);

  
  %drop the extra levels we used for the phase fitting, and prepare to return
  Vars = {'F1','F2','F3','k','l','A','m'};
  for iVar=1:1:numel(Vars)
    if strcmp(Vars{iVar},'F3') | strcmp(Vars{iVar},'m'); V = 1./Lambda;
    else V = Store.(Vars{iVar}); 
    end
    
    if isfield(Extra,'Bottom'); V = V(:,:,numel(Extra.Bottom.ret_z)+1:end); end
    if isfield(Extra,'Top');    V = V(:,:,1:end-numel(Extra.Top.ret_z)); end
    
    NewFields.(Vars{iVar}) = V;
  end; clear iVar V

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D+1 st, version 2. This uses the 3DST horizontal wavenumber fits.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NewFields = get_2dp1_lambdaz_v2(Airs,ST3D,Extra,Input)

  %glue the extra levels we retained to the top and bottom
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
  
  
  %take the 2D S-Transform for each level, compute complex cospectra, and
  %retain the information we need later the logic here is a bit convoluted 
  %this is to save memory, as the high-order arrays we're using
  %absolutely consume the stuff.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %identify the 2D scales that were identified as important at any height 
  %in the 3DST analysis, and analyse only these in the loop below
  Scales = {unique(ST3D.scales{2}),unique(ST3D.scales{1})};  
  
  %create storage array for 2DST output
  %the '3' is because, to reduce memory use, we only retain the 
  %3 levels we need at any point to compute the needed cospectra
  Store.ST = complex(NaN(numel(Scales{1}),numel(Scales{2}), ...
                     size(Airs.Tp,1),size(Airs.Tp,2),   ...
                     3,                                 ...
                     'single'));
  Store.ST = complex(Store.ST);
  
  %also create array for the reduced complex cospectra that we will retain
  CC = complex(NaN(size(Airs.Tp,1),size(Airs.Tp,2),   ...
               numel(Airs.ret_z),'single'));
  
 
  for iLevel=1:1:numel(Airs.ret_z)
    
    %do the 2DST
    ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel),'FullST',true,'c',Input.c(1:2),'Scales',Scales);
              
    %store the data we need for cospectral computation
    if iLevel <=3;
      %just store the data
      Store.ST(:,:,:,:,iLevel) = ST.ST;
    else
      %shunt everything down and then replace the top level
      Store.ST = circshift(Store.ST,-1,5);
      Store.ST(:,:,:,:,3) = ST.ST;
    end
    freqs = ST.freqs;    
    
    
    %if this is the first loop, use the frequencies to match the indices in the 3DST output 
    %to the indices in the 2DST output. Should be identical in value, just ordered differently
    if iLevel == 1;
      
      %find the matches
      F1_3D = ST3D.F1;  idx_F1_3D = NaN.*F1_3D;
      F2_3D = ST3D.F2;  idx_F2_3D = NaN.*F2_3D;
      for iPoint = 1:1:numel(F1_3D);
        idx_F1_3D(iPoint) = closest(ST.freqs{1},F1_3D(iPoint));
        idx_F2_3D(iPoint) = closest(ST.freqs{2},F2_3D(iPoint));
      end; clear iPoint F1_3D F2_3D
      
      %add on levels corresponding to the Extra levels we supplied for fitting
      sz = size(idx_F1_3D);
      if isfield(Extra,'Bottom');
        idx_F1_3D = cat(3,NaN(sz(1),sz(2),numel(Extra.Bottom.ret_z)),idx_F1_3D);
        idx_F2_3D = cat(3,NaN(sz(1),sz(2),numel(Extra.Bottom.ret_z)),idx_F2_3D);
      end
      if isfield(Extra,'Top');
        idx_F1_3D = cat(3,idx_F1_3D,NaN(sz(1),sz(2),numel(Extra.Top.ret_z)));
        idx_F2_3D = cat(3,idx_F2_3D,NaN(sz(1),sz(2),numel(Extra.Top.ret_z)));
      end 
      clear sz
    end
  
    %only proceed to compute cospectra if we have enough levels
    if iLevel < 3; continue; end
    
    %produce a complex cospectrum of the level below and level above the level of interest
    LevCC = Store.ST(:,:,:,:,1) .* conj(Store.ST(:,:,:,:,3)); %this is for level ******iLevel-1*****

    %extract just the points that correspond to the dominant waves at each point
    sz = size(LevCC);
    a = Store.ST(:,:,:,:,1);
    b = Store.ST(:,:,:,:,3);
    for iX=1:1:sz(3);
      for iY=1:1:sz(4)
        a(1,1,iX,iY) = a(idx_F1_3D(iX,iY,iLevel-1),...
                         idx_F2_3D(iX,iY,iLevel-1),...
                         iX,iY);
        b(1,1,iX,iY) = b(idx_F1_3D(iX,iY,iLevel-1),...
                         idx_F2_3D(iX,iY,iLevel-1),...
                         iX,iY);
        
% % %         LevCC(1,1,iX,iY) = LevCC(idx_F1_3D(iX,iY,iLevel-1),...
% % %                                  idx_F2_3D(iX,iY,iLevel-1),...
% % %                                  iX,iY);
      end
    end
    
    
    %and store
    CC(:,:,iLevel-1) = a(1,1,:,:) .* conj(b(1,1,:,:));%LevCC(1,1,:,:);
    clear sz Lev
     
  end; clear iLevel ST iVar idx_F1_3D idx_F2_3D iLevel 

  %drop the extra levels we used for the phase fitting
  sz = size(Airs.ret_z); if sz(2) > sz(1); Airs.ret_z = Airs.ret_z'; end
  dZ = [diff(Airs.ret_z)+ circshift(diff(Airs.ret_z),1);0];
  if isfield(Extra,'Bottom');
    if numel(Extra.Bottom.ret_z) > 0;
      CC = CC(:,:,numel(Extra.Bottom.ret_z)+1:end);
      dZ = dZ(numel(Extra.Bottom.ret_z)+1:end);
    end
  end
  if isfield(Extra,'Top');
    if numel(Extra.Top.ret_z) > 0;
      CC = CC(:,:,1:end-numel(Extra.Top.ret_z));
      dZ = dZ(1:end-numel(Extra.Top.ret_z));
    end
  end

  %retain the phase difference of these fits
  sz = size(CC);
  CC        = reshape(CC,sz(1)*sz(2),sz(3));
  PhiStore  = NaN(sz(1)*sz(2),sz(3)); AmpStore = PhiStore;
  for iLevel=1:1:sz(3)
    for iPoint=1:1:sz(1)*sz(2);
      Phi = CC(iPoint,iLevel);
      PhiStore(iPoint,iLevel) = Phi;
      AmpStore(iPoint,iLevel) = sqrt(abs(Phi));
    end
  end
  PhiStore = reshape(PhiStore,sz(1),sz(2),sz(3));
  AmpStore = reshape(AmpStore,sz(1),sz(2),sz(3));

  clear sz CC iPoint Phi CC


  %convert each levels CC values to covarying phase between the two voxels at peak wavelength
  AlldP = angle(PhiStore);  clear PhiStore
  
  %convert phase change to wavelength
  sz = size(AlldP);
  Lambda = abs(permute(repmat(dZ,1,sz(1),sz(2)),[2,3,1])./AlldP.*2*pi);
  clear sz AlldP dZ
  


  %and return
  NewFields.F3 = 1./Lambda;
  NewFields.m  = 1./Lambda;
  NewFields.A  = AmpStore;

return



