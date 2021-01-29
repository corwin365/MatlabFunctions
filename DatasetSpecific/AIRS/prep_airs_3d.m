function [Airs,Spacing,Error,ErrorInfo,MinorErrorInfo] = prep_airs_3d(DateNum,GranuleId,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%generalised function to load AIRS granules produced by the Hoffmann and
%Alexander (2009) retrieval and prepare them for spectral analysis
%
%Corwin Wright, c.wright@bath.ac.uk, 06/AUG/2019
% 
%includes the following functions written by others:
% Neil Hindley:    getnet, nph_haversine (n.hindley@bath.ac.uk)
% Anthony Kendall: date2doy (https://uk.mathworks.com /matlabcentral/fileexchange/18316-date-to-decimal-day-of-year)
% Anil Gannepali:  smoothn, ndgaussian, gridnd, padreplicate (https://uk.mathworks.com/matlabcentral/fileexchange/725-smoothn)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%The function assumes a file tree of the form:
%
%  FullDataDir/YYYY/DDD/airs_YYYY_DDD_GGG.nc
%
%  where:
%     FullDataDir is an input variable, with a default set below
%     YYYY is the year (4-digit) 
%     DDD the day-of-year (3-digit)
%     GGG the granule number (3-digit)
%
%UNLESS a data structure in the right format is supplied (see InputStruct in 'special' inputs, below)
%
%%%%%%%%%%%
%outputs:
%%%%%%%%%%%
%
%  Airs:  AIRS data, provided no critical errors hit. Empty struct otherwise.
%  Spacing: XT, AT and z spacing in km. See Note 1, below.
%  Error: error state; 0: no error; 1: critical error; 2: minor error
%  ErrorInfo: critical error details 
%  MinorErrorInfo: minor error details (cell array of all encountered)
%
%Note 1:
%    for XT and AT, this will always be the mean, and should be reasonably
%    accurate. For z, it will be the *modal* value if only 2D interpolation
%    is used, and the mean if 3D interpolation (not yet implemented) is used.
%
%%%%%%%%%%%
%inputs:
%%%%%%%%%%%
%
 %required (case-sensitive):
 %%%%%%%%%%%
%    DateNum:  (numeric)  Matlab-format date of the granule desired
%    Granule:  (numeric)  number of the granule in that day
%    See Note 2 on FullDataDir, below - this is REQUIRED under most system configs, but fed in as an optional variable
%
 %special:
 %%%%%%%%%%%
%    RelDataDir: (string)    path relative to 'LocalDataDir' of the AIRS data storage tree (Note 2)
%    FullDataDir:  (string)  (OVERRIDES RELDATADIR) path relative to root of the AIRS data storage tree (Note 2)
%    InputStruct:   (string) (OVERRIDES RELDATADIR AND FULLDATADIR) an input struct of Airs data (Note 3)
%
% Note 2: 
%     FullDataDir is ***required*** if you don't have a LocalDataDir file of
%     the type I use to manage running my code across various systems. If you
%     don't know what this comment means, you probably don't have one, in which
%     case you MUST set this variable unless you use an InputStruct
%
% Note 3:
%     if InputStruct is used, Datenum and Granule are ignored but required
%     set them to dummy numeric values just to pass the initial check
%
%
 %optional (all case-insensitive)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    VarName       (type,  default)  description
%    -----------------------------------------------------------------------------
%    NXT           (numeric,    90)  number of cross-track points to interpolate to.
%    NAT           (numeric,   135)  number of along-track points to interpolate to.
%    dXT           (numeric,   NaN)  (OVERRIDES NXT) distance between each XT point to interpolate to 
%    dAT           (numeric,   NaN)  (OVERRIDES NAT) distance between each AT point to interpolate to 
%    NoDetrend     (logical, false)  do not detrend the data.
%    DayNightFlag  (logical, false)  compute day/night flags.
%    KeepOldTime   (logical, false)  keep copy of AIRS time in original format.
%    Interpolant   (string, linear)  type of 2D interpolant to use for regridding AIRS.
%    Extrapolant   (string,   none)  type of 2D extrapolant to use for regridding AIRS.
%    PreSmooth     (array, [1,1,1])  size of boxcar smoother to pre-apply to the data (after interp, before detrend).
%    DetrendMethod (numeric,     2)  detrending method to use (2 is faster)
%    LoadOnly      (logical, false)  load and return data only
%    NoISCheck     (logical, false)  don't check the validity of an input structure (useful eg for merged multiple granules where the size is wrong)
%    Python        (logical, false)  Python can't handle recursive structs - if called from there, then remove the MetaData field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%EXAMPLE ONE:
%
% Airs = prep_airs_3d(datenum(2010,6,8),53)
%
%will:
%  - load the AIRS data for granule 53 on the 8th of June 2010
%  - add fields Tp (temperature perturbations) and BG (background T)
%  - use 90 cross-track  and 135 along-track evenly-spaced elements 
%  - not produce error information 
%
%EXAMPLE TWO:
%
% [Airs,Error,EI,MEI] = prep_airs_3d(datenum(2002,1,242),225,'dXT',10,'DayNightFlag',true,'Interpolant','nearest')
%
%will:
%  - load the AIRS data for granule 225 on day 242 of 2002
%  - use an even 10km across-track spacing
%  - use 135 evenly-spaced along-track elements (i.e. default)
%  - use a nearest-neighbour interpolant for regridding the AIRS data
%  - add variable DayNightFlag, equal to 1 for day and 0 for night
%  - add fields Tp (temperature perturbations) and BG (background T)
%  - provide information on all errors 
%
%
%EXAMPLE THREE:
%
% [Airs,Error,EI,MEI]  = prep_airs_3d(0,0,'InputStruct',AirsIn,'NoDetrend',true)
%
%will:
%  - use data from the input structure AirsIn, ignoring the DateNum/GranuleId inputs
%  - not generate Tp or BG
%  - interpolate data to 90 cross-track  and 135 along-track evenly-spaced elements 
%  - provide information on all errors
%
%
%%EXAMPLE FOUR:
%
% Airs = prep_airs_3d(datenum(2008,6,10),153,'PreSmooth',[3,3,1],'fulldatadir','/local/data/airs/')
%
%will:
%  - use data from file paths /local/data/airs/YYYY/DDD/airs_YYYY_DDD_GGG.nc 
%  - load the AIRS data for granule 153 on the 10th of June 2008
%  - add fields Tp (temperature perturbations) and BG (background T)
%  - use 90 cross-track  and 135 along-track evenly-spaced elements 
%  - will presmooth the data by 3 voxels in the XT and AT direction (after interpolation, but before detrending)
%  - not produce error information 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set default outputs
Airs = struct();
Error = 0;
ErrorInfo = '';
MinorErrorInfo = {};
Spacing = [NaN,NaN,NaN];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% general input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try

  %create input parser
  %%%%%%%%%%%%%%%%%%%%%%%
  p = inputParser;
  
  %inputs - required
  %%%%%%%%%%%%%%%%%%%
  addRequired(p,'DateNum',  @isnumeric);
  addRequired(p,'GranuleId',@isnumeric);
  
  %inputs - pseudo-optional (i.e. path to files)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %first, check if we have a LocalDataDir function
  if exist(LocalDataDir,'file') == 7;
    %we have a function. We can use either DataDir or FullDataDir as
    %optional flags
    addParameter(p,'RelDataDir',[LocalDataDir,'/AIRS/3d_airs/'],@ischar); %where AIRS data live on my system
    addParameter(p,'FullDataDir','notset',                      @ischar); %set this if you don't have a LocalDataDir function. See Headers, Note 1.
  else
    %we do not have a function. FullDataDir is now a required variable
    %easiest to do this a bit clunkily. 
    if sum(ismember(lower(varargin),'fulldatadir')) ~=1;
      ErrorInfo = 'No LocalDataDir function available and FullDataDir not specified';
      Error = 1;
      return
    else
      %check FORMAT of FullDataDir only
      addParameter(p,'FullDataDir','notset', @ischar)
    end
  end
  
 
  %inputs - fully optional
  %%%%%%%%%%%%%%%%%%%%%%%%%
    

  %NXT must be a positive integer
  CheckNT = @(x) validateattributes(x,{'numeric'},{'nonnegative','integer'});  
  addParameter(p,'NXT', 90,CheckNT);  %assume 90 XT rows unless specified  
  addParameter(p,'NAT',135,CheckNT);  %assume 135 AT rows unless specified  
   
  %dXT must be positive 
  CheckdT = @(x) validateattributes(x,{'numeric'},{'nonnegative'});
  addParameter(p,'dXT',NaN,CheckdT);  %assume we're using NXT instead by default
  addParameter(p,'dAT',NaN,CheckdT);  %assume we're using NAT instead by default

  %NoDetrend,DayNightFlag, LoadOnly and KeepOldTime must be logical
  addParameter(p,'NoDetrend',   'false',@islogical);  %assumes we want to detrend
  addParameter(p,'KeepOldTime', 'false',@islogical);  %assumes we want to discard
  addParameter(p,'DayNightFlag','false',@islogical);  %assumes we want to not compute
  addParameter(p,'LoadOnly',    'false',@islogical);  %assumes we want to preprocess the data  
  addParameter(p,'NoISCheck',   'false',@islogical);  %assumes we want to check any input structures 
  addParameter(p,'Python',      'false',@islogical);  %assumes we aren't calling this from Python
  
  %do we want to interpolate to a constant vertical scale?
  %this will be from the middle of the scale outwards, as AIRS
  %resolution is best there
  CheckZI = @(x) validateattributes(x,{'numeric'},{'nonnegative'}); 
  addParameter(p,'VerticallyInterpolate','true',@islogical);  %assumes we want to 
  addParameter(p,'VerticalSpacing',3,CheckZI);                %assume 3km if we do
  
  
  %Interpolant must be a string from the approved list
  ExpectedInterpolants = {'linear','nearest','natural'};
  addParameter(p,'Interpolant','linear',@(x) any(validatestring(x,ExpectedInterpolants)));  %assumes linear
  
  %Extrapolant must be a string from the approved list
  ExpectedExtrapolants = {'linear','nearest','none'};
  addParameter(p,'Extrapolant','none',@(x) any(validatestring(x,ExpectedExtrapolants)));  %assumes linear  
  
  %InputStruct must be a struct. Validation of contents done later...
  addParameter(p,'InputStruct',struct(),@isstruct);
  
  %PreSmooth must be a three-element vector, with all values odd
  CheckSmooth = @(x) validateattributes(x,{'numeric'},{'size',[1,3],'nonnegative','integer','odd'});
  addParameter(p,'PreSmooth',[1,1,1],CheckSmooth); %assume no smoothing.
  
  %DetrendMethod must be either 1 or 2
  CheckDetrendMethod = @(x) validateattributes(x,{'numeric'},{'>=',1,'<=',2});
  addParameter(p,'DetrendMethod',2,CheckDetrendMethod);  %assumes fast method   

  %do we have the right day?
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %if granuleId is outside the legit range (-1-240), then add or subtract
  %days as appropriate
  if GranuleId < 1 || GranuleId > 240;
    DayShift  = GranuleId./240;
    Sign      = DayShift/abs(DayShift);
    if Sign <= 0; DayShift = -ceil(abs(DayShift));
    else;         DayShift = floor(DayShift);
    end
    DateNum   = DateNum + DayShift;
    GranuleId = mod(GranuleId,240) ;   
  end
  
  %parse the inputs, and tidy up
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  parse(p,DateNum,GranuleId,varargin{:})
  clear CheckNT CheckdT CheckSmooth ExpectedInterpolants CheckDetrendMethod
  

  
  %pull out the contents into struct "Inputs", used throughout rest of routine
  Input = p.Results;
  clearvars -except Input Airs Error ErrorInfo MinorErrorInfo Spacing
  
  
  %special case: dXT overrides NXT 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %in the function airs_regularise(), a positive XT input is a dXT value and a
  %negative input is an NXT value
  %so, to override NXT, we simply copy dXT to NXT and multiply by -1
  if ~isnan(Input.dXT);  Input.NXT = -Input.dXT;end
  Input = rmfield(Input,'dXT');
  
  %special case: dAT overrides NAT (same logic)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ~isnan(Input.dAT);  Input.NAT = -Input.dAT;end
  Input = rmfield(Input,'dAT');
  
  %special case: FullDataDir overrides datadir
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  if strcmp(Input.FullDataDir,'notset') ~= 1
    Input.RelDataDir = Input.FullDataDir;
  end
  Input = rmfield(Input,'FullDataDir');
  
  
catch ERR
  ErrorInfo = getReport( ERR, 'extended', 'hyperlinks', 'on' );
  Error = 1;
  return
end
  





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the file, or check the struct is valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(fieldnames(Input.InputStruct)) == 0;
  
  %no struct supplied - get the data
  [Err,Airs,FilePath] = get_airs_granule(Input.RelDataDir,Input.DateNum,Input.GranuleId);
  
  if Err ~= 0;
    Error = 1;
    ErrorInfo = ['Problem finding granule ',FilePath];
    return
  end
  
  Airs.Source = FilePath;
  
  
else

  %test the validity of the structure, and use it is valid
  if Input.NoISCheck == 0;
    [Airs,Check.Error,Check.ErrorInfo] = check_input_struct(Input);
    
    if Check.Error ~= 0;
      %oh dear. exit.
      Error = Check.Error;
      ErrorInfo = Check.ErrorInfo;
      return
    end
  else
    Airs = Input.InputStruct; %BE CAREFUL!!
  end
  
  %all fine - tidy up and continue
  Input = rmfield(Input,'InputStruct');
  clear Check
  
  Airs.Source = 'InputStruct';
  
end
  



%convert time to Matlab time
if Input.KeepOldTime == 1;
  %make a copy of the original time, if wanted
  Airs.OriginalTime = Airs.l1_time;
end
Airs.l1_time = datenum(2000,1,1,0,0,Airs.l1_time);


if Input.LoadOnly == 1;
  %we only wanted to load the data. end programme now
  return
end
  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regularise the data in the horizontal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  %flip the temp array around into the order XT x AT x z
  Airs.ret_temp = permute(Airs.ret_temp,[2,3,1]);
  
  %regularise and find spacing
  [Airs,Spacing] = regularise_airs(Airs,[Input.NXT,Input.NAT],Input.Interpolant,Input.Extrapolant);

catch
  Error = 1;
  ErrorInfo = 'Problem regularising data';
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the data to a regular vertical grid?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Input.VerticallyInterpolate
  
  %reshape data to do interpolation in one 1d pass
  sz = size(Airs.ret_temp);
  T = reshape(Airs.ret_temp,sz(1)*sz(2),sz(3));
  
  %find middle level. Force this to be a level, not an average of two others
  Middle = median(Airs.ret_z);
  if ~ismember(Middle,Airs.ret_z); 
    [~,zidx] = min(abs(Airs.ret_z -Middle));
    Middle = Airs.ret_z(zidx); clear zidx;
  end
  
  %find middle-relative z-range
  ZRange = [min(Airs.ret_z),max(Airs.ret_z)] - Middle;
  
  %create new scale
  NewZ = (min(ZRange):Input.VerticalSpacing:max(ZRange))' + Middle;
  
  %interpolate
  T = interp1(Airs.ret_z,T',NewZ)';
  
  %reshape back
  T = reshape(T,sz(1),sz(2),numel(NewZ));
  
  %and store
  Airs.ret_temp = T;
  Airs.ret_z    = NewZ;
  Spacing(3)    = Input.VerticalSpacing;
  
  clear T NewZ Middle ZRange sz

else
  
  %vertical spacing is not regular, so set it to NaN
  Spacing(3) = NaN;
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% detrend the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%smooth the data
try
  Airs.ret_temp = smoothn(Airs.ret_temp,Input.PreSmooth);
catch;
  Error =1;
  ErrorInfo = 'Unidentified problem pre-smoothing data';
  return
end



try
  if Input.NoDetrend ~= 1;
    [Airs.Tp,Airs.BG] = airs_4dp_detrend(Airs.ret_temp,1,Input.DetrendMethod);
  end
catch
  Error = 1;
  ErrorInfo = 'Problem detrending data';
  return
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute day/night flag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  if Input.DayNightFlag == 1;
    Airs.DayNightFlag = which_airs_retrieval(Airs.l1_lon,Airs.l1_lat,Airs.l1_time);
  end
catch
  
  %check if lat/lon/time are the same size 
  if  numel(Airs.l1_time) == numel(Airs.l1_lon) ...
   && numel(Airs.l1_time) == numel(Airs.l1_lat) ...
   && numel(Airs.l1_time) == numel(Airs.ret_temp(:,:,1));
     MinorErrorInfo{end+1} = 'Unidentified problem computing day/night flags';
  else
     MinorErrorInfo{end+1} = 'To compute day/night flags, l1_lon, l1_lat and l1_time must all be the same size as each other and as ret_temp.';
  end
  
  Error = 2; 
  return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% done! return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if Python, remove MetaData struct
if Input.Python == 1
  Airs = rmfield(Airs,'MetaData'); 
  Airs = rmfield(Airs,'Source');
end

return


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NON-INBUILT FUNCTIONS USED (except LocalDataDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% airs_4dp_detrend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tp,BG] = airs_4dp_detrend(T,Dim,Method)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%detrends data by fitting a fourth-order polynomial to the data in
  %%the dimensions specified by Dim for each row/col/etc individually
  %
  %uses a stripped-down polynomial fitter (Method 2) for speed
  %
  %Corwin Wright, c.wright@bath.ac.uk, 07/AUG/2019
  %
  %inputs:
  %
  %  T - array to be detrended
  %  Dim - dimension to apply the filter in
  %  Method - 1 for slow standard method, 2 for stripped-down fast method (default 2)
  %
  % in tests, I see literally no difference at all betwene the results of 1
  % and 2, but 1 is retained just in case this is wrong...
  %
  %outputs:
  %  Tp - detrended data
  %  BG - trend which was removed (background)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %reshape the data so only the desired dimension is left
  sz = size(T);
  dims = 1:1:numel(sz);
  order = unique([Dim,dims],'stable');  
  T2 = permute(T,order);
  T2 = reshape(T2,[sz(order(1)),prod(sz(order(2:end)))]);

  %create backround array
  BG = T2.*NaN;    
  
  %%now we need to remove the 4th order polynomial
  %we can do this with either full polyfit or a massively stripped-down version with no error or quality checking
  % in tests, I see no difference whatsoever, so the fast method is the default (method 2)
  %but retain the option of doing it the old way just in case
  if ~exist('Method','var'); Method = 2; end
   
  if Method == 1; 
    %"classical" way, with all error checking using built-in

    for iRow=1:1:size(BG,2)
      p = polyval(polyfit(1:1:size(BG,1),T2(:,iRow)',4),1:1:size(BG,1));
      BG(:,iRow) = p;
    end
    
  elseif Method == 2; 
    %fast way - most stuff stripped out. Seems to work...
    %see comment by Matt Tearle on
    %https://uk.mathworks.com/matlabcentral/answers/1836-multiple-use-of-polyfit-could-i-get-it-faster
    
    x = repmat((1:1:size(BG,1))',1,prod(sz(order(2:end))));
    n = 4;
    m = size(x,2);
    p = zeros(n+1,m);
    for k = 1:m
      M = repmat(x(:,k),1,n+1);
      M = bsxfun(@power,M,0:n);
      p(:,k) = M\T2(:,k);
    end
    p = flip(p,1);   
    for iRow=1:1:size(BG,2); BG(:,iRow) = polyval(squeeze(p(:,iRow)),1:1:size(BG,1)); end

  end
  
  %back to common code
  
  %put back into original shape
  BG = reshape(BG,[sz(order(1)),(sz(order(2:end)))]);
  BG = permute(BG,order);
 
  %compute perturbations and return
  Tp = T - BG;

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regularise_airs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% put AIRS on a regular distance grid, same size as input.
% expect t to be XTxAT or XTxATxZ


function  [Airs,Spacing] = regularise_airs(Airs,Spacing,Interpolant,Extrapolant);


%rename vars to save typing
Lon  = Airs.l1_lon;
Lat  = Airs.l1_lat;
Time = Airs.l1_time;
T    = Airs.ret_temp;

%find granule size. use mean of all rows/cols, as there are slight variations between them
AlongTrackSize  = min(nph_haversine([Airs.l1_lat(:,  1),Airs.l1_lon(  :,1)], ...
                                     [Airs.l1_lat(:,end),Airs.l1_lon(:,end)]));
AcrossTrackSize = min(nph_haversine([Airs.l1_lat(1,  :)',Airs.l1_lon(  1,:)'], ...
                                     [Airs.l1_lat(end,:)',Airs.l1_lon(end,:)']));


%define a new grid relative to granule bottom left
dXT = Spacing(1);
dAT = Spacing(2);



if dXT > 0; 
  NewXT = linspace(0,AcrossTrackSize,dXT);
else
  NewXT = 0:abs(dXT):AcrossTrackSize;
end

if dAT > 0; 
  NewAT = linspace(0,AlongTrackSize,dAT);
else
  NewAT = 0:abs(dAT):AlongTrackSize;
end

[lonout,latout,ti, ...
 a,b,timeout]           = put_airs_on_regular_grid(Lon,Lat,Time,T, ...
                                                   mean(diff(NewXT(:))),mean(diff(NewAT(:))), ...
                                                   Interpolant,Extrapolant);


%overwrite the variables, and return
Airs.l1_lon    = lonout;
Airs.l1_lat    = latout;
Airs.l1_time   = timeout;
Airs.ret_temp  = ti;
Spacing        = [a,b];


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lightly-modified form of Neil's function to do this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = put_airs_on_regular_grid(lon,lat,time,t,dXT,dAT,Interpolant,Extrapolant)

xt_mid = ceil(size(t,1)/2);

% first get ALONG TRACK spacing and azimuth:
[~,az_at] = distance(lat(xt_mid,1:end-1),lon(xt_mid,1:end-1),lat(xt_mid,2:end),lon(xt_mid,2:end));
az_at(end+1) = az_at(end);
% d_at(end+1) = d_at(end);
at_spacing = dAT;%mean(deg2km(d_at));
d_at = 0:at_spacing:(at_spacing*(size(t,2)-1));

% now get CROSS TRACK spacing:
[d_xt,~] = distance(repmat(lat(xt_mid,:),size(t,1),1),repmat(lon(xt_mid,:),size(t,1),1),lat,lon);
d_xt = deg2km(mean(d_xt,2))';
d_xt(1:(xt_mid-1)) = -d_xt(1:(xt_mid-1));

% define new grid you want:
xt_vec = linspace(min(d_xt),max(d_xt),size(t,1));
at_vec = d_at;
[XT,~] = ndgrid(xt_vec,at_vec);
xt_spacing = dXT;%mean(diff(xt_vec));

% use reckon to find new lats and lons:
[latout,lonout] = reckon(repmat(lat(xt_mid,:),size(t,1),1),repmat(lon(xt_mid,:),size(t,1),1),km2deg(XT),repmat(az_at+90,size(t,1),1));

% interp each level:
ti = nan(size(t));
for z = 1:size(t,3)
    F = griddedInterpolant({d_xt,d_at},t(:,:,z),Interpolant,Extrapolant);
    ti(:,:,z) = F({xt_vec,at_vec});
end

%and time
F = griddedInterpolant({d_xt,d_at},time,Interpolant,Extrapolant);
timeout = F({xt_vec,at_vec});

% send to outputs:
varargout{1} = lonout;
varargout{2} = latout;
varargout{3} = ti;
varargout{4} = xt_spacing;
varargout{5} = at_spacing;
varargout{6} = timeout;

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get_airs_granule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Error,Data,FilePath] = get_airs_granule(DataDir,DateNum,GranuleId)


  [y,~,~] =datevec(DateNum);
  dn = date2doy(DateNum);

  FilePath = [DataDir,'/',sprintf('%04d',y),'/',sprintf('%03d',dn),'/', ...
              'airs_',sprintf('%04d',y),'_',sprintf('%03d',dn),'_',sprintf('%03d',GranuleId),'.nc'];
  
  if ~exist(FilePath,'file'); Error = 1; Data = struct(); return;
  else
    Data = cjw_readnetCDF(FilePath);
    Error = 0;
  end
    
              

return











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check_input_struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Airs,Error,ErrorInfo] = check_input_struct(Input)

  %default outputs
  Airs = struct();
  Error = 0;
  ErrorInfo = '';


  %check supplied struct has the right fields
  A = isfield(Input.InputStruct,'l1_time');
  B = isfield(Input.InputStruct,'l1_lon'); 
  C = isfield(Input.InputStruct,'l1_lat'); 
  D = isfield(Input.InputStruct,'ret_temp');   
  
  if  A+B+C+D ~= 4;
    Error = 1;
    ErrorInfo = 'Input data structure invalid';
    return
  end
  
  %check the temperature data are formatted correctly
  E = size(Input.InputStruct.ret_temp);
  if E(1) ~= 27 || E(2) ~= 90;
    Error = 1;
    ErrorInfo = 'Error in ret_temp: should be 27 x 90 x N';
    return
  end
  
  %check the lat,lon,time data are formatted correctly
  %writing it out this longhand looks clunky, but according to the Matlab
  %style guide provides a (microscopic) speed boost over a shorter command
  if size(Input.InputStruct.l1_time,1) ~= 90 ...
  || size(Input.InputStruct.l1_lat ,1) ~= 90 ...
  || size(Input.InputStruct.l1_lon ,1) ~= 90;  
    Error = 1;
    ErrorInfo = 'Error in geolocation arrays - should all be 90 x N';
    return
  end
  
  %and the rest aren't used. So we're all fine!
  %well, the data could be bad, eg complex or NaNs, but there reaches a point
  %where it's just overkill to check every possibility...
  clear A B C D E F
  Airs = Input.InputStruct;

  
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
%% getnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function NC = getnet(filepath,varargin)

%% Open netcdf file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nci = ncinfo(filepath);
ncid = netcdf.open(filepath,'nc_nowrite');

% Assign output:
NC = nci;

%% Attributes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NC.Attributes = struct;
NC.Attributes.Global = nci.Attributes;

% ALSO ASSIGN ATTRIBUTES FROM VARIABLES TO THE ATTRIBUTES FIELD:
for i = 1:length(nci.Variables)
    NC.Attributes.(nci.Variables(i).Name) = nci.Variables(i).Attributes;
end

%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NC.Data = struct;

for i = 1:length(nci.Variables)
    if any(strcmpi(varargin,'single'))
        NC.Data.(nci.Variables(i).Name) = single(netcdf.getVar(ncid,i-1));
    else
        NC.Data.(nci.Variables(i).Name) = double(netcdf.getVar(ncid,i-1));
    end
        
end

%% Apply scale factors and offsets if they exist %%%%%%%%%%%%%%%%%%%%%%%%%%
% note these names are for the ECMWF netcdf formats, not clear if they'll
% work for other formats, but you've got the attributes anyway now so you
% can do it afterward.
for i = 1:length(nci.Variables)
    
    % Fill values to NaN:
    if any(strcmpi({nci.Variables(i).Attributes.Name},'_FillValue'))
        ind = find(strcmpi({nci.Variables(i).Attributes.Name},'_FillValue'));
        fill_value = nci.Variables(i).Attributes(ind).Value;
        NC.Data.(nci.Variables(i).Name)(NC.Data.(nci.Variables(i).Name) == fill_value) = NaN;
    end
    % missing values to NaN:
    if any(strcmpi({nci.Variables(i).Attributes.Name},'missing_value'))
        ind = find(strcmpi({nci.Variables(i).Attributes.Name},'missing_value'));
        missing_value = nci.Variables(i).Attributes(ind).Value;
        NC.Data.(nci.Variables(i).Name)(NC.Data.(nci.Variables(i).Name) == missing_value) = NaN;
    end
    % if any of the variables have the term 'scale_factor':
    if any(strcmpi({nci.Variables(i).Attributes.Name},'scale_factor'))
        ind = find(strcmpi({nci.Variables(i).Attributes.Name},'scale_factor'));
        NC.Data.(nci.Variables(i).Name) = NC.Data.(nci.Variables(i).Name) .* nci.Variables(i).Attributes(ind).Value;
    end
    % if any of the variables have the term 'add_offset':
    if any(strcmpi({nci.Variables(i).Attributes.Name},'add_offset'))
        ind = find(strcmpi({nci.Variables(i).Attributes.Name},'add_offset'));
        NC.Data.(nci.Variables(i).Name) = NC.Data.(nci.Variables(i).Name) + nci.Variables(i).Attributes(ind).Value;
    end
    
end

%% Close netcdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

netcdf.close(ncid); clear ncid;

return












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cjw_readnetCDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FileContents = cjw_readnetCDF(FileName)
  
  Data = getnet(FileName);
  
  MetaFields = {'Filename','Name','Dimensions','Variables','Attributes','Groups','Format'};
  FileContents = Data.Data;
  
  for iField=1:1:numel(MetaFields)
    FileContents.MetaData.(MetaFields{iField}) = Data.(MetaFields{iField});
  end
  
return












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% which_airs_retrieval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mask = which_airs_retrieval(Lon,Lat,Time)

%identified whether AIRS-3D is using the day or night retrieval, using the same
%logic as Lars' retrieval does itself (shown in C at bottom of function)
%
%Corwin Wright, c.wright@bath.ac.uk, 14/AUG/2018
%
%input is time in Matlab format, latitude, and longitude (array-safe, must all
%be same shape)
%output is a binary mask, with 1 for daytime and 0 for nighttime, of the same
%size as the inputs

%Number of days and fraction with respect to 2000-01-01T12:00Z
D = Time-datenum(2000,1,1,12,0,0);

%Geocentric apparent ecliptic longitude [rad]
g = (357.529 + 0.98560028 .* D) .* pi ./ 180;
q = 280.459 + 0.98564736 .* D;
L = (q + 1.915 .* sin(g) + 0.020 .* sin(2 * g)) .* pi ./ 180;

%Mean obliquity of the ecliptic [rad]
e = (23.439 - 0.00000036 * D) .* pi ./ 180;

% Declination [rad]
dec = asin(sin(e) .* sin(L));

% Right ascension [rad]...
ra = atan2(cos(e) .* sin(L), cos(L));
 
% Greenwich Mean Sidereal Time [h]...
GMST = 18.697374558 + 24.06570982441908 .* D;


%Local Sidereal Time [h]... 
LST = GMST + Lon ./ 15;

% Hour angle [rad]... 
h = LST ./ 12 .* pi - ra;
 
% Convert latitude... 
Lat = Lat .* pi ./ 180;

%hence, solar zenith angle [deg].
SZAd = acos(sin(Lat) .* sin(dec) + cos(Lat) .* cos(dec) .* cos(h)) .* 180 / pi;

%create mask
Mask = SZAd .*NaN;
Mask(SZAd <  96) = 1; %day
Mask(SZAd >= 96) = 0; %night

return












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% date2doy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [doy,fraction] = date2doy(inputDate)
%DATE2DOY  Converts the date to decimal day of the year.
%   [doy,fraction] = date2doy(inputDate)
%
%   Descriptions of Input Variables:
%   inputDate:  Input date as a MATLAB serial datenumber
%
%   Descriptions of Output Variables:
%   doy: Decimal day of the year. For instance, an output of 1.5 would 
%       indicate that the time is noon on January 1.
%   fraction: Outputs the fraction of a year that has elapsed by the input
%       date.
%
%   Example(s):
%   >> [doy,frac] = date2doy(datenum('1/25/2004'));
%
%   See also:

% Author: Anthony Kendall
% Contact: anthony [dot] kendall [at] gmail [dot] com
% Created: 2008-03-11
% Copyright 2008 Michigan State University.

%Want inputs in rowwise format
[doy,fraction] = deal(zeros(size(inputDate)));
inputDate = inputDate(:);

%Parse the inputDate
[dateVector] = datevec(inputDate);

%Set everything in the date vector to 0 except for the year
dateVector(:,2:end) = 0;
dateYearBegin = datenum(dateVector);

%Calculate the day of the year
doyRow = inputDate - dateYearBegin;

%Optionally, calculate the fraction of the year that has elapsed
flagFrac = (nargout==2);
if flagFrac
    %Set the date as the end of the year
    dateVector(:,1) = dateVector(:,1) + 1;
    dateYearEnd = datenum(dateVector);
    fracRow = (doyRow - 1) ./ (dateYearEnd - dateYearBegin);
end

%Fill appropriately-sized output array
doy(:) = doyRow;
if flagFrac
    fraction(:) = fracRow;
end

return











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% smoothn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = smoothn(X,sz,filt,std)

%SMOOTHN Smooth N-D data
%   Y = SMOOTHN(X, SIZE) smooths input data X. The smoothed data is
%       retuirned in Y. SIZE sets the size of the convolution kernel
%       such that LENGTH(SIZE) = NDIMS(X)
%
%   Y = SMOOTHN(X, SIZE, FILTER) Filter can be 'gaussian' or 'box' (default)
%       and determines the convolution kernel.
%
%   Y = SMOOTHN(X, SIZE, FILTER, STD) STD is a vector of standard deviations 
%       one for each dimension, when filter is 'gaussian' (default is 0.65)

%     $Author: ganil $
%     $Date: 2001/09/17 18:54:39 $
%     $Revision: 1.1 $
%     $State: Exp $

if nargin == 2;
  filt = 'b';
elseif nargin == 3;
  std = 0.65;
elseif nargin>4 || nargin<2
  error('Wrong number of input arguments.');
end

% check the correctness of sz
if ndims(sz) > 2 || min(size(sz)) ~= 1
  error('SIZE must be a vector');
elseif length(sz) == 1
  sz = repmat(sz,ndims(X));
elseif ndims(X) ~= length(sz)
  error('SIZE must be a vector of length equal to the dimensionality of X');
end

% check the correctness of std
if filt(1) == 'g'
  if length(std) == 1
    std = std*ones(ndims(X),1);
  elseif ndims(X) ~= length(std)
    error('STD must be a vector of length equal to the dimensionality of X');
  end
  std = std(:)';
end

sz = sz(:)';

% check for appropriate size
padSize = (sz-1)/2;
if ~isequal(padSize, floor(padSize)) || any(padSize<0)
  error('All elements of SIZE must be odd integers >= 1.');
end

% generate the convolution kernel based on the choice of the filter
filt = lower(filt);
if (filt(1) == 'b')
  smooth = ones(sz)/prod(sz); % box filter in N-D
elseif (filt(1) == 'g')
  smooth = ndgaussian(padSize,std); % a gaussian filter in N-D
else
  error('Unknown filter');
end


% pad the data
X = padreplicate(X,padSize);

% perform the convolution
Y = convn(X,smooth,'valid');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = ndgaussian(siz,std)

% Calculate a non-symmetric ND gaussian. Note that STD is scaled to the
% sizes in SIZ as STD = STD.*SIZ


ndim = length(siz);
sizd = cell(ndim,1);

for i = 1:ndim
  sizd{i} = -siz(i):siz(i);
end

grid = gridnd(sizd);
std = reshape(std.*siz,[ones(1,ndim) ndim]);
std(find(siz==0)) = 1; % no smoothing along these dimensions as siz = 0
std = repmat(std,2*siz+1);


h = exp(-sum((grid.*grid)./(2*std.*std),ndim+1));
h = h/sum(h(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function argout = gridnd(argin)

% exactly the same as ndgrid but it accepts only one input argument of 
% type cell and a single output array

nin = length(argin);
nout = nin;

for i=nin:-1:1;
  argin{i} = full(argin{i}); % Make sure everything is full
  siz(i) = numel(argin{i});
end
if length(siz)<nout, siz = [siz ones(1,nout-length(siz))]; end

argout = [];
for i=1:nout;
  x = argin{i}(:); % Extract and reshape as a vector.
  s = siz; s(i) = []; % Remove i-th dimension
  x = reshape(x(:,ones(1,prod(s))),[length(x) s]); % Expand x
  x = permute(x,[2:i 1 i+1:nout]);% Permute to i'th dimension
  argout = cat(nin+1,argout,x);% Concatenate to the output 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b=padreplicate(a, padSize)
%Pad an array by replicating values.
numDims = length(padSize);
idx = cell(numDims,1);
for k = 1:numDims
  M = size(a,k);
  onesVector = ones(1,padSize(k));
  idx{k} = [onesVector 1:M M*onesVector];
end

b = a(idx{:});

return