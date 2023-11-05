function Data =  get_limbsounders(TimeRange,Instrument,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%general-purpose function to load and format limb sounder data
%loads a specific set of instruments in the formats I store them, so
%may not work on your system!
%
%Note that the guts of the programme is handled by an external module
%file for each instrument that loads and formats the specific data.
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/08/15
%
%changes:
%  2023/09/13 added ACE-FTS as a valid instrument
%  2023/09/13 added ability to select profiles by lat/lon
%  2023/09/19 added MIPAS and SOFIE
%  2023/10/30 added filenames and profile numbers for backtracking to raw data 
%  2023/11/05 split off specific instrument cases into module files
%
%inputs:
%  required:
%    TimeRange [double] - time range, in Matlab units. Details of how this is handled can be specified using 'TimeHandling' flag below.
%    Instrument [string] - instrument name to load, from specified list (see InstInfo struct below)
%
%  optional:
%
%    VarName         (type,           default)  description
%    -----------------------------------------------------------------------------
%    AdditionalVars  (cell,                {})  list of additional variables to extract, if available  
%    OriginalZ       (logical,          false)  return data on original vertical grid rather than interpolated to common scale.
%    KeepOutliers    (logical,          false)  don't remove outliers from the data. NOTE THAT BY DEFAULT THEY WILL BE REMOVED.
%    HeightScale     (numeric,      18:0.5:60)  heightscale to interpolate the data onto, in km, if OriginalZ is not set
%    HeightRange     (numeric,      [0,99e99])  height range to clip data to. Combines with HeightScale, but is more useful with OriginalZ.
%    LatRange        (numeric,       [-90,90])  latitude  range to select. Maximally permissive - allows profiles which enter the box at any height.
%    LonRange        (numeric,     [-180,180])  longitude range to select. Also maximally permissive.
%    FileSource      (logical,          false)  pass out original point locations as file list plus for each point a file and profile number
%    TimeHandling    (numeric,              3)  see list below
%
%TimeHandling options:
% 1. absolutely strictly - (e.g.) datenum(2010,1,[1,2])     will include all of 2010/01/01 and the first second of 2010/01/02
% 2. generously          - (e.g.) datenum(2010,1,[1.5,2.1]) will include all of 2010/01/01 and 2010/01/02, but not the first second of 2010/01/03 
% 3. fuzzy-ended         - (e.g.) datenum(2010,1,[1,2])     will include all of 2010/01/01 and 2010/01/02, but not the first second of 2010/01/03 
%(3) is the default because it behaves almost as intuitively as (1) but reduces rutime by not loading a whole day of data to grab one second
%
%
%outputs:
%   Data: struct containing all variables, on a [profiles x height] grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% instruments we can use this on, and their properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ACE-FTS
%only currently loads temperature and pressure, as generated from GLC data
InstInfo.ACE.TimeRange      = [datenum(2004,1,35),datenum(9999,999,999)]; %still running at time of writing
InstInfo.ACE.HeightRange    = [18,125]; %data are available outside this range but are entirely a priori
InstInfo.ACE.Path           = [LocalDataDir,'/ACE/raw/'];

%GNSS
InstInfo.GNSS.TimeRange     = [datenum(2002,1,1),datenum(9999,999,999)]; %still running at time of writing
InstInfo.GNSS.HeightRange   = [0,40];
InstInfo.GNSS.Path          = [LocalDataDir,'/GNSS/raw/'];

%HIRDLS
InstInfo.HIRDLS.TimeRange   = [datenum(2005,1,29),datenum(2008,1,77)];
InstInfo.HIRDLS.HeightRange = [0,80];
InstInfo.HIRDLS.Path        = [LocalDataDir,'/HIRDLS/raw/'];

%MIPAS
InstInfo.MIPAS.TimeRange   = [datenum(2002,7,2),datenum(2012,4,8)];
InstInfo.MIPAS.HeightRange = [0,80];
InstInfo.MIPAS.Path        = [LocalDataDir,'/MIPAS/raw/'];

%MLS
%only currently loads temperature files, due to the way I have them stored
InstInfo.MLS.TimeRange      = [datenum(2004,1,275),datenum(9999,999,999)]; %still running at time of writing
InstInfo.MLS.HeightRange    = [0,100];
InstInfo.MLS.Path           = [LocalDataDir,'/MLS/raw/'];

%SABER
InstInfo.SABER.TimeRange   = [datenum(2002,1,1),datenum(9999,999,999)]; %still running at time of writing
InstInfo.SABER.HeightRange = [0,120];
InstInfo.SABER.Path        = [LocalDataDir,'/SABER/raw/'];

%SOFIE
InstInfo.SOFIE.TimeRange   = [datenum(2007,1,135),datenum(9999,999,999)]; %still running at time of writing
InstInfo.SOFIE.HeightRange = [10,110]; %this is the range of the mission time/height cross-sections on their website
InstInfo.SOFIE.Path        = [LocalDataDir,'/SOFIE/raw/'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%create parser object
%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

%validation of inputs
%%%%%%%%%%%%%%%%%%%%%%

%date
CheckDates  = @(x) validateattributes(x,{'numeric'},{'>=',datenum(1979,1,1)}); %must be in the satellite era
addRequired(p,'TimeRange',CheckDates); %time range

%instrument(s)
CheckInst  = @(x) validateStringParameter(x,fieldnames(InstInfo),mfilename,Instrument);
addRequired(p,'Instrument',CheckInst)

%height range and height scale
CheckHeights = @(x) validateattributes(x,{'numeric'},{'>=',0}); 
addParameter(p,'HeightScale',18:0.5:60,CheckHeights)
addParameter(p,'HeightRange',[0,99e99],CheckHeights)

%latrange and lonrange
addParameter(p,'LatRange',[ -90, 90],@(x) validateattributes(x,{'numeric'},{'>=', -90,'<=', 90}))
addParameter(p,'LonRange',[-180,180],@(x) validateattributes(x,{'numeric'},{'>=',-180,'<=',180}))

%additional variables
addParameter(p,'AdditionalVars',   {},@iscell)
addParameter(p,     'OriginalZ',false,@islogical)
addParameter(p,  'KeepOutliers',false,@islogical)
addParameter(p,    'FileSource',false,@islogical)
addParameter(p,    'StrictTime',false,@islogical)
addParameter(p,  'TimeHandling',    3,@isnumeric)


%parse inputs and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,TimeRange,Instrument,varargin{:})

%pull out the contents into struct "Settings", used throughout rest of routine
Settings = p.Results;
Settings.TimeRange = TimeRange;
clearvars -except InstInfo Settings

%extract just the metadata for the instrument we want
InstInfo  = InstInfo.(Settings.Instrument);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time range handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%additional check on dates - must be either a single value (a day to load all of) or two values (start and end time)
if numel(Settings.TimeRange)  > 2;  
  error('TimeRange must be either one value (a day of interest) or two values (a time range)');
end

%additional check on dates - must be in valid range for instrument
if   min(Settings.TimeRange) > InstInfo.TimeRange(2) ...
   | max(Settings.TimeRange) < InstInfo.TimeRange(1)
  error(['Data for ',Instrument,' is only available from ',datestr(InstInfo.TimeRange(1)),' to ',datestr(InstInfo.TimeRange(2))]);
end

%duplicate times if single value given
if numel(Settings.TimeRange) == 1; Settings.TimeRange = [1,1].*Settings.TimeRange; end

%handle time strictness
switch Settings.TimeHandling
  case 1; %do nothing, output should be literally what we asked for
  case 2;
    %expand out to include a most generous range of days the user could have meant
    Settings.TimeRange(1) = floor(Settings.TimeRange(1));
    if mod(Settings.TimeRange(2),1) == 0; Settings.TimeRange(2) = Settings.TimeRange(2)+1-1e-8;
    else                                  Settings.TimeRange(2) = ceil(Settings.TimeRange(2))-1e-8;
    end
  case 3;
    %if the last entry is an integer, feather it slightly to avoid loading an extra day
    if mod(Settings.TimeRange(2),1) == 0; Settings.TimeRange(2) = Settings.TimeRange(2)+1-1e-8; end
  otherwise
    error(['Invalid time handling option chosen'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading - instrument specific, see modules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%at the end of this we want a struct called Data containing the following
%variables on a grid of [profiles x levels]
%
%Lat  - latitude
%Lon  - longitude, using -180 to 180 range
%Time - time, in Matlab units
%Temp - temperature, in K
%Pres - pressure, in hPa
%Alt  - height, in km
% + any additional vars requested if valid for that instrument

%list of variables
Vars = [{'Lat','Lon','Time','Temp','Pres','Alt','SourceProf','SourceFile'},Settings.AdditionalVars];

%get the data for this instrument using the appropriate module
switch Settings.Instrument
  case {'ACE','GNSS'}; [Data,FileList] = module_load_ACE_GNSS(Settings,InstInfo,Vars);
  case 'HIRDLS';       [Data,FileList] = module_load_HIRDLS(  Settings,InstInfo,Vars);
  case 'MIPAS';        [Data,FileList] = module_load_MIPAS(   Settings,InstInfo,Vars);
  case 'MLS';          [Data,FileList] = module_load_MLS(     Settings,InstInfo,Vars);
  case 'SABER';        [Data,FileList] = module_load_SABER(   Settings,InstInfo,Vars);
  case 'SOFIE';        [Data,FileList] = module_load_SOFIE(   Settings,InstInfo,Vars);
  otherwise
    disp(['Instrument ',Settings.Instrument,' not currently handled by this function, terminating'])
    return
end
clear FileCount;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the data to chosen height scale
%  and make lons -180 to 180
%  select actual time range specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time range: do first to reduce size for later steps
%do at profile level to avoid breaking profiles
Data = reduce_struct(Data,inrange(nanmean(Data.Time,2),Settings.TimeRange),[],1);


if Settings.OriginalZ == false

  %height scale
  sz = [size(Data.Lat,1),numel(Settings.HeightScale)];
  Data2 = struct();
  for iVar=1:1:numel(Vars);
    a = NaN(sz);
    b = Data.(Vars{iVar});
    for iProf=1:1:sz(1)

      %some extra handling here to deal with bad data, but all we're actualy doing is linear interpolation
      Good = find(isfinite(Data.Alt(iProf,:)) ~=0);
      [~,uidx] = unique(Data.Alt(iProf,:));
      Good = intersect(Good,uidx);
      if numel(Good) > 2;  a(iProf,:) = interp1(Data.Alt(iProf,Good),b(iProf,Good),Settings.HeightScale); end

    end;
    Data2.(Vars{iVar}) = a;
  end
  Data = Data2;
  clear iProf iVar Data2 a b sz Good uidx
end

%longitude
Data.Lon(Data.Lon > 180) = Data.Lon(Data.Lon > 180)-360;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postprocessing based on input options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove outliers?
%this section does the time filtering as well - if it isn't run you'll get whole days only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.KeepOutliers == 0;

  %create one big array of badness, then remove from ALL variables in one pass
  Bad = [];

  %latitude and longitude have physical limits
  Bad = [Bad;find(Data.Lat <  -90 | Data.Lat >  90 | Data.Lon < -180 | Data.Lon > 180)];
  
  %altitude should always be in the specified range
  if Settings.OriginalZ == false;  Bad = [Bad;find(Data.Alt  < min(Settings.HeightScale) | Data.Alt  > max(Settings.HeightScale))]; end

  %time should always be in the specified range
  Bad = [Bad;find(Data.Time < min(Settings.TimeRange  ) | Data.Time > max(Settings.TimeRange  ))];  

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

  clear Bad Fields F iF
end


%latitude and longitude filtering?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(Settings.LatRange) ~=  -90 || max(Settings.LatRange) ~=  90 ...
 | min(Settings.LonRange) ~= -180 || max(Settings.LonRange) ~= 180

  %find the min and max of lon and lat for each profile
  Limits = [min(Data.Lon,[],2),max(Data.Lon,[],2), ...
            min(Data.Lat,[],2),max(Data.Lat,[],2)];

  %generate a list of profiles where any point falls inside the limits
  Bad = [];
  Bad = [Bad;find(Limits(:,1) > max(Settings.LonRange))];
  Bad = [Bad;find(Limits(:,2) < min(Settings.LonRange))];
  Bad = [Bad;find(Limits(:,3) > max(Settings.LatRange))];
  Bad = [Bad;find(Limits(:,4) < min(Settings.LatRange))];

  Good = 1:1:size(Limits,1); Good(Bad) = [];
  Data = reduce_struct(Data,Good,[],1);

end


%do we want to know where the data came from?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.FileSource == 1;
  %add the filelist to the structure
  Data.OriginalFiles = FileList;
else
  %remove tracker variables
  Data = rmfield(Data,{'SourceFile','SourceProf'});
end



clear FileList;


%end of programme
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to validate list of allowed instruments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validateStringParameter(varargin)
validatestring(varargin{:});
end
