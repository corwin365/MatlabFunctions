function Data =  get_limbsounders(TimeRange,Instrument,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%general-purpose function to load and format limb sounder data
%loads a specific set of instruments in the formats I store them, so
%may not work on your system!
%
%Note that the guts of the programme is a large number of individual
%instrument cases, so it may be tricky to read
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/08/15
%
%changes:
%  2023/09/13 added ACE-FTS as a valid instrument
%  2023/09/13 added ability to select profiles by lat/lon
%  2023/09/19 added MIPAS and SOFIE
%
%
%
%
%inputs:
%  required:
%    TimeRange [double] - time range, in Matlab units. Usually takes 2-element array
%                         can also take one number, assumed to be a whole day. NOTE 
%                         THAT IF YOU GIVE A RANGE IT TAKES IT LITERALLY, i.e. if you
%                         request from datenum(2020,1,20) to datenum(2020,1,23) it will
%                         end at the midnight at the BEGINNING of 23/1/2020, not the end
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
InstInfo.ACE.Path           = [LocalDataDir,'/ACE/v5.2/'];

%GNSS
InstInfo.GNSS.TimeRange     = [datenum(2002,1,1),datenum(9999,999,999)]; %still running at time of writing
InstInfo.GNSS.HeightRange   = [0,40];
InstInfo.GNSS.Path          = [LocalDataDir,'/GNSS-RO/'];

%HIRDLS
InstInfo.HIRDLS.TimeRange   = [datenum(2005,1,29),datenum(2008,1,77)];
InstInfo.HIRDLS.HeightRange = [0,80];
InstInfo.HIRDLS.Path        = [LocalDataDir,'/HIRDLS'];

%MIPAS
InstInfo.MIPAS.TimeRange   = [datenum(2002,7,2),datenum(2012,4,8)];
InstInfo.MIPAS.HeightRange = [0,80];
InstInfo.MIPAS.Path        = [LocalDataDir,'/MIPAS'];

%MLS
%only currently loads temperature files, due to the way I have them stored
InstInfo.MLS.TimeRange      = [datenum(2004,1,275),datenum(9999,999,999)]; %still running at time of writing
InstInfo.MLS.HeightRange    = [0,100];
InstInfo.MLS.Path           = [LocalDataDir,'/MLS/T/'];

%SABER
InstInfo.SABER.TimeRange   = [datenum(2002,1,1),datenum(9999,999,999)]; %still running at time of writing
InstInfo.SABER.HeightRange = [0,120];
InstInfo.SABER.Path        = [LocalDataDir,'/SABER/rawnc-v2/'];

%SOFIE
InstInfo.SOFIE.TimeRange   = [datenum(2007,1,135),datenum(9999,999,999)]; %still running at time of writing
InstInfo.SOFIE.HeightRange = [10,110]; %this is the range of the mission time/height cross-sections on their website
InstInfo.SOFIE.Path        = [LocalDataDir,'/SOFIE/'];


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


%parse inputs and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,TimeRange,Instrument,varargin{:})

%additional check on dates - must be either a single value (a day to load all of) or two values (start and end time)
if     numel(TimeRange)  > 2;  error('TimeRange must be either one value (a day of interest) or two values (a time range)');
elseif numel(TimeRange) == 1; TimeRange = [floor(TimeRange),floor(TimeRange)+1-1e-8]; end % strictly the latter bound is very slight less than a day...

%additional check on dates - must be in valid range for instrument
if   min(TimeRange) > InstInfo.(Instrument).TimeRange(2) ...
   | max(TimeRange) < InstInfo.(Instrument).TimeRange(1)
  error(['Data for ',Instrument,' is only available from ',datestr(InstInfo.(Instrument).TimeRange(1)),' to ',datestr(InstInfo.(Instrument).TimeRange(2))]);
end


%pull out the contents into struct "Inputs", used throughout rest of routine
Input = p.Results;
Input.TimeRange = TimeRange;
clearvars -except InstInfo Input

%extract just the metadata for the instrument we want
InstInfo  = InstInfo.(Input.Instrument);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading - instrument specific
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
Vars = [{'Lat','Lon','Time','Temp','Pres','Alt'},Input.AdditionalVars];

%empty structure to store data
Data = struct();
for iVar=1:1:numel(Vars); Data.(Vars{iVar}) = []; end


switch Input.Instrument

  case {'ACE','GNSS'}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ACE-FTS v5 and GNSS have very similar file formats, because I
    %made the formats we store them in.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %temporary fix for some anomalous prs values
    if strcmp(Input.Instrument,'GNSS')
      warning('for GNSS, temporarily using alt -> pres instead of pres due to bad stored data')
    end


    for DayNumber=floor(min(Input.TimeRange)):1:floor(max(Input.TimeRange));

      %work out year and day number and hence filepath
      [y,~,~] = datevec(DayNumber); dn = date2doy(DayNumber);
      switch Input.Instrument
        case 'GNSS'; File = [InstInfo.Path,'/',sprintf('%04d',y),'/merged_ro_',sprintf('%04d',y),'_',sprintf('%03d',dn),'.mat'];
        case  'ACE'; File = [InstInfo.Path,'/',sprintf('%04d',y),'/ace_v52_',  sprintf('%04d',y),'d',sprintf('%03d',dn),'.mat'];
      end
      if ~exist(File,'file'); clear y dn File; continue; end

      %load the data
      InstData = load(File);
      clear y dn File

      %get the variables we want
      Store = struct();
      for iVar=1:1:numel(Vars)

        if ~strcmp(Vars{iVar},'Alt')  ...
         & ~strcmp(Vars{iVar},'Time') ...
         & ~strcmp(fieldnames(InstData),Vars{iVar});
          disp(['Variable ',Vars{iVar},' not found, terminating'])
          return
        end

        switch Vars{iVar}
          case 'Time'; 
            switch Input.Instrument
              case 'GNSS'; Store.Time = repmat(InstData.MetaData.time,[1,size(InstData.Lat,2)]);
              case 'ACE';  Store.Time = InstData.Time;
            end
          case 'Alt'; 
            switch Input.Instrument
              case 'GNSS'; Store.Alt  = repmat(InstData.MSL_alt,[size(InstData.Lat,1),1]);
              case 'ACE';  Store.Alt = InstData.Alt;            
            end
          otherwise;   Store.(Vars{iVar}) = InstData.(Vars{iVar});
        end

      end; clear iVar

    %temporary fix for some anomalous prs values
    if strcmp(Input.Instrument,'GNSS')
      Store.Pres = h2p(Store.Alt);
    end

      %store in main repository
      Data = cat_struct(Data,Store,1);

      clear Store iVar

    end; clear DayNumber

  case 'HIRDLS'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIRDLS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for DayNumber=floor(min(Input.TimeRange)):1:floor(max(Input.TimeRange));

      %work out year and day number and hence filepath
      [y,~,~] = datevec(DayNumber); dn = date2doy(DayNumber);
      File = wildcardsearch(InstInfo.Path,['_',sprintf('%04d',y),'d',sprintf('%03d',dn)]);
      if numel(File) == 0; clear y dn File; continue; end

      %load variables we need
      Store = struct();
      for iVar=1:1:numel(Vars)
        switch Vars{iVar}
          case 'Temp'; Store.Temp         = get_HIRDLS(File{1},'Temperature')';
          case 'Lat';  Store.Lat          = get_HIRDLS(File{1},'Latitude');
          case 'Lon';  Store.Lon          = get_HIRDLS(File{1},'Longitude');
          case 'Alt';  Store.Alt          = get_HIRDLS(File{1},'Altitude')'./1000;
          case 'Pres'; Store.Pres         = get_HIRDLS(File{1},'Pressure'); 
          case 'Time'; Store.Time         = datenum(1993,1,1,0,0,get_HIRDLS(File{1},'Time'));
          otherwise;  
            try;   Store.(Vars{iVar}) = get_HIRDLS(File{1},Vars{iVar}); 
            catch; disp(['Variable ',Vars{iVar},' not found, terminating']); return
            end
        end      
      end

      %reshape 1d variables
      Store.Lat  = repmat(Store.Lat,  [1,size(Store.Temp,2)]);
      Store.Lon  = repmat(Store.Lon,  [1,size(Store.Temp,2)]);
      Store.Time = repmat(Store.Time, [1,size(Store.Temp,2)]);
      Store.Pres = repmat(Store.Pres',[size(Store.Temp,1),1]);

      %HIRDLS stores geolocation at the 30km level, but travels while scanning up and down
      %see Wright et al (ACP, 2015) for the logic of what we're going to do here to reverse
      %this choice and get 'true' lat and lon values

      %find the 30km level
      [~,zidx] = min(abs(nanmean(Store.Alt,1)-30));

      %for each profile compute the directional azimuth
      Theta = azimuth(Store.Lat(1:end-1,zidx),Store.Lon(1:end-1,zidx), ...
                      Store.Lat(2:end,  zidx),Store.Lon(2:end,  zidx),'degrees');      
      Theta(end+1) = Theta(end); %close enough for last point

      %hence find the true(ish) location of each point in 2D space
      KmAlongtrackPerKmVertical = 0.6;
      dx = Store.Lat.*NaN;
      for iLev=1:1:size(Store.Lat,2)
        dx(:,iLev) =  KmAlongtrackPerKmVertical.*(iLev-zidx);
        [Store.Lat(:,iLev),Store.Lon(:,iLev)] = reckon(Store.Lat(:,zidx),Store.Lon(:,zidx),km2deg(dx(:,iLev)),Theta);  
      end
      clear zidx Theta KmAlongtrackPerKmVertical dx iLev

      %store in main repository
      Data = cat_struct(Data,Store,1);
      clear Store iVar y dn File


    end; clear DayNumber

  case 'MIPAS'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MIPAS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for DayNumber=floor(min(Input.TimeRange)):1:floor(max(Input.TimeRange));

      %work out filepath
      [y,m,d] = datevec(DayNumber); 
      File = wildcardsearch(InstInfo.Path,['_',sprintf('%04d',y),sprintf('%02d',m),sprintf('%02d',d),'_nom']);
      if numel(File) == 0; clear y dn File; continue; end

      %load data and pull out vars
      Working = rCDF(File{1});
      
      %get variables
      Store = struct();
      for iVar=1:1:numel(Vars)
        switch Vars{iVar}
          case 'Temp'; Store.Temp  = Working.TEM';
          case 'Lat';  Store.Lat   = Working.latitude;
          case 'Lon';  Store.Lon   = Working.longitude;
          case 'Alt';  Store.Alt   = Working.hgt';
          case 'Pres'; Store.Pres  = Working.PRE';
          case 'Time'; Store.Time  = DayNumber+Working.time./24;
          otherwise;  
            try;   Store.(Vars{iVar}) = Working.(Vars{iVar}); 
            catch; disp(['Variable ',Vars{iVar},' not found, terminating']); return
            end
        end      
      end
      clear iVar Working y dn

      %reshape
      NLevs = size(Store.Alt,2); NProfs = numel(Store.Lat);
      Store.Lat  = repmat(Store.Lat, 1,NLevs);
      Store.Lon  = repmat(Store.Lon, 1,NLevs);
      Store.Time = repmat(Store.Time,1,NLevs);
      clear NLevs NProfs


      %store in main repository
      Data = cat_struct(Data,Store,1);

    end; clear DayNimber


  case 'MLS'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MLS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for DayNumber=floor(min(Input.TimeRange)):1:floor(max(Input.TimeRange));

      %work out year and day number and hence filepath
      [y,~,~] = datevec(DayNumber); dn = date2doy(DayNumber);
      File = wildcardsearch(InstInfo.Path,['_',sprintf('%04d',y),'d',sprintf('%03d',dn)]);
      if numel(File) == 0; clear y dn File; continue; end

      %load data for day
      Working = get_MLS(File{1},'Temperature');
      %load variables we need
      Store = struct();
      for iVar=1:1:numel(Vars)
        switch Vars{iVar}
          case 'Temp'; Store.Temp         = Working.L2gpValue';
          case 'Lat';  Store.Lat          = Working.Latitude;
          case 'Lon';  Store.Lon          = Working.Longitude;
          case 'Alt';  Store.Alt          = p2h(Working.Pressure);
          case 'Pres'; Store.Pres         = Working.Pressure;
          case 'Time'; Store.Time         = datenum(1993,1,1,0,0,Working.Time);
          otherwise;  
            try;   Store.(Vars{iVar}) = Working.(Vars{iVar}); 
            catch; disp(['Variable ',Vars{iVar},' not found, terminating']); return
            end
        end      
      end

      %reshape 1d variables
      Store.Lat  = repmat(Store.Lat,  [1,size(Store.Temp,2)]);
      Store.Lon  = repmat(Store.Lon,  [1,size(Store.Temp,2)]);
      Store.Time = repmat(Store.Time, [1,size(Store.Temp,2)]);
      Store.Pres = repmat(Store.Pres',[size(Store.Temp,1),1]);
      Store.Alt  = repmat(Store.Alt', [size(Store.Temp,1),1]);
      for iField=1:1:numel(Input.AdditionalVars);
        F = Store.(Input.AdditionalVars{iField});
        if     numel(F) == size(Store.Temp,2); F = repmat(F',[size(Store.Temp,1),1]);
        elseif numel(F) == size(Store.Temp,1); F = repmat(F, [1,size(Store.Temp,2)]);
        end
        Store.(Input.AdditionalVars{iField}) = F;
      end

      %store in main repository
      Data = cat_struct(Data,Store,1);

      clear Store iVar y dn File Working


    end; clear DayNumber


  case 'SABER';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SABER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %monthly files, ugh
    OldData.Data = struct();
    OldData.Name = '';

    for DayNumber=floor(min(Input.TimeRange))-1:1:floor(max(Input.TimeRange));  %we need to use the day before the start as well as the underlying data format is orbit-based, not day-based

      %identify file and load, but only if we didn't on a previous loop
      [y,~,~] = datevec(DayNumber);
      m = month(datetime(DayNumber,'ConvertFrom','datenum'),'name'); m = m{1};
      File = wildcardsearch(InstInfo.Path,[m,num2str(y)]);

      if numel(File) == 0; clear y m File; continue; end

      if ~strcmp(OldData.Name,File{1});

        %load file
        OldData.Name = File{1};
        OldData.Data = rCDF(File{1});

        %put dates in a useful unit. We'll handle time within day later, as this can be quite computationally expensive.
        for iEl=1:1:numel(OldData.Data.date)
          id = num2str(OldData.Data.date(iEl));
          OldData.Data.Date(iEl) = datenum(str2double(id(1:4)),1,str2double(id(5:7)),0,0,0);
        end

      end

      clear File y m

      %select the day we actually want
      Working = OldData.Data;
      OnThisDay = find(floor(Working.Date) == DayNumber);
      Working = reduce_struct(Working,OnThisDay,{'orbit','date','MetaData'},2);
      Working.orbit = Working.orbit(OnThisDay); Working.date = Working.date(OnThisDay);
      clear OnThisDay

      %now, extract the values we want
      Store.Lat = Working.tplatitude';
      Store.Lon = Working.tplongitude';
      Store.Pres = Working.pressure';
      Store.Alt  = Working.tpaltitude';
      Store.Temp = Working.ktemp';

      %remove bad heights, as they break the interpolation below
      Store.Alt(Store.Alt > 1000) = NaN;

      %time is a bit more awkward. This is horrible, but it is what it is.
      Working.date = repmat(Working.date,[1,size(Working.time,1)])';
      year = floor(Working.date./1000);
      dn   = Working.date - year.*1000;
      s    = Working.time./1000; s(s < 0) = NaN;
      Store.Time = datenum(year,1,dn,1,1,s)';

      %store in main repository, and tidy
      Data = cat_struct(Data,Store,1);

    end; clear  dn id iEl iVar OldData s Store year working

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SOFIE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'SOFIE';
    for DayNumber=floor(min(Input.TimeRange)):1:floor(max(Input.TimeRange));
      %work out year and day number and hence filepath
      [y,~,~] = datevec(DayNumber); dn = date2doy(DayNumber);
      File = wildcardsearch(InstInfo.Path,['Level2_',sprintf('%04d',y),sprintf('%03d',dn)]);
      if numel(File) == 0; clear y dn File;continue; end

      %load file
      Working = rCDF(File{1});
      

      %get variables
      Store = struct();
      for iVar=1:1:numel(Vars)
        switch Vars{iVar}
          case 'Temp'; Store.Temp  = Working.Temperature';
          case 'Lat';  Store.Lat   = Working.Latitude_83km;
          case 'Lon';  Store.Lon   = Working.Longitude_83km;
          case 'Alt';  Store.Alt   = Working.Altitude';
          case 'Pres'; Store.Pres  = Working.Pressure';
          case 'Time'; Store.Time  = DayNumber+Working.Time_UT./24;
          otherwise;  
            try;   Store.(Vars{iVar}) = Working.(Vars{iVar}); 
            catch; disp(['Variable ',Vars{iVar},' not found, terminating']); return
            end
        end      
      end
      clear iVar Working y dn

      %reshape
      NLevs = numel(Store.Alt); NProfs = numel(Store.Lat);
      Store.Lat  = repmat(Store.Lat, 1,NLevs);
      Store.Lon  = repmat(Store.Lon, 1,NLevs);
      Store.Time = repmat(Store.Time,1,NLevs);
      Store.Alt  = repmat(Store.Alt, NProfs,1);
      clear NLevs NProfs


      %store in main repository, and tidy
      Data = cat_struct(Data,Store,1);



    end; clear DayNumber
  otherwise

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% invalid instrument!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    disp(['Instrument ',Input.Instrument,' not currently handled by this function, terminating'])
    return

 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the data to chosen height scale
%  and make lons -180 to 180
%  select actual time range specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time range: do first to reduce size for later steps
%do at profile level to avoid breaking profiles
Data = reduce_struct(Data,inrange(nanmean(Data.Time,2),Input.TimeRange),[],1);


if Input.OriginalZ == false

  %height scale
  sz = [size(Data.Lat,1),numel(Input.HeightScale)];
  Data2 = struct();
  for iVar=1:1:numel(Vars);
    a = NaN(sz);
    b = Data.(Vars{iVar});
    for iProf=1:1:sz(1)

      %some extra handling here to deal with bad data, but all we're actualy doing is linear interpolation
      Good = find(isfinite(Data.Alt(iProf,:)) ~=0);
      [~,uidx] = unique(Data.Alt(iProf,:));
      Good = intersect(Good,uidx);
      if numel(Good) > 2;  a(iProf,:) = interp1(Data.Alt(iProf,Good),b(iProf,Good),Input.HeightScale); end

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
%%%%%%%%%%%%%%%%%%

if Input.KeepOutliers == 0;

  %create one big array of badness, then remove from ALL variables in one pass
  Bad = [];

  %latitude and longitude have physical limits
  Bad = [Bad;find(Data.Lat <  -90 | Data.Lat >  90 | Data.Lon < -180 | Data.Lon > 180)];
  
  %time and altitude should always be in the specified range
  Bad = [Bad;find(Data.Time < min(Input.TimeRange  ) | Data.Time > max(Input.TimeRange  ))];
  if Input.OriginalZ == false;  Bad = [Bad;find(Data.Alt  < min(Input.HeightScale) | Data.Alt  > max(Input.HeightScale))]; end

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

if min(Input.LatRange) ~=  -90 || max(Input.LatRange) ~=  90 ...
 | min(Input.LonRange) ~= -180 || max(Input.LonRange) ~= 180

  %find the min and max of lon and lat for each profile
  Limits = [min(Data.Lon,[],2),max(Data.Lon,[],2), ...
            min(Data.Lat,[],2),max(Data.Lat,[],2)];

  %generate a list of profiles where any point falls inside the limits
  Bad = [];
  Bad = [Bad;find(Limits(:,1) > max(Input.LonRange))];
  Bad = [Bad;find(Limits(:,2) < min(Input.LonRange))];
  Bad = [Bad;find(Limits(:,3) > max(Input.LatRange))];
  Bad = [Bad;find(Limits(:,4) < min(Input.LatRange))];

  Good = 1:1:size(Limits,1); Good(Bad) = [];
  Data = reduce_struct(Data,Good,[],1);

end



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
