function [IAGOS,Error] = prep_iagos(FilePath,varargin)


%load an IAGOS netCDF file and add:
%1. a matlab time index
%2. fill any bad data with NaNs
%3. interpolate to constant time sampling (units are in matlabtime, i.e. days)
%4. work out spatial distance between each pair of points



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create input parser
%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;


%inputs - required
%%%%%%%%%%%%%%%%%%%
addRequired(p,'FilePath',@ischar);


%inputs - fully optional
%%%%%%%%%%%%%%%%%%%%%%%%%

CheckPositive = @(x) validateattributes(x,{'numeric'},{'>',0});

%sampling rate to interpolate to (days)
addParameter(p,'SamplingRate',1./24./60./60.*5,CheckPositive);  %5-second sampling

%max height change in a cruise
addParameter(p,'CruiseDz',100,CheckPositive);  %100m

%assume we want to apply the error flags unless stated otherwise
addParameter(p,'ApplyFlags',true,@islogical);

%assume we want to sample evenly in time rather than space
addParameter(p,'SpaceSamplingRate',false,@islogical);

%time window for the above change (in units of SamplingRate, above)
CheckWindow = @(x) validateattributes(x,{'numeric'},{'>',0,'odd'});
addParameter(p,'CruiseWindow',25,CheckWindow);  %25 time steps


%parse the inputs, and tidy up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse inputs
parse(p,FilePath,varargin{:})

%pull out the contents into struct "Inputs", used throughout rest of routine
Input = p.Results;

%tidy up
clear CheckWindow CheckPositive p varargin


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check file exists, and fail if not
if ~exist(FilePath)
  Error = 'No file found';
  IAGOS = [];
  return
end

IAGOS = rCDF(FilePath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute and add timestamps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get the metadata variable that gives the time basis for this file
Variables = {IAGOS.MetaData.Variables(:).Name};
VarIdx = find(strcmp(Variables,'UTC_time') ~= 0);

%go inside this variable and find the units attribute
Variables = {IAGOS.MetaData.Variables(VarIdx).Attributes(:).Name};
VarIdx2 = find(strcmp(Variables,'units') ~= 0);
Units = IAGOS.MetaData.Variables(VarIdx).Attributes(VarIdx2).Value;
clear Variables VarIdx VarIdx2

%hence, convert the dates
IAGOS.Time = datenum(Units(15:end),'yyyy-mm-dd HH:MM:SS') + IAGOS.UTC_time./60./60./24;
clear Units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% deduplicate lats and lons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%some records (especially early ones) have quite imprecise lat and lon
%data, which results in points which are at the "same place"

%let's fix them, by linearly inteprolating to the next unique lat/lon

[~,llidx] = unique(IAGOS.lon + 1000.*IAGOS.lat);
llidx = sort(llidx,'asc'); %this is necessary for the internal if loop to work correctly
if numel(llidx) ~= unique(IAGOS.lon);
  
  for iX = 1:1:numel(llidx)
    
    %find all points with this lat/lon pair
    Lon = IAGOS.lon(llidx(iX));
    Lat = IAGOS.lat(llidx(iX));
    
    idx = find(IAGOS.lon == Lon & IAGOS.lat == Lat);
    
    %hence, interpolate the lats and lons that are missing
    %obviously we cannot do this if it is the end of the flight, so check
    %this first
    
    if max(idx) < numel(IAGOS.lon)
      NextLon = IAGOS.lon(max(idx)+1);
      NextLat = IAGOS.lat(max(idx)+1);
      
      NewLon = linspace(Lon,NextLon,numel(idx));
      NewLat = linspace(Lat,NextLat,numel(idx));
      
      IAGOS.lon(idx) = NewLon;
      IAGOS.lat(idx) = NewLat;
      
    else
      %extrapolate if we're at the end
      %this isn't great, but we won't be using landing data anyway and it
      %should be very rare
      x = min(idx)-5:min(idx)-1;
      xq = idx;
      IAGOS.lon(idx) = interp1(x,IAGOS.lon(x),xq,'linear','extrap');
      IAGOS.lat(idx) = interp1(x,IAGOS.lat(x),xq,'linear','extrap');
      clear x xq
    end
    

    
  end
  
end
clear llidx iX Lon Lat idx NextLon NextLat NewLon NewLat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NaNify actively bad data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iVar=1:1:numel(IAGOS.MetaData.Variables)
  
  %find the fill value
  Attrib= {IAGOS.MetaData.Variables(iVar).Attributes(:).Name};
  idx = find(strcmp(Attrib,'missing_value') ~= 0);
  if numel(idx) == 0; continue; end
  FillVal = IAGOS.MetaData.Variables(iVar).Attributes(idx).Value;
  
  %convert fill value to NaN
  Var = IAGOS.(IAGOS.MetaData.Variables(iVar).Name);
  Var(Var == FillVal) = NaN;
  IAGOS.(IAGOS.MetaData.Variables(iVar).Name) = Var;
  
end
clear Attrib idx FillVal iVar Var

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% apply error flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Input.ApplyFlags == 1;
  
  %list of vars which have validation flags
  Vars = {'baro_alt_AC','radio_alt_AC','gps_alt_AC','air_press_AC',...
          'air_speed_AC','ground_speed_AC','air_temp_AC','air_stag_temp_AC', ...
          'wind_dir_AC','wind_speed_AC','zon_wind_AC','mer_wind_AC', ...
          'air_temp_PM','air_stag_temp_PM'};
  
  for iVar=1:1:numel(Vars)
    
    Flags = IAGOS.([Vars{iVar},'_val']);
    if nansum(Flags) == 0; continue; end %no problem
    
    %identify bad data and set it to NaN
    Bad = find(Flags ~= 0 & ~isnan(Flags));
    Var = IAGOS.(Vars{iVar});
    Var(Bad) = NaN;
    IAGOS.(Vars{iVar}) = Var;
    clear Bad Var
    
    
    
  end; clear iVar Vars
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate to constant time sampling or constant space sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Input.SpaceSamplingRate ~= 1;
  %%%%%%%%%%%%%%%%%%
  %time resampling
  %%%%%%%%%%%%%%%%%%
    
  OldTime = IAGOS.Time;
  NewTime = min(IAGOS.Time): Input.SamplingRate : max(IAGOS.Time);
  
  Fields = fieldnames(IAGOS);
  for iField=1:1:numel(Fields)
    if strcmp(Fields{iField},'MetaData'); continue; end
    
    IAGOS.(Fields{iField}) = interp1(OldTime,IAGOS.(Fields{iField}),NewTime);
    
  end
  
  IAGOS.OriginalTime = OldTime;

  clear OldTime NewTime iField Fields

else
  %%%%%%%%%%%%%%%%%%
  %space resampling
  %%%%%%%%%%%%%%%%%%
  
  %work out distance between each point
  x = [IAGOS.lat,IAGOS.lon];
  y = circshift(x,1,1);
  dx = nph_haversine(x,y);
  dx([1,end]); %assume first point spacing is same as second, as can't compute otherwise
  
  %and hence stepwise distance along-track
  dxS = dx; for iX=2:1:numel(dxS); dxS(iX) = dxS(iX-1)+dxS(iX); end
  
  %work out new spacing
  dx2 = 0:Input.SamplingRate:max(dxS);
  
  %interpolate
  Fields = fieldnames(IAGOS);  
  for iField=1:1:numel(Fields)
    if strcmp(Fields{iField},'MetaData'); continue; end
    
    IAGOS.(Fields{iField}) = interp1(dxS,IAGOS.(Fields{iField}),dx2);
    
  end
  
  clear iField Fields dxS dx2
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% travel distance (in km)
%% using lat/lon, but cross-checked with speed/time for a few hundred flights and look pretty similar
%if we did space sampling then this is just repeating what we did and
%storing the answer, but the runtime is very small so this is fine, and
%makes the computational logic cleaner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = [IAGOS.lat;IAGOS.lon];
b = circshift(a,1,2);

IAGOS.dx = nph_haversine(a',b');
IAGOS.dx([1,end])  = 0;
clear a b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify "cruises" in the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


         
%derivative of altitude
dz = [diff(IAGOS.baro_alt_AC),0];


%produce a running windowed *sum* of this value
dz = smooth(dz,Input.CruiseWindow).*Input.CruiseWindow;

%take absolute value, as we're looking for discontinuities, not signs
dz = abs(dz);

%create a binary mask - 0 isgood, 1 is bad (during or near a height change)
Track = zeros(size(dz)); Track(dz >= Input.CruiseDz) = 1;

%and then divide the good bits into separate sequences seperated by
%the bad bits. Remember that the plane will START AND END in a BAD
%state, as it has to get to and from ground level.
Starts = find(diff(Track) == -1); Starts = Starts(1:end-1);%last one is the end of the final descent
Ends   = find(diff(Track) ==  1);
if numel(Starts) ~= numel(Ends); Ends = Ends(2:end); end %cross-check for sanity, as some have a bobble up and down at the start that triggers the end detector

%these are the starts and ends of the cruises
%split the data up
csize = Ends-Starts+1;

Cruises = NaN(numel(Starts),max(csize));
for iCruise=1:1:numel(Starts);
  Cruises(iCruise,1:csize(iCruise)) = Starts(iCruise):1:Ends(iCruise);
end
clear iCruise Starts Ends Track dz csize

IAGOS.Cruises = Cruises; clear Cruises

end

