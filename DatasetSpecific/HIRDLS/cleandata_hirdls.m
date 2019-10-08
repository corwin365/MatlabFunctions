function HirdlsData = cleandata_hirdls(HirdlsData,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clean HIRDLS data, using a set of flags to specify exactly how
%assume all methods if no flags specified
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%04/MAY/2014
%
%inputs
%---------
%
%HirdlsData - struct containing HIRDLS data
%varargin (optional) - flags showing methods to NOT apply, as follows:
%                      - OmitDateConvert - converts date format to matlab
%                      - OmitGridGeolocation - puts all geolocation on T grid
%                      - OmitAddPres - restores pressure data, since I forgot to put it in the files...                     
%                      - OmitCamelCase - CamelCases all the variables               
%                      - OmitDouble - convert all numbers to double format  
%                      - OmitOutRem - removes outliers (T<100 or T>400)  
%                      - OmitInterp - inteprolates data to a uniform 1km vertical scale. removes NaNs.
%                      - OmitHResScale - corrects for lat/lon difference as function of height
%                    -also has options (positive choices):
%                      - FastInterp - assumes height profile is same for every profile to height-interpolate rapidly
%
%*****NOTE THAT OMITTING SOME STEPS MAY WELL BREAK OTHERS - STRONGLY RECOMMEND USING ALL, I.E. NO FLAGS)
%
%outputs
%---------
%
% HirdlsData - cleaned data to return
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%double() variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitDouble',varargin);
  HirdlsData.TEMP       = double(HirdlsData.TEMP  );
  HirdlsData.HEIGHT     = double(HirdlsData.HEIGHT);
  HirdlsData.LAT        = double(HirdlsData.LAT   ); 
  HirdlsData.LON        = double(HirdlsData.LON   ); 
  HirdlsData.TIME       = double(HirdlsData.TIME  ); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%put geolocation on same grid as temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitGridGeolocation',varargin);

  NewLat  = HirdlsData.TEMP .* NaN;
  NewLon  = NewLat;
  NewTime = NewLat;
  
  for i=1:1:numel(NewLat(:,1));
    NewLat( i,:) = HirdlsData.LAT;
    NewLon( i,:) = HirdlsData.LON;
    NewTime(i,:) = HirdlsData.TIME;
  end

  HirdlsData.LAT  = NewLat;  clear NewLat;
  HirdlsData.LON  = NewLon;  clear NewLon;
  HirdlsData.TIME = NewTime; clear NewTime;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert date format to matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitDateConvert',varargin);
  %HIRDLS format is TAI93
  HirdlsData.MatlabTime = (datenum(1993,1,1, 0, 0, 0) + HirdlsData.TIME./86400);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%restore pressure axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitAddPres',varargin);
  %HIRDLS L2 is computed on a standard P grid. This is it. (121 levels)
  Prs = [1000.00,908.518,825.404,749.894,681.292,618.966,562.341,510.897,464.159,...
       421.697,383.119,348.070,316.228,287.298,261.016,237.137,215.443,195.734,...
       177.828,161.560,146.780,133.352,121.153,110.069,100.000,90.8517,82.5404,...
       74.9894,68.1292,61.8966,56.2341,51.0897,46.4159,42.1697,38.3119,34.8070,...
       31.6228,28.7298,26.1016,23.7137,21.5443,19.5734,17.7828,16.1560,14.6780,...
       13.3352,12.1153,11.0069,10.0000,9.08517,8.25404,7.49894,6.81292,6.18966,...
       5.62341,5.10897,4.64159,4.21696,3.83119,3.48070,3.16228,2.87298,2.61016,...
       2.37137,2.15443,1.95734,1.77828,1.61560,1.46780,1.33352,1.21153,1.10069,...
       1.00000,0.908517,0.825404,0.749894,0.681292,0.618966,0.562341,0.510897,...
       0.464159,0.421697,0.383119,0.348070,0.316228,0.287298,0.261016,0.237137,...
       0.215443,0.195734,0.177828,0.161560,0.146780,0.133352,0.121153,0.110069,...
       0.100000,0.0908517,0.0825404,0.0749894,0.0681292,0.0618966,0.0562341,...
       0.0510897,0.0464159,0.0421696,0.0383118,0.0348070,0.0316228,0.0287298,...
       0.0261016,0.0237137,0.0215443,0.0195734,0.0177828,0.0161560,0.0146780,...
       0.0133352,0.0121153,0.0110069,0.0100000]; 
  
  HirdlsData.Prs = HirdlsData.TEMP.*NaN;
  for i=1:1:numel(HirdlsData.TEMP(1,:)); HirdlsData.Prs(:,i) = Prs; end
  clear Prs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove outliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitOutRem',varargin);
  HirdlsData.TEMP(HirdlsData.TEMP <  50) = NaN;
  HirdlsData.TEMP(HirdlsData.TEMP > 500) = NaN;  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inteprolate to 1km scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitInterp',varargin);
  
  %find top and bottom, to nearest km, in m
  LowestHeight  = floor(nanmin(HirdlsData.HEIGHT(HirdlsData.HEIGHT > 0))./1000).*1000;
  HighestHeight = ceil( nanmax(HirdlsData.HEIGHT(HirdlsData.HEIGHT > 0))./1000).*1000;
  
  %find number of new levels, and NProfiles
  NNewLevels = (HighestHeight-LowestHeight)./1000;
  NProfiles  = numel(HirdlsData.TEMP(1,:));
  
  %produce new variables
  NewTemp       = NaN(NNewLevels,NProfiles);
  NewHeight     = NaN(NNewLevels,NProfiles);
  NewLat        = NaN(NNewLevels,NProfiles);
  NewLon        = NaN(NNewLevels,NProfiles);
  NewTime       = NaN(NNewLevels,NProfiles);
  NewPrs        = NaN(NNewLevels,NProfiles);
  NewMatlabTime = NaN(NNewLevels,NProfiles);
  
  
  %interpolate  
  if ~ismember('FastInterp',varargin);
    %interpolate every profile properly
    
    for iProfile=1:1:NProfiles;
      A = (1:1:NNewLevels).*1000;
      B = inpaint_nans(squeeze(HirdlsData.HEIGHT(:,iProfile)));
      
      %occasionally get multiple zeros at bottom of profile which mess up the interp
      %just add negative single meters to fix this - not interested in this region anyway
      NZerosAtBottom = numel(find(B == 0));
      if NZerosAtBottom > 1; for i=NZerosAtBottom-1:-1:1; B(i) = B(i+1)-1;end;  end;
      
      %sometimes it's just all bad datas...
      if nansum(B) < 0 ;continue;end
      
      %interpolate!
      NewHeight(    :,iProfile) = A;
      NewTemp(      :,iProfile) = interp1(B,inpaint_nans(squeeze(HirdlsData.TEMP(      :,iProfile))),A);
      NewLat(       :,iProfile) = interp1(B,inpaint_nans(squeeze(HirdlsData.LAT(       :,iProfile))),A);
      NewLon(       :,iProfile) = interp1(B,inpaint_nans(squeeze(HirdlsData.LON(       :,iProfile))),A);
      NewTime(      :,iProfile) = interp1(B,inpaint_nans(squeeze(HirdlsData.TIME(      :,iProfile))),A);
      NewPrs(       :,iProfile) = interp1(B,inpaint_nans(squeeze(HirdlsData.Prs(       :,iProfile))),A);
      NewMatlabTime(:,iProfile) = interp1(B,inpaint_nans(squeeze(HirdlsData.MatlabTime(:,iProfile))),A);
    end
  else
    %assume every profile is the same height-wise, and interpolate accordingly
    %slightly less accurate, but significantly faster
      A = (1:1:NNewLevels).*1000;
      B = inpaint_nans(HirdlsData.HEIGHT(:,1));
      
     %interpolate!
     %%
      for iProfile=1:1:NProfiles;NewHeight(    :,iProfile) = A;end
      NewTemp       = interp1(B,inpaint_nans(squeeze(HirdlsData.TEMP)),      A);
      NewLat        = interp1(B,inpaint_nans(squeeze(HirdlsData.LAT)),       A);
      NewLon        = interp1(B,inpaint_nans(squeeze(HirdlsData.LON)),       A);
      NewTime       = interp1(B,inpaint_nans(squeeze(HirdlsData.TIME)),      A);
      NewPrs        = interp1(B,inpaint_nans(squeeze(HirdlsData.Prs)),       A);
      NewMatlabTime = interp1(B,inpaint_nans(squeeze(HirdlsData.MatlabTime)),A);

  end
  
  %overwrite originals
  HirdlsData.TEMP       = NewTemp;
  HirdlsData.HEIGHT     = NewHeight;
  HirdlsData.LAT        = NewLat;
  HirdlsData.LON        = NewLon;
  HirdlsData.TIME       = NewTime;
  HirdlsData.Prs        = NewPrs;
  HirdlsData.MatlabTime = NewMatlabTime;
  
  %wipe copies and temp variables
  clear NewTemp NewHeight NewLat NewLon NewTime NewPrs NewMatlabTime
  clear A B HighestHeight LowestHeight NNewLevels NProfiles iProfile
  clear NZerosAtBottom
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%correct for lat/lon difference as function of height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitHResScale',varargin);
  
  %from IDL:
%   ;these correction values are HIRDLS SPECIFIC
%   GeoLocationHeight = 30000.   ;m - approx height where lat,lon actually equal the stored value
%   ShiftPerLevel     = 0.8      ;horizontal shift factor, meters horizontal per meter vertical, due to scan
%   ShiftPerLevel     = ShiftPerLevel/1000.*StLevelSize

%   Direction     = L2ScanUpFlag[iProfile]-L2ScanUpFlag[iProfile+1]
%   ZShift        = (StHeightAxis[jLevel]-GeoLocationHeight)
%   LocationShift = Direction*ZShift*ShiftPerLevel
%   DeltaX       -= LocationShift
  
  
  AllNewLats = HirdlsData.LAT;
  AllNewLons = HirdlsData.LON;

  for iProfile=1:1:numel(HirdlsData.TIME(1,:))-1;
  
    %get lat and lon of this profile and the next 
    LatOI = HirdlsData.LAT(1,iProfile:iProfile+1);
    LonOI = HirdlsData.LON(1,iProfile:iProfile+1);
    
    %beware of NaNs
    if ~isfinite(sum(LatOI)) || ~isfinite(sum(LonOI)); continue ; end;
    
    %check if we cross the dateline, and shift the lower value into the 180-360 range if so
    if range(LonOI) > 180;
      [Small,SmallIndex] = min(LonOI);
      LonOI(SmallIndex) = LonOI(SmallIndex)+360;
    end
    
    %compute the bearing of the satellite vector
    Bearing = radtodeg(atan2((LatOI(2)-LatOI(1)),(LonOI(2)-LonOI(1))));
    
    %compute the separation between the profiles, in m 
    SeparationM = distdim(distance(LatOI(1),LonOI(1),LatOI(2),LonOI(2)),'deg','meters');

    %scan up or scan down?
    Direction = HirdlsData.SUF(iProfile)-HirdlsData.SUF(iProfile+1);
    
    %for each level, project the position the appropriate distance
    %along the scan track
    GeolocationHeight = 30000;
    ShiftPerLevel     = 0.8; %horizontal shift factor, meters horizontal per meter vertical, due to scan
    for iLevel=1:1:numel(HirdlsData.TIME(:,1));
     ZShift = HirdlsData.HEIGHT(iLevel,iProfile)-GeolocationHeight;
     LocationShift = Direction*ZShift*ShiftPerLevel;
     FracLocationShift = LocationShift./SeparationM;
      
     %shift values appropirate distance along track
     NewLon = LonOI(1) + ((LonOI(2)-LonOI(1)).*FracLocationShift);
     NewLat = LatOI(1) + ((LatOI(2)-LatOI(1)).*FracLocationShift);
     
     AllNewLats(iLevel,iProfile) = NewLat;
     AllNewLons(iLevel,iProfile) = NewLon;
    end

  end
  
  %return to original vars
  HirdlsData.LAT = AllNewLats;
  HirdlsData.LON = AllNewLons;
  
  clear AllNewLats AllNewLons Bearing Direction FracLocationShift
  clear GeolocationHeight LatOI LonOI NewLat NewLon SeparationM ShiftPerLevel
  clear ZShift i iLevel iProfile LocationShift 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%camelcase variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitCamelCase',varargin);
  HirdlsData.Temp   = HirdlsData.TEMP;   HirdlsData = rmfield(HirdlsData,'TEMP');
  HirdlsData.Height = HirdlsData.HEIGHT; HirdlsData = rmfield(HirdlsData,'HEIGHT');
  HirdlsData.Lat    = HirdlsData.LAT;    HirdlsData = rmfield(HirdlsData,'LAT');
  HirdlsData.Lon    = HirdlsData.LON;    HirdlsData = rmfield(HirdlsData,'LON');
  HirdlsData.Time   = HirdlsData.TIME;   HirdlsData = rmfield(HirdlsData,'TIME');
end

