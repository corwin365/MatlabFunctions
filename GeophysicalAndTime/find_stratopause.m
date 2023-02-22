function [StratopausePressure,ValidMask] = find_stratopause(TimeGrid,LonGrid,LatGrid,Verbose)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%find stratopause in ERA5 data
% method of https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2011JD016893
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.DataDir = [LocalDataDir,'/ERA5'];

%create list of dates to loop over
Settings.DayScale = floor(min(TimeGrid(:))):1:ceil(max(TimeGrid(:)));

%verbosity
if ~exist('Verbose'); Settings.Verbose == 0;
else
  if Verbose == 0; Settings.Verbose = 0; 
  else             Settings.Verbose = 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if the inputs were vectors, meshgrid them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz1 = size(TimeGrid);
sz2 = size(LonGrid);
sz3 = size(LatGrid);

if     sz1(1) == 1 &&     sz2(1) == 1 &&     sz3(1) == 1 ...
&& numel(sz1) == 2 && numel(sz2) == 2 && numel(sz3) == 2;
  %we have vectors - meshgrid
  [TimeGrid,LonGrid,LatGrid] = meshgrid(TimeGrid,LonGrid,LatGrid);
elseif isequal(sz1,sz2) && isequal(sz1,sz3)
  %we have grids - keep them and do nothing
else
  %we have some combination - reject
  disp('For find_stratopause, grid arrays must be provided as either all vectors or as valid 3D grids, stopping')
  Height = NaN;
  Pressure = NaN;
  return
end

clear sz1 sz2 sz3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processing - compute and extract at native resolution of
%ERA5 files used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.DayScale)
  
 
  %find and load the ERA5 data for this day
  %
  %file naming scheme is the one I use, and 
  %allows for both ERA5 and ERA5T. On another
  %system you probably want to change this
  %bit
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %get data
  [yy,~,~] = datevec(Settings.DayScale(iDay));
  dd = date2doy(Settings.DayScale(iDay));
  File = [Settings.DataDir,'/',sprintf('%04d',yy),'/',...
          'era5_',sprintf('%04d',yy),'d',sprintf('%03d',dd),'.nc'];
  if ~exist(File,'file');
    File = [Settings.DataDir,'/',sprintf('%04d',yy),'/',...
            'era5t_',sprintf('%04d',yy),'d',sprintf('%03d',dd),'.nc'];
    if ~exist(File,'file'); continue; end
  end
  Data = rCDF(File);

  %compute pressure. We'll ignore the effects on lnsp on the ERA5
  %column spacing as we're at the stratopause so this effect is
  %negligible
  [Prs, OldZ] = ecmwf_prs_v3(numel(Data.level));
  clear yy dd File  

  %store geoloc for later
  if ~exist('Results');
    Lat = Data.latitude;
    Lon = Data.longitude;
  end
  
  %extract temperature
  T = Data.t;
  clear Data
  
  %make the data ascend vertically
  [~,idx] = sort(Prs,'desc');
  Prs = Prs(idx); OldZ = OldZ(idx);
  T = T(:,:,idx,:);
  clear idx
  
  %same for latitude
  [~,idx] = sort(Lat,'asc');
  Lat = Lat(idx);
  T = T(:,idx,:,:);
  clear idx


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% find stratopause
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  %interpolate the data to a regular 1km height grid between
  %25 and 80km altitude
  NewZ = 25:1:80;
  sz = size(T);  
  T = reshape(permute(T,[1,2,4,3]),sz(1)*sz(2)*sz(4),sz(3));
  NewPrs = interp1(OldZ,Prs,NewZ);
  T = interp1(OldZ,T',NewZ)';

  %smooth with an 11km boxcar
  Ts = smoothn(T,[1,11]);

  %find maximum in each profile
  [~,idx] = max(Ts,[],2);

  %check 5 levels above and below:
    %5 levels above must have -ve lapse rate
    %5 levels below must have +ve lapse rate
  dTdZ = diff(Ts,1,2);
  dTdZ = cat(2,dTdZ(:,1),dTdZ); %so points line up, rather than half-levels

  Passed = zeros(sz(1)*sz(2)*sz(4),1);
  for iProf=1:1:size(dTdZ,1)
    
    Above = idx(iProf)+1:1:idx(iProf)+5; Above = Above(Above > 0 & Above < size(NewZ,2));
    Below = idx(iProf)-5:1:idx(iProf)-1; Below = Below(Below > 0 & Below < size(NewZ,2));    
    
    Above = -dTdZ(iProf,Above); Below = dTdZ(iProf,Below); %note - sign on Above
     
    if min(Above) > 0 & min(Below) > 0;
      %label profile as to use
      Passed(iProf) = 1;

      %also remove anything outside +/- 15 km from peak
      T(iProf,NewZ < NewZ(idx(iProf))-15) = NaN;
      T(iProf,NewZ > NewZ(idx(iProf))+15) = NaN;
    end
    
  end; clear iProf Above Below

  %for all profiles that pass the check, find maximum in unsmoothed data
  Stratopause = Passed.*NaN;
  [~,idx] = max(T(Passed == 1,:),[],2);
  Stratopause(Passed == 1) = h2p(NewZ(idx));
  a = reshape(Stratopause,[sz(1),sz(2),sz(4)]);


  %store
  if ~exist('StratoStore'); StratoStore = repmat(Stratopause,[1,numel(Settings.DayScale)]).*NaN; end
  StratoStore(:,iDay) = Stratopause;


  
  if Settings.Verbose == 1; disp(['Stratopause computed for ',datestr(Settings.DayScale(iDay))]); end
end; clear iDay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate onto the requested spatiotemporal grid
%interpolation is linear in space/time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%there are some NaNs but they cover small areas where conditions weren't met above 
StratoStore = inpaint_nans(StratoStore);



%make time linear
StratoStore = reshape(StratoStore,[sz(1),sz(2),sz(4).*numel(Settings.DayScale)]);
thours = linspace(0,24,sz(4)+1); thours=thours(1:end-1)./24;
[td,th] = meshgrid(Settings.DayScale,thours);
t = td(:)+th(:); clear td th thours

%create interpolant (remember ERA5 lat is monotonically descending)
I = griddedInterpolant({Lon,Lat,t},StratoStore);

%and interpolate
StratopausePressure = permute(I(LatGrid,LonGrid,TimeGrid),[2,3,1]); %dimension order is same as inputs

%occasionally the interpoolation produces some small negative numbers. These are unphysical.
StratopausePressure  = abs(StratopausePressure);
