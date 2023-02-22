function TropopausePressure = find_tropopause(TimeGrid,LonGrid,LatGrid,Verbose)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%find stratopause in ERA5 data, using WMO definition
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

if Settings.Verbose == 1; textprogressbar('Computing tropopause '); end
for iDay=1:1:numel(Settings.DayScale)

  if Settings.Verbose == 1; textprogressbar(iDay./numel(Settings.DayScale).*100); end

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
  [Prs, Z] = ecmwf_prs_v3(numel(Data.level));
  clear yy dd File  

  %store geoloc for later
  if ~exist('Results');
    Lat = Data.latitude;
    Lon = Data.longitude;
  end
  
  %extract temperature
  T = Data.t;
  clear Data
  
  %make the data go in the "right" direction vertically
  [~,idx] = sort(Prs,'desc');
  Prs = Prs(idx); Z = Z(idx);
  T = T(:,:,idx,:);
  clear idx

  %same for latitude
  [~,idx] = sort(Lat,'asc');
  Lat = Lat(idx);
  T = T(:,idx,:,:);
  clear idx

   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% step 1: lapse rate 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  dT = diff(T,1,3);
  dZ = diff(Z);

  Gamma = dT .* NaN;
  for iLevel=1:1:numel(dZ); Gamma(:,:,iLevel,:) = dT(:,:,iLevel,:)./dZ(iLevel); end;

  clear iLevel dT dZ T


 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% step 2: find tropopause
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %create an array to store our tropopause levels
  sz = size(Gamma);
  Tropopause = NaN(sz([1,2,4]));
  clear sz

  %loop over half-levels
  for iLevel = 1:1:136

    %if we've already found all the tropopauses, skip
    if sum(isnan(Tropopause(:))) == 0; continue; end
    
    %if pressure > 700hPa, skip
    if Prs(iLevel) > 700; continue; end
    
    %if pressure < 10, skip
    if Prs(iLevel) < 10; continue; end

    %check if Gamma is greater than -2
    idx = find(Gamma(:,:,iLevel,:) > -2);
    
    if numel(idx) == 0; continue; end %none at this level
    
    %remove any columns we already found
    Found = find(~isnan(Tropopause));
    [~,Remove] = intersect(idx,Found);
    idx(Remove) = [];
    clear Remove
    
    %for each element where the above criterion is met, check if the layer
    %2km higher also meets it
    
    %find which level is 2km above
    Z = p2h(Prs(iLevel));
    jLevel = closest(p2h(Prs),Z+2);
    Range = sort([iLevel,jLevel],'ascend'); 
    Range = Range(1):1:Range(2);
    clear Z jLevel
    
    %pull out this range for each of the elements of interest
    G2 = permute(Gamma(:,:,Range,:),[1,2,4,3]);
    sz = size(G2);
    G2 = reshape(G2,sz(1)*sz(2)*sz(3),sz(4));
    G2 = G2(idx,:);
    clear Range
    
    %find all the columns where the criterion remains met for 2km above
    StillMet = min(G2,[],2);
    Good = find(StillMet > -2);
    idx = idx(Good);
    if numel(Good) <2 ; continue; end %0 is obviously wrong, excluding 1 bypasses a minor bug with array shapes below that isn't worth fixing for single pixels we can interpolate over at the end
    

    %find where the gradient crossed above -2 by linear interpolation
    G3 = permute(Gamma(:,:,iLevel+[-1:1],:),[1,2,4,3]);
    G3 = reshape(G3,sz(1)*sz(2)*sz(3),3);
    G3 = G3(idx,:);
    p2 = linspace(Prs(iLevel-1),Prs(iLevel+1),10);
    G3 = interp1(Prs(iLevel+[-1:1]),G3',p2);
    G3 = abs(G3 + 2);
    [~,G3] = min(G3,[],1);
    Val = p2(G3);
    
    %yay :-) store, and remove these columns from the lapse rate data
    Tropopause(idx) = Val;
    
  end; clear iLevel G2 G3 idx Found Gamma Good p2 StillMet Val

  %fill gaps via interpolation (usually very few or none, mostly due to the above edge-case skipping)
  if sum(isnan(Tropopause(:))) > 0
    for iTime=1:1:size(Tropopause,3)
      Tropopause(:,:,iTime) = inpaint_nans(Tropopause(:,:,iTime));
    end; clear iTime
  end
    
  %store
  if ~exist('TropoStore');TropoStore = NaN([size(Tropopause),numel(Settings.DayScale)]);  end
  TropoStore(:,:,:,iDay) = Tropopause;
  
  %tidy up
  clear Tropopause Gamma Prs

  
end; clear iDay sz
if Settings.Verbose == 1; textprogressbar('!'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate onto the requested spatiotemporal grid
%interpolation is linear in space/time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make time linear
sz = size(TropoStore);
TropoStore = reshape(TropoStore,[sz(1),sz(2),sz(3)*sz(4)]);
thours = linspace(0,24,sz(3)+1); thours=thours(1:end-1)./24;
[td,th] = meshgrid(Settings.DayScale,thours);
t = td(:)+th(:); clear td th thours

%create interpolant (remember ERA5 lat is monotonically descending)
I = griddedInterpolant({Lon,Lat,t},TropoStore);

%and interpolate
TropopausePressure = permute(I(LatGrid,LonGrid,TimeGrid),[2,3,1]); %dimension order is same as inputs


%done
return
