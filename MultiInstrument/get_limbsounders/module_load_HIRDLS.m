function [Data,FileList] = module_load_HIRDLS(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from HIRDLS for get_limbsounders()
%
%Corwin Wright, c.wright@bath.ac.uk, 05/NOV/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%empty structure to store data
Data = struct();
for iVar=1:1:numel(Vars); Data.(Vars{iVar}) = []; end
FileCount = 0;
FileList = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for DayNumber=floor(min(Settings.TimeRange)):1:floor(max(Settings.TimeRange));

  %work out year and day number and hence filepath
  [y,~,~] = datevec(DayNumber); dn = date2doy(DayNumber);
  File = wildcardsearch(InstInfo.Path,['_',sprintf('%04d',y),'d',sprintf('%03d',dn)]);
  if numel(File) == 0; clear y dn File; continue; end

  %store file information
  FileCount = FileCount+1;
  f = strsplit(File{1},filesep); FileList{end+1} = f{end}; clear f;


  %load variables we need
  Store = struct();
  for iVar=1:1:numel(Vars)
    switch Vars{iVar}
      case 'Temp';       Store.Temp         = get_HIRDLS(File{1},'Temperature')';
      case 'Lat';        Store.Lat          = get_HIRDLS(File{1},'Latitude');
      case 'Lon';        Store.Lon          = get_HIRDLS(File{1},'Longitude');
      case 'Alt';        Store.Alt          = get_HIRDLS(File{1},'Altitude')'./1000;
      case 'Pres';       Store.Pres         = get_HIRDLS(File{1},'Pressure');
      case 'Time';       Store.Time         = datenum(1993,1,1,0,0,get_HIRDLS(File{1},'Time'));
      case 'SourceProf'; Store.SourceProf   = single((1:1:numel(Store.Lat))');
      case 'SourceFile'; Store.SourceFile   = ones(size(Store.SourceProf)).*FileCount;
      otherwise;
        try;   Store.(Vars{iVar}) = get_HIRDLS(File{1},Vars{iVar});
        catch; disp(['Variable ',Vars{iVar},' not found, terminating']); return
        end
    end
  end

  %reshape 1d variables
  Store.Lat        = repmat(Store.Lat,         [1,size(Store.Temp,2)]);
  Store.Lon        = repmat(Store.Lon,         [1,size(Store.Temp,2)]);
  Store.Time       = repmat(Store.Time,        [1,size(Store.Temp,2)]);
  Store.Pres       = repmat(Store.Pres',       [size(Store.Temp,1),1]);
  Store.SourceProf = repmat(Store.SourceProf,  [1,size(Store.Temp,2)]);
  Store.SourceFile = repmat(Store.SourceFile,  [1,size(Store.Temp,2)]);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return