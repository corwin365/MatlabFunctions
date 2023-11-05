function [Data,FileList] = module_load_SOFIE(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from SOFIE for get_limbsounders()
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
  File = wildcardsearch(InstInfo.Path,['Level2_',sprintf('%04d',y),sprintf('%03d',dn)]);
  if numel(File) == 0; clear y dn File;continue; end

  %load file
  Working = rCDF(File{1});

  %store file information
  FileCount = FileCount+1;
  f = strsplit(File{1},filesep); FileList{end+1} = f{end}; clear f;

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
      case 'SourceProf'; Store.SourceProf = repmat((1:1:size(Store.Lat,1))',1,size(Store.Alt,2));
      case 'SourceFile'; Store.SourceFile = ones(size(Store.SourceProf)).*FileCount;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return