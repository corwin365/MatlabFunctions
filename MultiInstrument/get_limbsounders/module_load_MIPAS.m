function [Data,FileList] = module_load_MIPAS(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from MIPAS for get_limbsounders()
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

  %work out filepath
  [y,m,d] = datevec(DayNumber);
  File = wildcardsearch(InstInfo.Path,['_',sprintf('%04d',y),sprintf('%02d',m),sprintf('%02d',d),'_nom']);
  if numel(File) == 0; clear y dn File; continue; end

  %load data and pull out vars
  Working = rCDF(File{1});

  %store file information
  FileCount = FileCount+1;
  f = strsplit(File{1},filesep); FileList{end+1} = f{end}; clear f;

  %get variables
  Store = struct();
  for iVar=1:1:numel(Vars)
    switch Vars{iVar}
      case 'Temp';       Store.Temp       = Working.TEM';
      case 'Lat';        Store.Lat        = Working.latitude;
      case 'Lon';        Store.Lon        = Working.longitude;
      case 'Alt';        Store.Alt        = Working.hgt';
      case 'Pres';       Store.Pres       = Working.PRE';
      case 'Time';       Store.Time       = DayNumber+Working.time./24;
      case 'SourceProf'; Store.SourceProf = (1:1:numel(Store.Lat))';
      case 'SourceFile'; Store.SourceFile   = ones(size(Store.SourceProf)).*FileCount;
      otherwise;
        try;   Store.(Vars{iVar}) = Working.(Vars{iVar});
        catch; disp(['Variable ',Vars{iVar},' not found, terminating']); return
        end
    end
  end
  clear iVar Working y dn

  %reshape
  NLevs = size(Store.Alt,2); NProfs = numel(Store.Lat);
  Store.Lat        = repmat(Store.Lat,       1,NLevs);
  Store.Lon        = repmat(Store.Lon,       1,NLevs);
  Store.Time       = repmat(Store.Time,      1,NLevs);
  Store.SourceProf = repmat(Store.SourceProf,1,NLevs);
  Store.SourceFile = repmat(Store.SourceFile,1,NLevs);
  clear NLevs NProfs


  %store in main repository
  Data = cat_struct(Data,Store,1);

end; clear DayNimber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return