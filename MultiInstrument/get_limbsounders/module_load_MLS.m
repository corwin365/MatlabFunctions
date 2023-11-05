function [Data,FileList] = module_load_MLS(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from MLS for get_limbsounders()
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

  %load data for day
  Working = get_MLS(File{1},'Temperature');

  %store file information
  FileCount = FileCount+1;
  f = strsplit(File{1},filesep); FileList{end+1} = f{end}; clear f;

  %load variables we need
  Store = struct();
  for iVar=1:1:numel(Vars)
    switch Vars{iVar}
      case 'Temp';       Store.Temp         = Working.L2gpValue';
      case 'Lat';        Store.Lat          = Working.Latitude;
      case 'Lon';        Store.Lon          = Working.Longitude;
      case 'Alt';        Store.Alt          = p2h(Working.Pressure);
      case 'Pres';       Store.Pres         = Working.Pressure;
      case 'Time';       Store.Time         = datenum(1993,1,1,0,0,Working.Time);
      case 'SourceProf'; Store.SourceProf   = (1:1:numel(Store.Lat))';
      case 'SourceFile'; Store.SourceFile   = ones(size(Store.SourceProf)).*FileCount;
      otherwise;
        try;   Store.(Vars{iVar}) = Working.(Vars{iVar});
        catch; disp(['Variable ',Vars{iVar},' not found, terminating']); return
        end
    end
  end

  %reshape 1d variables
  Store.Lat        = repmat(Store.Lat,  [1,size(Store.Temp,2)]);
  Store.Lon        = repmat(Store.Lon,  [1,size(Store.Temp,2)]);
  Store.Time       = repmat(Store.Time, [1,size(Store.Temp,2)]);
  Store.SourceProf = repmat(Store.SourceProf,  [1,size(Store.Temp,2)]);
  Store.SourceFile = repmat(Store.SourceFile,  [1,size(Store.Temp,2)]);
  Store.Pres       = repmat(Store.Pres',[size(Store.Temp,1),1]);
  Store.Alt        = repmat(Store.Alt', [size(Store.Temp,1),1]);
  for iField=1:1:numel(Settings.AdditionalVars);
    F = Store.(Settings.AdditionalVars{iField});
    if     numel(F) == size(Store.Temp,2); F = repmat(F',[size(Store.Temp,1),1]);
    elseif numel(F) == size(Store.Temp,1); F = repmat(F, [1,size(Store.Temp,2)]);
    end
    Store.(Settings.AdditionalVars{iField}) = F;
  end

  %store in main repository
  Data = cat_struct(Data,Store,1);

  clear Store iVar y dn File Working


end; clear DayNumber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return