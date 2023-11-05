function [Data,FileList] = module_load_ACE_GNSS(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from ACE-FTS and GNSS for get_limbsounders()
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
  switch Settings.Instrument
    case 'GNSS'; File = [InstInfo.Path,'/',sprintf('%04d',y),filesep,'merged_ro_',sprintf('%04d',y),'_',sprintf('%03d',dn),'.mat'];
    case  'ACE'; File = [InstInfo.Path,'/',sprintf('%04d',y),filesep,'ace_v52_',  sprintf('%04d',y),'d',sprintf('%03d',dn),'.mat'];
  end
  if ~exist(File,'file'); clear y dn File; continue; end

  %load the data
  InstData = load(File);

  %store file information
  FileCount = FileCount+1;
  f = strsplit(File,filesep); FileList{end+1} = f{end}; clear f;
  clear y dn File

  %get the variables we want
  Store = struct();
  for iVar=1:1:numel(Vars)

    if ~strcmp(Vars{iVar},'Alt')  ...
        & ~strcmp(Vars{iVar},'Time') ...
        & ~strcmp(Vars{iVar},'SourceProf') ...
        & ~strcmp(Vars{iVar},'SourceFile') ...
        & ~strcmp(fieldnames(InstData),Vars{iVar});
      disp(['Variable ',Vars{iVar},' not found, terminating'])
      return
    end

    switch Vars{iVar}
      case 'Time';
        switch Settings.Instrument
          case 'GNSS'; Store.Time = repmat(InstData.MetaData.time,[1,size(InstData.Lat,2)]);
          case 'ACE';  Store.Time = InstData.Time;
        end
      case 'Alt';
        switch Settings.Instrument
          case 'GNSS'; Store.Alt  = repmat(InstData.MSL_alt,[size(InstData.Lat,1),1]);
          case 'ACE';  Store.Alt = InstData.Alt;
        end
      case 'SourceProf';Store.SourceProf = single(repmat(1:1:size(Store.Lat,1),size(Store.Lat,2),1)');
      case 'SourceFile';Store.SourceFile = ones(size(Store.SourceProf)).*FileCount;
      otherwise;   Store.(Vars{iVar}) = InstData.(Vars{iVar});
    end

  end; clear iVar

  %store in main repository
  Data = cat_struct(Data,Store,1);

  clear Store iVar

end; clear DayNumber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return