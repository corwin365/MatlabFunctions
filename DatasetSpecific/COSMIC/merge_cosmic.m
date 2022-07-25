clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce merged daily COSMIC data files
%Corwin Wright, c.wright@bath.ac.uk
%2022/05/14
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%wget -r -np -nc -A "atmPrf*" -nd https://data.cosmic.ucar.edu/gnss-ro/cosmic1/repro2013/level2/


Settings.TimeScale = datenum(2006,1,1):1:datenum(2022,5,10);
Settings.InDir     = '/media/DB/COSMIC/temp';%[LocalDataDir,'/COSMIC/temp'];
Settings.OutDir    = '/media/DB/COSMIC/daily_atmPrf';%[LocalDataDir,'/COSMIC/temp'];%Settings.InDir;



for iDay=1:1:numel(Settings.TimeScale)

  %locate file for this day
  disp(['Processing ',datestr(Settings.TimeScale(iDay))]);

  [y,~,~] = datevec(Settings.TimeScale(iDay));
  dn =date2doy(Settings.TimeScale(iDay));
  OutFile = [Settings.OutDir,'/',sprintf('%04d',y),'/cosmic_',sprintf('%04d',y),'_',sprintf('%03d',dn),'.mat'];
  
  if exist(OutFile); disp('Already done'); continue; end



%  try
  %untar file to get individual profiles
  Error = extract_cosmic_tar(Settings.TimeScale(iDay),1,Settings.InDir);
  if Error ~= 0; disp('Input data not found, skipping day'); continue; end

  %get list of profiles
  FileList = wildcardsearch([Settings.InDir,'/temp_extraction'],'*_nc');
  disp('Extracted files')

  %loop over list of profiles and merge into single struct

  MasterStore = struct();

  textprogressbar('Processing individual profiles ')
  for iFile=1:1:numel(FileList)

%      try
      %load file
      Data = rCDF(FileList{iFile});


      %compute time of profile
      Globals = Data.MetaData.Attributes.Global;
      Store = NaN(6,1);
      for iGlobal = 1:1:numel(Globals)
        ID = Globals(iGlobal).Name;
        if strcmp(ID,  'year'); Store(1) = Globals(iGlobal).Value; end
        if strcmp(ID, 'month'); Store(2) = Globals(iGlobal).Value; end
        if strcmp(ID,   'day'); Store(3) = Globals(iGlobal).Value; end
        if strcmp(ID,  'hour'); Store(4) = Globals(iGlobal).Value; end
        if strcmp(ID,'minute'); Store(5) = Globals(iGlobal).Value; end
        if strcmp(ID,'second'); Store(6) = Globals(iGlobal).Value; end
      end
      Time = datenum(Store(1),Store(2),Store(3),Store(4),Store(5),Store(6));
      clear Store Globals ID iGlobal

      %store data
      Fields = {'Lat','Lon','MSL_alt','Ref','Azim','Pres','Bend_ang','Impact_height','Temp'};
      for iF=1:1:numel(Fields);
        MasterStore(iFile).(Fields{iF}) = single(Data.(Fields{iF}));
      end;
      clear iF Fields
      MasterStore(iFile).Time = Time;
      clear Time
%      catch; end


    textprogressbar(iFile./numel(FileList).*100)
  end; clear iFile FileList
  textprogressbar(100);textprogressbar('!')


  %save to new file

  Data = MasterStore;
  save(OutFile,'Data')
  clear MasterStore y dn OutFile
%  catch; end

  %clear up tarball mess
  Error = extract_cosmic_tar(Settings.TimeScale(iDay),2,Settings.InDir);
  clear FileList Error
  disp('Saved and extracted profiles deleted')


end
