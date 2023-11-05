function [Data,FileList] = module_load_SABER(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from SABER for get_limbsounders()
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

%monthly files, ugh
OldData.Data = struct();
OldData.Name = '';

for DayNumber=floor(min(Settings.TimeRange))-1:1:floor(max(Settings.TimeRange));  %we need to use the day before the start as well as the underlying data format is orbit-based, not day-based

  %identify file and load, but only if we didn't on a previous loop
  [y,~,~] = datevec(DayNumber);
  m = month(datetime(DayNumber,'ConvertFrom','datenum'),'name'); m = m{1};
  File = wildcardsearch(InstInfo.Path,[m,num2str(y)]);

  if numel(File) == 0; clear y m File; continue; end

  if ~strcmp(OldData.Name,File{1});

    %load file
    OldData.Name = File{1};
    OldData.Data = rCDF(File{1});

    %store file information
    FileCount = FileCount+1;
    f = strsplit(File{1},filesep); FileList{end+1} = f{end}; clear f;

    %put dates in a useful unit. We'll handle time within day later, as this can be quite computationally expensive.
    for iEl=1:1:numel(OldData.Data.date)
      id = num2str(OldData.Data.date(iEl));
      OldData.Data.Date(iEl) = datenum(str2double(id(1:4)),1,str2double(id(5:7)),0,0,0);
    end

    %create source profile indices
    OldData.Data.SourceProf = repmat(1:1:numel(OldData.Data.orbit),size(OldData.Data.time,1),1);
    OldData.Data.SourceFile = ones(size(OldData.Data.SourceProf)).*FileCount;

  end


  clear File y m

  %select the day we actually want
  Working = OldData.Data;
  OnThisDay = find(floor(Working.Date) == DayNumber);
  Working = reduce_struct(Working,OnThisDay,{'orbit','date','MetaData'},2);
  Working.orbit = Working.orbit(OnThisDay); Working.date = Working.date(OnThisDay);
  clear OnThisDay

  %now, extract the values we want
  Store.Lat = Working.tplatitude';
  Store.Lon = Working.tplongitude';
  Store.Pres = Working.pressure';
  Store.Alt  = Working.tpaltitude';
  Store.SourceProf  = Working.SourceProf';
  Store.SourceFile  = Working.SourceFile';
  Store.Temp = Working.ktemp';

  %remove bad heights, as they break the interpolation below
  Store.Alt(Store.Alt > 1000) = NaN;

  %time is a bit more awkward. This is horrible, but it is what it is.
  Working.date = repmat(Working.date,[1,size(Working.time,1)])';
  year = floor(Working.date./1000);
  dn   = Working.date - year.*1000;
  s    = Working.time./1000; s(s < 0) = NaN;
  Store.Time = datenum(year,1,dn,1,1,s)';

  %store in main repository, and tidy
  Data = cat_struct(Data,Store,1);

end; clear  dn id iEl iVar OldData s Store year working
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return