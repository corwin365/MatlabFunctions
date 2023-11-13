function [Data,FileList] = module_load_AIRS(Settings,InstInfo,Vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%module to prepare data from AIRS for get_limbsounders()
%
%this is NOT RECOMMENDED for use as AIRS is a nadir sounder, but it will work if absolutely needed
%this setting will load it as if it was a limb sounder
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

  for iGranule=InstInfo.Granules;

    %work out year and day number and hence filepath
    [y,~,~] = datevec(DayNumber); dn = date2doy(DayNumber);
    File = wildcardsearch([InstInfo.Path,'/',sprintf('%04d',y),'/'],['_',sprintf('%04d',y),'_',sprintf('%03d',dn),'_',sprintf('%03d',iGranule)]);

    if numel(File) == 0; clear y dn File; continue; end

    %load granule
    Working = rCDF(File{1});

    %store file information
    FileCount = FileCount+1;
    f = strsplit(File{1},filesep); FileList{end+1} = f{end}; clear f;

    %thin the data out
    Thinned = struct();
    Thinned.T    = Working.ret_temp(:,1:InstInfo.ThinFactor(1):end,:); Thinned.T    = Thinned.T( :,:,1:InstInfo.ThinFactor(2):end);
    Thinned.Lon  = Working.l1_lon(    1:InstInfo.ThinFactor(1):end,:); Thinned.Lon  = Thinned.Lon( :,1:InstInfo.ThinFactor(2):end);
    Thinned.Lat  = Working.l1_lat(    1:InstInfo.ThinFactor(1):end,:); Thinned.Lat  = Thinned.Lat( :,1:InstInfo.ThinFactor(2):end);
    Thinned.Time = Working.l1_time(   1:InstInfo.ThinFactor(1):end,:); Thinned.Time = Thinned.Time(:,1:InstInfo.ThinFactor(2):end);
    


    %load variables we need
    Store = struct();
    sz = size(Thinned.T);
    for iVar=1:1:numel(Vars)
      switch Vars{iVar}
        case 'Temp';       Store.Temp         = reshape(Thinned.T,sz(1),sz(2)*sz(3));
        case 'Lat';        Store.Lat          = reshape(Thinned.Lat,    sz(2)*sz(3),1);
        case 'Lon';        Store.Lon          = reshape(Thinned.Lon,    sz(2)*sz(3),1);
        case 'Alt';        Store.Alt          = Working.ret_z;
        case 'Pres';       Store.Pres         = h2p(Working.ret_z);
        case 'Time';       Store.Time         = datenum(2000,1,1,0,0,reshape(Thinned.Time,sz(2)*sz(3),1));
        case 'SourceProf'; Store.SourceProf   = (1:1:numel(Store.Lat))';
        case 'SourceFile'; Store.SourceFile   = ones(size(Store.SourceProf)).*FileCount;
        otherwise;
          try;   Store.(Vars{iVar}) = Working.(Vars{iVar});
          catch; disp(['Variable ',Vars{iVar},' not found, terminating']); return
          end
      end
    end

    %reshape vars
    Store.Lat  = repmat(Store.Lat, [1,sz(1)]);
    Store.Lon  = repmat(Store.Lon, [1,sz(1)]);
    Store.Time = repmat(Store.Time,[1,sz(1)]);
    Store.Temp = Store.Temp';
    Store.Pres = repmat(Store.Pres,[1,sz(2)*sz(3)])';
    Store.Alt  = repmat(Store.Alt,[1,sz(2)*sz(3)])';
    Store.SourceProf  = repmat(Store.SourceProf, [1,sz(1)]);
    Store.SourceFile  = repmat(Store.SourceFile, [1,sz(1)]);

    %store in main repository
    Data = cat_struct(Data,Store,1);
    
  end; clear iGranule

end; clear DayNumber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
