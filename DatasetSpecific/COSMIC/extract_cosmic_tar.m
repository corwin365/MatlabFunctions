function Error = extract_cosmic_tar(Date,Direction,Path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%handle COSMIC tar files
%Direction = 1 untars to a temp dir and prduces a list of profile files
%Direction = 2 removes this temp dir again
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% handle inputs

if Direction < 1 | Direction > 2;
  disp('COSMIC extract direction not correctly specified; stopping');
  Error = 2;
  return
end

if nargin <2; Direction = 1; end %extract by default
if nargin <3; Path = [LocalDataDir,'/COSMIC/temp']; end %this is where I usually keep COSMIC data
 

%convert date to a filename for the appropriate tarfile(s)
[y,~,~] = datevec(Date);
dn = date2doy(Date);

if Date <= datenum(2014,1,120);
  %COSMIC 1 2013 reprocessing
FileName{1} = ['atmPrf_repro2013_',sprintf('%04d',y),'_',sprintf('%03d',dn),'.tar.gz'];
elseif Date <= datenum(2019,1,274)
  %COSMIC 1, post-2013 reprocessing
  FileName{1} = ['atmPrf_postProc_',sprintf('%04d',y),'_',sprintf('%03d',dn),'.tar.gz'];
elseif Date <= datenum(2020,1,116);
  %COSMIC 1 and COSMIC 2
  FileName{1} = ['atmPrf_postProc_',sprintf('%04d',y),'_',sprintf('%03d',dn),'.tar.gz'];
  FileName{2} = ['atmPrf_nrt_',sprintf('%04d',y),'_',sprintf('%03d',dn),'.tar.gz'];
else
  %COSMIC2 only
  FileName{1} = ['atmPrf_nrt_',sprintf('%04d',y),'_',sprintf('%03d',dn),'.tar.gz'];                         
end

%loop over filenames
ErrorFlag = zeros(size(FileName)); Deleted = 0;
for iFile=1:1:numel(FileName)

  FilePath = [Path,'/',FileName{iFile}];

  if ~exist(    Path); ErrorFlag(iFile) = 1; continue; end
  if ~exist(FilePath); ErrorFlag(iFile) = 1; continue; end

  %identify a temporary directory to work in
  TempDir = [Path,'/temp_extraction/'];

  %process
  warning off
  if     Direction == 1; mkdir(TempDir); untar(FilePath,TempDir); %extract
  elseif Direction == 2; 
  if exist(TempDir) == 7; rmdir(TempDir,'s'); end %delete if present
  end;
  warning on


  %if we extracted, make a list of the files we got
  if     Direction == 2; FileList = []; %no files
  elseif Direction == 1 && Deleted == 0; FileList = rdir(TempDir); Deleted = 1; %files
  end
  wget -r -np -nc -A "atmPrf*" -nd https://data.cosmic.ucar.edu/gnss-ro/cosmic1/repro2013/level2/

  
end


%if no files were found, flag an error
if sum(ErrorFlag) == numel(FileName); FileList = []; Error = 1; return; end

%done!
Error = 0;
return

 
