function [FileList,Error] = extract_cosmic_tar(Date,Direction,Path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%handle COSMIC tar files
%Direction = 1 untars to a temp dir and prduces a list of profile files
%Direction = 2 removes this temp dir again
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% handle inputs

if Direction < 1 | Direction > 2;
  disp('COSMIC extract direction not specified; stopping');
  FileList = [];
  return
end

if nargin <2; Direction = 1; end %extract by default
if nargin <3; Path = [LocalDataDir,'/COSMIC']; end %this is where I usually keep COSMIC data
 

%% process

%convert date to a filename for the appropriate tarfile
[y,~,~] = datevec(Date);
dn = datevec2doy(datevec(Date));
FileName = ['cosmic2013_atmPrf_',sprintf('%04d',y),'.',sprintf('%03d',dn),'.tar'];
FilePath = [Path,'/',FileName];


if ~exist(FilePath); FileList = []; Error = 1; return; end

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
elseif Direction == 1; FileList = rdir(TempDir); %files
end

%done!
Error = 0;
return


