function FilePath = era5_path(DateNum,RootPath,Check)

%takes a Matlab datenumand generates a relative ERA5 path using my file
%storage conventions
%Corwin Wright, c.wright@bath.ac.uk, 14/OCT/2019

%modified 2021/01/11 to add additional option to check file exists and, 
%if not, also check for ERA5T data and return that instead if found


if nargin == 1; RootPath = LocalDataDir; end
if nargin < 3; Check = 0; end

%identify year and doy
[y,~,~] = datevec(DateNum);
dn = date2doy(floor(DateNum));

%generate filepath
FilePath = [RootPath,'/ERA5/',sprintf('%04d',y),'/', ...
                     '/','era5_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.nc'];
                   
if Check == 1
  if ~exist(FilePath,'file');
    FilePath = [RootPath,'/ERA5/',sprintf('%04d',y),'/', ...
                         '/','era5t_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.nc'];
    if ~exist(FilePath,'file'); Filepath = '';
    end
  end
end
  
                   
                   
return                  


