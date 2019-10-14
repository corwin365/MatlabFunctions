function FilePath = era5_path(DateNum,RootPath)

%takes a Matlab datenumand generates a relative ERA5 path using my file
%storage conventions
%Corwin Wright, c.wright@bath.ac.uk, 14/OCT/2019


if nargin == 1; RootPath = LocalDataDir; end

%identify year and doy
[y,~,~] = datevec(DateNum);
dn = date2doy(floor(DateNum));

%generate filepath
FilePath = [RootPath,'/ERA5/',sprintf('%04d',y),'/', ...
                     '/','era5_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.nc'];
                   
                   
return                  


