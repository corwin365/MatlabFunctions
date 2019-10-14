function FilePath = era5_path(RootPath,DateNum)

%takes a Matlab datenumand generates a relative ERA5 path using my file
%storage conventions
%Corwin Wright, c.wright@bath.ac.uk, 14/OCT/2019


[y,~,~] = datevec(DateNum);
dn = date2doy(floor(DateNum));

FilePath = [RootPath,'/',sprintf('%04d',y),'/', ...
                     '/','era5_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.nc'];
return                   


