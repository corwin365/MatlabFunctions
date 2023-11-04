function [Data] = read_imerg(File)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to read IMERG native HDF5 format and return the key data as
%a Matlab struct
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/04/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%list of variables in an IMERG HDF file
Vars = {'lat','lon','precipitation','probabilityLiquidPrecipitation','randomError','time',};

%bad data value
BadFlag = -9999;

%load vars
Data = struct();
for iVar=1:1:numel(Vars)

  %get variable
  Var = h5read(File,['/Grid/',Vars{iVar}]);

  %replace bad data with NaNs
  Bad = find(Var <= BadFlag);
  Var(Bad) = NaN;

  %store
  Data.(Vars{iVar}) = Var;

end

%convert time to Matlab
Data.time = datenum(1980,1,6,0,0,double(Data.time));

%produce some metadata to help new reader
MetaData = struct();
MetaData.lat                            = "degrees north at box centre";
MetaData.lon                            = "degrees east at box centre";
MetaData.precipitation                  = "Complete merged microwave-infrared (gauge-adjusted), mm/hr";
MetaData.probabilityLiquidPrecipitation = "Probability of liquid precipitation estimated with a diagnostic parameterization using ancillary data, percent";
MetaData.randomError                    = "Root-mean-square error estimate for complete merged microwave-infrared (gauge-adjusted) precipitation, mm/hr";
MetaData.time                           = "Matlab units";


%done!
return
end