function Output = get_HIRDLS(FilePath,Variable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract a single variable from a HIRDLS HDF5 data file
%returns NaN if variable does not exist
%
%adapted from a previous version developed by N Hindley, U Bath, May 2014
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%11/MAY/2014
%
%inputs
%---------
%
%FilePath: path to file
%Variable: variable wanted from file
%
%outputs
%---------
%
%Output: variable contents (or NaN if failed)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------
%list possible variables in an HDF5-EOS-HIRDLS file, by type
%------------------------------------------------------------------------------------

GeolocationFields = {'Altitude';'Latitude';'LocalSolarTime';'Longitude'; ...
                     'OrbitAscendingFlag';'OrbitNumber';'Pressure';      ...
                     'ProfileID';'ScanAzimuthAtNominalAltitude';         ...
                     'ScanElevationAtNominalAltitude';'ScanTable';       ...
                     'ScanUpFlag';'ScienceScanMode';'SecondsInDay';      ...
                     'SolarZenithAngle';'SpacecraftAltitude';            ...
                     'SpacecraftLatitude';'SpacecraftLongitude';         ...
                     'TangentHeightAtNominalAltitude';'Time';            ...
                     'ViewDirectionAtNominalAltitude'};

DataFields = {'10.8MicronCloudAerosolFlag';'10.8MicronExtinction';             ...
              '10.8MicronExtinctionNormChiSq';'10.8MicronExtinctionPrecision'; ...
              '10.8MicronExtinctionQuality';'12.1MicronCloudAerosolFlag';      ...
              '12.1MicronExtinction';'12.1MicronExtinctionNormChiSq';          ...
              '12.1MicronExtinctionPrecision';'12.1MicronExtinctionQuality';   ...
              '17.4MicronCloudAerosolFlag';'17.4MicronExtinction';             ...
              '17.4MicronExtinctionNormChiSq';'17.4MicronExtinctionPrecision'; ...
              '17.4MicronExtinctionQuality';'7.1MicronCloudAerosolFlag';       ...
              '7.1MicronExtinction';'7.1MicronExtinctionNormChiSq';            ...
              '7.1MicronExtinctionPrecision';'7.1MicronExtinctionQuality';     ...
              '8.3MicronCloudAerosolFlag';'8.3MicronExtinction';               ...
              '8.3MicronExtinctionNormChiSq';'8.3MicronExtinctionPrecision';   ...
              '8.3MicronExtinctionQuality';'CFC11';'CFC11NormChiSq';           ...
              'CFC11Precision';'CFC11Quality';'CFC12';'CFC12NormChiSq';        ...
              'CFC12Precision';'CFC12Quality';'CH4';'CH4NormChiSq';            ...
              'CH4Precision';'CH4Quality';'ClONO2';'ClONO2NormChiSq';          ...
              'ClONO2Precision';'ClONO2Quality';'CloudTopPressure';'GPH';      ...
              'GPHPrecision';'H2O';'H2ONormChiSq';'H2OPrecision';'H2OQuality'; ...
              'HNO3';'HNO3NormChiSq';'HNO3Precision';'HNO3Quality';'N2O';      ...
              'N2O5';'N2O5NormChiSq';'N2O5Precision';'N2O5Quality';            ...
              'N2ONormChiSq';'N2OPrecision';'N2OQuality';'NO2';'NO2NormChiSq'; ...
              'NO2Precision';'NO2Quality';'O3';'O3NormChiSq';'O3Precision';    ...
              'O3Quality';'RawGPH';'RawGPHPrecision';'Temperature';            ...
              'TemperatureNormChiSq';'TemperaturePrecision';'TemperatureQuality'};

%------------------------------------------------------------------------------------
%identify type of desired variable
%------------------------------------------------------------------------------------
   
IsGeo  = numel(find(ismember(GeolocationFields,Variable)));
IsData = numel(find(ismember(DataFields,Variable       )));

%------------------------------------------------------------------------------------
%load data
%------------------------------------------------------------------------------------

%open hdf5 file
FileID = H5F.open (FilePath, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');


if IsGeo == 1; 
  FieldPath = strcat('HDFEOS/SWATHS/HIRDLS/Geolocation Fields/',Variable);
elseif IsData ==1; 
  FieldPath = strcat('HDFEOS/SWATHS/HIRDLS/Data Fields/',Variable);
else
  disp(['Error: ''',Variable,''' is not present in this file format']); 
  Output = NaN;
  return
end;

%extract data and close file
data_id = H5D.open (FileID,FieldPath);
Output = H5D.read(data_id,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT');
H5F.close(FileID)


return,Output
end








