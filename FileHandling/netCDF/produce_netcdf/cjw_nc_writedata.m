function Success = cjw_nc_writedata(NcFile,Data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to write data to a new netCDF file
%Corwin Wright, c.wright@bath.ac.uk, 05/January/2019
%
%inputs: NCFile: name of file to add dimensions to
%        Data:   cell structure - example given below header of format
%                ALL VALUES MUST BE SET, EVEN IF JUST AS AN EMPTY STRING
%                DIMENSION NUMBERS MUST CORRESPOND WITH THOSE IN THE FILE
%outputs: DimIDs: netCDF file identifiers for each dimension
%         Success: 1: successful
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %for a 1d var in dimension 1, you would input like this:
% % 
% % Data.Latitude.Dims     = [1]; 
% % Data.Latitude.VarName  = 'Latitude';
% % Data.Latitude.FullName = 'Latitude';
% % Data.Latitude.Units    = 'degrees N';
% % Data.Latitude.Fill     = -9999;
% % Data.Latitude.Data     = Store.Geo(:,2);
% % Data.Latitude.Type    = 'double';
% % 
% % %for a 1d var in dimensions 1 and 2, you would input like this:
% % Data.Temperature.Dims     = [1,2]; 
% % Data.Temperature.VarName  = 'Temperature';
% % Data.Temperature.FullName = 'Temperature';
% % Data.Temperature.Units    = 'Kelvin';
% % Data.Temperature.Fill     = -9999;
% % Data.Temperature.Data     = Store.Profiles;
% % Data.Temperature.Type     = 'double';


%open the netCDF file
FileId = netcdf.open(NcFile,'WRITE');


%loop over variables
Variables = fieldnames(Data);
for iVar=1:1:numel(Variables)
  
  
  %first, make sure there are no invalid characters in the name, by just dropping illegal characters
  %if this causes duplication, tough.
  VarName = matlab.lang.makeValidName(Data.(Variables{iVar}).VarName);
  
  
  %create the variable
  VarId = netcdf.defVar(FileId,                         ...
                        Data.(Variables{iVar}).VarName, ...
                        Data.(Variables{iVar}).Type,    ...
                        Data.(Variables{iVar}).Dims-1); %the -1 is important, as netCDF is zero-indexed and Matlab is 1-indexed
                      
  %fill in the attributes
  netcdf.putAtt(FileId,VarId,'name',     Data.(Variables{iVar}).VarName);
  netcdf.putAtt(FileId,VarId,'long_name',Data.(Variables{iVar}).FullName);  
  netcdf.putAtt(FileId,VarId,'units',    Data.(Variables{iVar}).Units);      
  
  %specify fill value
  netcdf.defVarFill(FileId,VarId,false,Data.(Variables{iVar}).Fill);
  
  %now, replace any bad values in the data with the fill value...
  v = Data.(Variables{iVar}).Data; 
  v(isnan(v)) = Data.(Variables{iVar}).Fill;
                      

  %and write it
  netcdf.putVar(FileId,VarId,v); 
  
  %done!
  
end


%and close the netCDF file
netcdf.close(FileId);


%done
Success = 1;



end

