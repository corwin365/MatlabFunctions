function Success = cjw_nc_makedims(NcFile,Dimensions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to set dimensions in a new netCDF file
%Corwin Wright, c.wright@bath.ac.uk, 05/January/2019
%
%inputs: NCFile:     name of file to add dimensions to
%        Dimensions: cell structure - example given below header of format
%                    ALL VALUES MUST BE SET, EVEN IF JUST AS AN EMPTY STRING
%outputs: Success: 1: successful
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for two axes, you would input them like this:

% % Dimensions(1).Name     = 'Profiles';
% % Dimensions(1).FullName = 'Profile Number in File';
% % Dimensions(1).Units    = '';
% % Dimensions(1).Fill     = -9999;
% % Dimensions(1).Axis     = 1:1:13;
% % Dimensions(1).Type     = 'single';
% % 
% % Dimensions(2).Name     = 'Altitude';
% % Dimensions(2).FullName = 'Altitude above surface';
% % Dimensions(2).Units    = 'km';
% % Dimensions(2).Fill     = -9999;
% % Dimensions(2).Axis     = 1:1:27;
% % Dimensions(2).Type     = 'double';


%open the netCDF file
FileId = netcdf.open(NcFile,'WRITE');


%create the axis data
for iDim=1:1:size(Dimensions,2)
  
  %first, make sure there are no invalid characters in the name, by just dropping illegal characters
  %if this causes duplication, tough.
  VarName = matlab.lang.makeValidName(Dimensions(iDim).Name);
  
  
  %create a variable representing it
  DimIDs.(VarName) = netcdf.defDim(FileId,                ... %file id
                                  Dimensions(iDim).Name, ... %use the true name here, not the matlab-safe one
                                  numel(Dimensions(iDim).Axis));
                                
% % %   %create the variable
% % %   VarId = netcdf.defVar(FileId,                ...
% % %                         Dimensions(iDim).Name, ...
% % %                         Dimensions(iDim).Type, ...
% % %                         DimIDs.(VarName));
% % %                                 
% % %   %fill in the metadata                    
% % %   netcdf.putAtt(FileId,VarId,'standard_name',Dimensions(iDim).Name);
% % %   netcdf.putAtt(FileId,VarId,'long_name',    Dimensions(iDim).FullName);
% % %   netcdf.putAtt(FileId,VarId,'units',        Dimensions(iDim).Units);
% % %   
% % %   %define the fill value
% % %   netcdf.defVarFill(FileId,VarId,false,Dimensions(iDim).Fill);
% % %   
% % %   %and fill in the axis values
% % %   netcdf.putVar(FileId,VarId,Dimensions(iDim).Axis);
  
end

%and close the netCDF file
netcdf.close(FileId);



%done
Success = 1;

end

