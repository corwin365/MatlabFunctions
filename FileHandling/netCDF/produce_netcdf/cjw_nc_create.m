function Success = cjw_nc_create(NCFile,MetaData,Clobber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to create a netCDF file and store metadata
%Corwin Wright, c.wright@bath.ac.uk, 05/January/2019
%
%inputs: NCFile:   name of file to create
%        MetaData: struct containing fields of strings, which will be transposed to netCDF global attributes
%        Clobber:  1 to overwrite any existing file, any other value to leave intact (default 0)
%
%outputs:  Success: 1: successful
%                   2: problem creating netcdf file
%                   3: problem setting global attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%
%create file
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3; Clobber = 0; end
if Clobber ~= 1; Clobber = 0; end

%always netCDF4
Mode = netcdf.getConstant('NETCDF4');

%clobber?
if Clobber == 1; Mode = bitor(Mode,netcdf.getConstant('CLOBBER'));
else             Mode = bitor(Mode,netcdf.getConstant('NOCLOBBER'));
end

%go!
try
  FileId = netcdf.create(NCFile,Mode);
catch
  disp('Error: unable to create netCDF file')
  Success = 2;
  return
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%set global attributes
%%%%%%%%%%%%%%%%%%%%%%%%%

%variable name for accessing global attributes
Global = netcdf.getConstant('NC_GLOBAL');

%time file was created
netcdf.putAtt(FileId,Global,'CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'))

%set other attributes (can include an override of the above if the var is set)
Attributes = fieldnames(MetaData);
for iAttribute = 1:1:numel(Attributes)
  
  %check the attributes is a character array
  if ~ischar(MetaData.(Attributes{iAttribute}))
    %if not, print out and continue
    disp(['Global variable ',MetaData.(Attributes{iAttribute}),' not set - not text'])
  else
    %set variable
    netcdf.putAtt(FileId,Global,Attributes{iAttribute},MetaData.(Attributes{iAttribute}));
  end
end

%close file
netcdf.close(FileId)
Success = 1;


end

