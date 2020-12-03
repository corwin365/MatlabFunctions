
function Result = cjw_writenetCDF(FileName,Data,Clobber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%write input structure to a specified netCDF file
%
%
%Warning: not heavily tested, likely to be buggy if you try anything odd!
%in particular, error handling is not great, and not all types are currently included
%
%to do:
% 1. error handling
% 2. better cell array handling (currently only handles 1d, by converting to 2d string)
% 3. allow for structures within structures
%
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%06/JAN/2014
%
%modified:
%
%CJW 13/JAN/2014 convert cells to string arrays (crudely)
%
%inputs
%---------
%
%FileName - (string   )  filename to write to
%Data     - (structure)  data to write to file
%Clobber  - (1 or 0   )  overwrite file if already exists (1 yes 0 no, default no)
%
%outputs
%---------
%
%Return - 1 if successful
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%handle imputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 3); Clobber = 0; end; %don't overwrite files by default



% create the netCDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Clobber == 0;
  
  %check if file exists
  if exist(FileName,'file')
    disp(' ');
    disp('*****NOCLOBBER selected and file exists. Writing aborted.*****');
    Result = 0;
    return;
  end

  %if it doesn't, create the file. Carefully...
  NetCdfId = netcdf.create(FileName,'NOCLOBBER');

else             NetCdfId = netcdf.create(FileName,'CLOBBER');
end;

%work out what variables we want to put in it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StructureContents = fieldnames(Data);
NVariables = numel(StructureContents);

DimCount = 1;
for iVar=1:1:NVariables;

  %check variable type
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %check variable type and rename to whatever netCDF has as the equivalent
  VarClass = class(eval(char(strcat('Data.',StructureContents(iVar)))));
  
  if strcmp(VarClass,'single'); VarClass = 'float'; end; %single --> float
  %%%                                                    %double --> double
  if strcmp(VarClass,'int8');   VarClass = 'byte';  end; %int8 --> byte  
  if strcmp(VarClass,'int32');  VarClass = 'int';   end; %int32 --> int  

  %if cell array, convert to string array and save as that
  %format  might look a bit weird if loaded into another language!
  if strcmp(VarClass,'cell');
    eval(char(strcat('Data.',StructureContents(iVar),'=char(','Data.',StructureContents(iVar),');')));
    VarClass = 'char';
  end; clear NAcross;
  
  
  %read structure variable into a working variable
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  eval(char(strcat('WorkingVariable = Data.',StructureContents(iVar),';')));
  
  %start producing variable definition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Vardef = '[';
  
  %define dimensions, and continue producing variable definition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for iDim=1:1:ndims(WorkingVariable);
    
    %find dimension size
    DimSize = size(WorkingVariable,iDim);
    
    %define dimension
    Dim  = netcdf.defDim(NetCdfId,num2str(DimCount), DimSize);
  
    %add to definition of variable
    Vardef = char(strcat([Vardef,' ',num2str(Dim),' ']));

    %next dimension!
    DimCount = DimCount + 1;
  end; clear iDim
  Vardef = char(strcat(Vardef,' ]'));
  
  
  %finish and fire off variable definition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  Vardef = char(strcat('Var',                          ...
                       num2str(iVar),                  ...
                       ' = netcdf.defVar(NetCdfId,''', ...
                       StructureContents(iVar),        ...
                       ''', ''',                       ...
                       VarClass,                       ...
                       ''', ',Vardef,                  ...
                       ');'));
  eval(Vardef);
end; clear iVar Vardef DimCount WorkingVariable;


%done defining the NetCdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 netcdf.endDef(NetCdfId);
 
 
%store data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for iVar=1:1:NVariables;
    
  %read structure variable into a working variable
  eval(char(strcat('netcdf.putVar(NetCdfId,Var',num2str(iVar),',Data.',StructureContents(iVar),');')));

end ; clear iVar;

%close netCDF file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
netcdf.close(NetCdfId);
  
Result = 1; %finished

