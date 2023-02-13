function Data = cjw_nc_prepop_data(Data,VarName,Class)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to generate prepopulated data fields for my netCDF wrappers with default settings
%purely exists to reduce typing
%
%Corwin Wright, c.wright@bath.ac.uk, 05/January/2019
%
%wrapper to several functions which are themselves wrappers to the built-in
%netCDF libraries! 
%
%inputs:   Data: the data structure we will be adding to
%          VarName: the name of the variable
%
%OPTIONAL: Class: class of the variable
%
%outputs:  Data: the resulting data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Data.(VarName).Dims = []; %this must always be specified outside
Data.(VarName).VarName  = VarName;
Data.(VarName).FullName = VarName; %better than nothing
Data.(VarName).Units    = ' '; %can be blank
Data.(VarName).Fill     = -9999; %default to filling with -9999
Data.(VarName).Data     = []; %this must always be specified
Data.(VarName).Type    = 'double'; %default to double type

if nargin ==3;
  switch Class
    case 'int8';   Data.(VarName).Type = 'NC_BYTE';
    case 'char';   Data.(VarName).Type = 'NC_CHAR'; 
    case 'string'; Data.(VarName).Type = 'NC_STRING';  
    case 'int16';  Data.(VarName).Type = 'NC_SHORT';  
    case 'int32';  Data.(VarName).Type = 'NC_INT';       
    case 'single'; Data.(VarName).Type = 'NC_FLOAT';
    case 'double'; Data.(VarName).Type = 'NC_DOUBLE';      
    otherwise      Data.(VarName).Type = 'NC_DOUBLE';
  end
end
return

end
