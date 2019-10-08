function [Row,Success] = cjw_nc_prepop_data(Data,VarName)

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
%outputs:  Data: the resulting data structure
%          Success: 1 if done
%                   0 if error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%
%incorrect request?
%%%%%%%%%%%%%%%%%%%%%%%%
if DataDim < 1 | DataOrDim > 2; Success = 0; Row = []; return; end


%%%%%%%%%%%%%%%%%%%%%%%%
%data field
%%%%%%%%%%%%%%%%%%%%%%%%

if DataDim == 1
  
  Field.Dims     = []; %this must always be specified
  Field.VarName  = 'VARNAME'; %should be overwritten, this is a filler so it doesn't crash
  Field.FullName = ''; %can be blank
  Field.Units    = ' '; %can be blank
  Field.Fill     = -9999; %default to filling with -9999
  Field.Data     = []; %this must always be specified
  Field.Type    = 'double'; %default to double type

end

%%%%%%%%%%%%%%%%%%%%%%%%
%dimension field
%%%%%%%%%%%%%%%%%%%%%%%%

if DataDim == 1
  Dimensions(1).Name     = 'Profiles';
  Dimensions(1).FullName = 'Profile Number in File';
  Dimensions(1).Units    = '';
  Dimensions(1).Fill     = -9999;
  Dimensions(1).Axis     = 1:1:size(Store.Profiles,1);
  Dimensions(1).Type     = 'double';
end

%%%%%%%%%%%%%%%%%%%%%%%%
%define dimensions
%%%%%%%%%%%%%%%%%%%%%%%%



Dimensions(2).Name     = 'Altitude';
Dimensions(2).FullName = 'Altitude above surface';
Dimensions(2).Units    = 'km';
Dimensions(2).Fill     = -9999;
Dimensions(2).Axis     = 1:1:size(Store.Profiles,2);
Dimensions(2).Type     = 'double';

cjw_nc_makedims(Settings.NcFile,Dimensions);

