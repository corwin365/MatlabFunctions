function Dimensions = cjw_nc_prepop_dim(Dimensions,VarName,Axis)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to generate prepopulated dimension fields for my netCDF wrappers with default settings
%purely exists to reduce typing
%
%Corwin Wright, c.wright@bath.ac.uk, 05/January/2019
%
%wrapper to several functions which are themselves wrappers to the built-in
%netCDF libraries! 
%
%inputs:   Dimensions: the Dimensions structure we will be adding to
%          VarName: the name of the variable
%          Axis: the axis values (can just be 1-N)
%outputs:  Dimensions: the resulting data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(fieldnames(Dimensions)) == 0
  %this is the first pass. 
  NDimensions = 0; %this wil trick the logic below.
else
  NDimensions = size(Dimensions,2);
end


Dimensions(NDimensions+1).Name     = VarName;
Dimensions(NDimensions+1).FullName = VarName; %better than nothing!
Dimensions(NDimensions+1).Units    = ''; %can be empty
Dimensions(NDimensions+1).Fill     = -9999; %default fill value
Dimensions(NDimensions+1).Axis     = Axis;
Dimensions(NDimensions+1).Type     = 'double'; %by default

return

end
