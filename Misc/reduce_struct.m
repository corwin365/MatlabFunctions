function Struct = reduce_struct(Struct,SubSetIndices,VarsToExclude,Dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%For a series of identical fields in a struct(), select a given set of 
%indices from every field
%
%For backwards-comptability with older versions, VarsToExclude needs to be
%ordered before Dim. If we want to operate on ALL fields, just set this to [].
%
%inputs:
%  Struct        - the struct to operate on
%  SubSetIndices - the list of indices to select from each field
%  VarsToExclude - fields to ignore when subsettings
%  Dim           - the dimension to operate on for each field
%
%outputs:
%  Struct        - the  reduced structure

%
%Corwin Wright, c.wright@bath.ac.uk, 2020/JUN/01
%updated 2023/08/12 to let the user choose a dimension to operate along
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('VarsToExclude','var'); VarsToExclude = {' '}; end  %assume applies to all vars
if ~exist(          'Dim','var'); Dim = 1;               end  %assume operation is along first dimension


Fields = fieldnames(Struct);
for iField=1:1:numel(Fields);
  
  %skip specified variables
  if any(strcmp(Fields{iField},VarsToExclude)); continue; end

  %reduce desired variables
  F = Struct.(Fields{iField});
  F = index_dim(F,SubSetIndices,Dim);
  Struct.(Fields{iField}) = F;
  
  %done!
  
end

return
