function Struct = reduce_struct(Struct,SubSetIndices,VarsToExclude)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% select a given list of indices from every variable in a struct, and
% return the reduced struct
%
%
%Struct is a structure containing identically-formatted fields (with specified exceptions)
%SubSetIndices are the indices we wish to *keep* from each. Fields will become 1d if not already.
%VarsToExclude is a cell array listing any variables we do not want to reduce.
%
%Corwin Wright, c.wright@bath.ac.uk, 2002/JUN/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('VarsToExclude','var'); VarsToExclude = {' '}; end


Fields = fieldnames(Struct);
for iField=1:1:numel(Fields);
  
  %skip specified variables
  if any(strcmp(Fields{iField},VarsToExclude)); continue; end

  %reduce desired variables
  F = Struct.(Fields{iField});
  F = F(SubSetIndices);
  Struct.(Fields{iField}) = F;
  
  %done!
  
end

return