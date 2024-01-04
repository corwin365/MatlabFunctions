function Vector = force_col(Vector)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matlab function to force a 1D vector to be a column
%
%
%
%Corwin Wright, c.wright@bath.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check we have a 1D data series
sz = size(Vector);
if    numel(sz) > 2 || (sz(1) ~= 1 & sz(2) ~= 1);
  warning('make_col() is being applied to a variable with > 1 dimension, skipping operation')
  return
end

%ok, operate
if ~iscolumn(Vector); Vector = Vector'; end


end