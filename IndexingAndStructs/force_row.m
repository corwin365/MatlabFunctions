function Vector = force_row(Vector,Silent)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matlab function to force a 1D vector to be a row
%
%will skip arrays with more than 1 dimension, giving a warning unless "Silent" is set to 1
%
%Corwin Wright, c.wright@bath.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check Silent var exists
if nargin < 2; Silent = 0; end

%check we have a 1D data series
sz = size(Vector);
if    numel(sz) > 2 || (sz(1) ~= 1 & sz(2) ~= 1);
  if Silent ~=1; warning('force_row() is being applied to a variable with > 1 dimension, skipping operation'); end
  return
end

%ok, operate
if ~isrow(Vector); Vector = Vector'; end


end