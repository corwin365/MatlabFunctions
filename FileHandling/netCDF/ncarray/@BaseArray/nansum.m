% Compute the sum (ignoring NaNs).
% [s,count] = sum (X, DIM)
% Compute the sum along dimension DIM.
% The variable s is the sum and count the number of values summed.

function [total,count] = nansum(self,varargin)

% s will be a cell element containing 
% {the sum, the count of elements different from NaN}
% this is necessary to avoid 2 calls to reduce 

[s,n] = reduce(self,...
               @(x,y) funred(x,y),...
               @(x) funelem(x),...
               varargin{:});

if isempty(s)
  s = 0;
else
  total = s{1};
  count = s{2};
end

end

function x = funelem(x)
  % make sure x is not an ncArray
  x = full(x);
  m = isnan(x);
  x(m) = 0;
  x = {x,~m};
end

function s = funred(x,y)
  s = {x{1} + y{1}, x{2} + y{2}} ;
end



% Copyright (C) 2013 Alexander Barth <barth.alexander@gmail.com>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.

