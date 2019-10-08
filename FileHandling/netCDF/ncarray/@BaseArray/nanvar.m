% Compute the variance (ignoring NaNs).
% V = var (X, OPT, DIM)
% Compute the variance along dimension DIM.
% If OPT is equal to 1, then the variance is bias-corrected.

function s = nanvar(self,opt,varargin)

if nargin == 1
  opt = 0;
elseif isempty(opt)
  opt = 0;
end

nm = nanmean(self,varargin{:});

[s,n] = reduce(self,@funred, ...
               @(x) funelem(x,nm),varargin{:});
if isempty(s)
  s = 0;
else
  total = s{1};
  count = s{2};

  if opt == 0
    s = total./(count-1);
  else
    s = total./count;
  end
end
end

function x = funelem(x,nm)
  % make sure x is not an ncArray
  x = full(x);
  mask = isnan(x) || isnan(nm);
  diff = zeros(size(x));
  diff(mask) = x(mask) - nm(mask);
  x = {diff.^2, ~mask};  
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

