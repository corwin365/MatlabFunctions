function display(self)

sz = size(self);
tmp = sprintf('%dx',sz);
fprintf('Size: %s\n',tmp(1:end-1))

c = coord(self);

fprintf('Coordinates:\n')

for i = 1:length(c)
  tmp = sprintf('%dx',size(c(i).val));
  fprintf('  Name: "%s" standard name: "%s" size %s\n',...
          c(i).name,c(i).standard_name,tmp(1:end-1));
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

