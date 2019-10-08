% parts = strsplit(name,sep);
%
% A simple strsplit implementation.
%
% Once a strsplit becomes widespread in octave and Matlab,
% this file can be deleted.
%
% This file is preferred over octave's strsplit, because it is compatible 
% with the Matlab interpreter.

function parts = strsplit(name,sep);

ind = find(name == sep);

parts = cell(length(ind)+1,1);

ind = [0 ind length(name)+1];

for i=1:length(parts)
    parts{i} = name(ind(i)+1:ind(i+1)-1);
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



