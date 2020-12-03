% Create an array that represent a concatenated NetCDF variables.
%
% C = ncCatArray(dim,filenames,varname)
% C = ncCatArray(dim,pattern,varname)
% C = ncCatArray(dim,filenamefun,varname,range)
% C = ncCatArray(...,'property',value)
%
% create a concatenated array from variables (varname) in a list of
% netcdf files along dimension dim.Individual elements can be accessed by
% subscribs, e.g. C(2,3) and the corrsponding subset of the appropriate file is loaded
%
% This list of netcdf files can be specified as a cell array (filenames),
% shell wildcard pattern (e.g. file_*.nc) or a function handle
% filenamefun. In this later case, this i-th filename is
% filenamefun(range(i)).
%
% Properties:
%   'SameAttributes': false or true (default). If SameAttribute is true, 
%     the variables' NetCDF attribute of all files are assumed to be the same.
%     Only the attributes of the first file is loaded in this case.
%
% Example:
%
% data = ncCatArray(3,{'file-20120708.nc','file-20120709.nc'},'SST')
%
% data = ncCatArray(3,'file-*.nc','SST')
%
% data = ncCatArray(3,@(t) ['file-' datestr(t,'yyyymmdd') '.nc'],'SST',...
%              datenum(2012,07,08):datenum(2012,07,09));
%
% Note: in Octave the glob function is used to determine files matching the
% shell wildcard pattern, while in Matlab rdir is used. The function rdir
% is available from Matlab exchange under BSD license
% (http://www.mathworks.com/matlabcentral/fileexchange/19550).
%
% see also cache_decompress, ncArray
% Web: http://modb.oce.ulg.ac.be/mediawiki/index.php/ncArray

% Author: Alexander Barth (barth.alexander@gmail.com)
%
function data = ncCatArray(dim,pattern,varname,varargin)

catdimname = '_cat_dim';
SameAttributes = true;
range = [];

[reg, prop] = parseparams(varargin);

if length(reg) == 1
    range = reg{1};
end


for i = 1:2:length(prop)
  if strcmp(prop{i},'SameAttributes')
    SameAttributes = prop{i+1};
  elseif strcmp(prop{i},'range')
    range = prop{i+1};
  else
    error(['unknown property value ' prop{i}]);
  end
end


% file names

if iscell(pattern)
    filenames = pattern;
    
elseif ischar(pattern)
    try
        filenames = glob(pattern);
    catch
        try
            d = rdir(pattern);
            filenames = {d(:).name};
        catch
                error(['The function rdir or glob (octave) is not available. '...
                'rdir can be installed from '...
                'http://www.mathworks.com/matlabcentral/fileexchange/19550']);
        end
    end
elseif isa(pattern, 'function_handle')
    filenames = cell(1,length(range));
    
    for i=1:length(range)
        filenames{i} = pattern(range(i));
    end
end

if isempty(range)
    range = 1:length(filenames);
end

% get all file information
finfos = cell(length(filenames),1);

if SameAttributes
  % assume all files have the same ncinfo as the first one
  tmp = ncinfo(cached_decompress(filenames{1}));
  for i=1:length(filenames)
    finfos{i} = tmp;
  end
  % octave allows the following, but not matlab
  %finfos(:) = tmp;
else
  for i=1:length(filenames)
    finfos{i} = ncinfo(cached_decompress(filenames{i}));
  end
end

var = arr(dim,filenames,varname,finfos);

[dims,coord] = nccoord(finfos{1},varname);

% add new dimension to coord if arrays are contatenated over a new dimension 
% and if coord information already exist 
if (dim > length(dims)) && ~isempty(coord)
    % concatenate is new dimension
    dims{dim} = catdimname;
    coord(dim).dims = {catdimname};
    coord(dim).val = range;
end


for i=1:length(coord)
    % the number of the dimension might be different
    % find in coord(i).dims the index of the dimension called  dims{dim}
    dimc = find(strcmp(coord(i).dims,dims{dim}));
    
    if isempty(dimc)      
      vinfo = varinfo(finfos{1},coord(i).name);
      coord(i).val = ncBaseArray(filenames{1},coord(i).name,'vinfo',vinfo);
    else    
      % coordinates do also depend on the dimension over which we concatenate
      coord(i).val = arr(dimc,filenames,coord(i).name,finfos);
    end
    
    if dim > length(coord(i).dims)
        coord(i).dims{dim} = catdimname;
    end
end

data = ncArray(var,dims,coord);

end


function CA = arr(dim,filenames,varname,finfos)
arrays = cell(1,length(filenames));

for i=1:length(filenames)
  vinfo = varinfo(finfos{i},varname);
  arrays{i} = ncBaseArray(filenames{i},varname,'vinfo',vinfo);
end

CA = CatArray(dim,arrays);
end

function vinfo = varinfo(fileinfo,varname)
  index = find(strcmp({fileinfo.Variables(:).Name},varname));
  vinfo = fileinfo.Variables(index);
end


% Copyright (C) 2012,2013 Alexander Barth <barth.alexander@gmail.com>
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

