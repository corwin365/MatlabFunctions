% Test ncBaseArray, ncCatArray and ncArray.
function test_ncarray_nan()

% for octave prior to 3.8.0
if isempty(which('isequaln'))
  isequaln = @(x,y) isequalwithequalnans(x,y);
end

varname = 'SST';

tmpdir = tempname;
mkdir(tmpdir);

tmpfname = tempname(tmpdir);
dataref = randn(220,144,3);
%dataref(rand(size(dataref)) > 0.7) = NaN;
dataref(50:80,30:90,1:2) = NaN;

for i = 1:3  
  files{i} = fullfile(tmpdir,sprintf('file%d.nc',i));
  ncarray_example_file(files{i},dataref(:,:,i));
end

data = ncCatArray(3,files,varname);
reddata = nanmean(data,3);
reddataref = nanmean(dataref,3);
assert(isequaln(reddata, reddataref))

reddata = nansum(data,3);
reddataref = nansum(dataref,3);
assert(isequaln(reddata, reddataref))

reddata = nanvar(data,[],3);
reddataref = nanvar(dataref,[],3);
diff = reddata - reddataref;
assert(max(diff(:)) < 1e-6)

reddata = nanstd(data,[],3);
reddataref = nanstd(dataref,[],3);
diff = reddata - reddataref;
assert(max(diff(:)) < 1e-6)

% clean-up
for i = 1:3  
  delete(files{i});
end
rmdir(tmpdir);




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

