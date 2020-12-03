function ncarray_example_file(filename,data)

if ~isempty(which('nccreate')) && ~isempty(which('ncwriteatt')) && ...
      ~isempty(which('ncwrite'))
  % use matlab netcdf high level interface

  %dtype = 'single';
  dtype = 'double';

  % Variables
  nccreate(filename,'lon','Format','classic','Datatype',dtype,...
           'Dimensions',{'x',220, 'y',144});
  ncwriteatt(filename,'lon','long_name','Longitude')
  ncwriteatt(filename,'lon','units','degrees_east')
  
  nccreate(filename,'lat','Datatype',dtype,'Dimensions',{'x',220, 'y',144});
  ncwriteatt(filename,'lat','long_name','Latitude')
  ncwriteatt(filename,'lat','units','degrees_north')
  
  nccreate(filename,'time','Datatype',dtype,'Dimensions',{'time',1});
  ncwriteatt(filename,'time','long_name','Time')
  ncwriteatt(filename,'time','units','days since 1858-11-17 00:00:00 GMT')
  
  nccreate(filename,'SST','Datatype',dtype,'Dimensions',...
           {'x',220, 'y',144, 'time',1});
  ncwriteatt(filename,'SST','missing_value',single(9999))
  ncwriteatt(filename,'SST','_FillValue',single(9999))
  ncwriteatt(filename,'SST','units','degC')
  ncwriteatt(filename,'SST','long_name','Sea Surface Temperature')
  ncwriteatt(filename,'SST','coordinates','lat lon')

  ncwrite(filename,'SST',data);  
else
  % use octcdf interface

  nc = netcdf(filename,'c');

  % dimensions
  
  nc('x') = size(data,1);
  nc('y') = size(data,2);
  nc('time') = size(data,3);

  % variables

  nc{'lon'} = ncdouble('y','x');  % 31680 elements 
  nc{'lon'}.long_name = ncchar('Longitude');
  nc{'lon'}.units = ncchar('degrees_east');
  
  nc{'lat'} = ncdouble('y','x');  % 31680 elements 
  nc{'lat'}.long_name = ncchar('Latitude');
  nc{'lat'}.units = ncchar('degrees_north');

  nc{'time'} = ncdouble('time');  % 1 elements 
  nc{'time'}.long_name = ncchar('Time');
  nc{'time'}.units = ncchar('days since 1858-11-17 00:00:00 GMT');
  
  nc{'SST'} = ncdouble('time','y','x');  % 31680 elements 
  nc{'SST'}.missing_value = ncdouble(9999);
  nc{'SST'}.FillValue_ = ncdouble(9999);
  nc{'SST'}.units = ncchar('degC');
  nc{'SST'}.long_name = ncchar('Sea Surface Temperature');
  nc{'SST'}.coordinates = ncchar('lat lon');

  % global attributes
  
  nc{'SST'}(:) = permute(data,[3 2 1]);
  close(nc)
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

