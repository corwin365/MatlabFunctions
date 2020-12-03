% Parse netCDF time unit.
%
% [t0,f] = nctimeunits(u)
%
% Parse the netCDF time unit u and returns the time offset (days since 31 
% December 1 BC, as datenum) and scaling factor f (in days).
% See the netCDF CF convention for the structure of the time units.

function [t0,f] = nctimeunits(u)

if strfind(u,'seconds')
  f = 1/(24*60*60);
elseif strfind(u,'hours')
  f = 1/24;
elseif strfind(u,'days')
  f = 1;
else
  error(['unknown units "' u '"']);
end
  
l = strfind(u,'since')+6;

try
  t0 = datenum(u(l:end),'yyyy-mm-dd HH:MM:SS');
catch
  try
    t0 = datenum(u(l:end),'yyyy-mm-ddTHH:MM:SS');
  catch    
    try
      t0 = datenum(u(l:end),'yyyy-mm-dd');
    catch
      error(['date format is not recogized ' u(l:end)])
    end
  end
end
