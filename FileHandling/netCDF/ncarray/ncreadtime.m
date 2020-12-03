% Read a time variable as serial day number.
% 
% t = ncreadtime(filename,varname)
%
% Read a time variable called varname from the file called filename as serial 
% day number (days since 31 December 1 BC, as datenum).

function t = ncreadtime(filename,varname)

t = ncread(filename,varname);
units = ncreadatt(filename,varname,'units');

[t0,f] = nctimeunits(units);

t = f*t + t0;
