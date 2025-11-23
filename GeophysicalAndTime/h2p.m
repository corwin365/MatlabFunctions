
function POut = h2p(alt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert pressures to altitudes
%
%inputs:
%   alt - altitude levels, in km (any size of array is fine)
%outputs:
%   p - atmospheric pressure values, in hPa
%
%Corwin Wright, c.wright@bath.ac.uk, 2025/04/02 (full rewrite of older 2023 version to properly vectorise,
%which was itself based on much older code by Neil)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%define pressure and altitude bounds for each level
p = [1013.25, 226.3210, 54.7489,  8.6802,  1.1091,  0.6694,  0.0396];
z = [   0,     11,      20,      32,      47,      51,      71];  

% calculate scale heights within each level
H = -diff(z)./log(p(2:end)./p(1:end-1));  % Scale heights (km)

%create output array
POut = NaN(size(alt)); 

%apply to each altitude **band**
for i = 1:length(z)-1
  idx = (alt >= z(i) & alt < z(i+1));                 % find altitudes in range
  POut(idx) = p(i) * exp(-(alt(idx) - z(i)) / H(i));  %compute output pressure
end

%deal with high altitudes (above 71km)
POut(alt >= z(end)) = p(end) * exp(-(alt(alt >= z(end)) - z(end)) / H(end));


end