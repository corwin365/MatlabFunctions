function Goes = goes_latlon(Goes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute latitude and longitude for each point in GOES data
%
%based on: https://makersportal.com/blog/2018/11/25/goes-r-satellite-latitude-and-longitude-grid-projection-algorithm
%
%Corwin Wright, c.wright@bath.ac.uk
%12/FEB/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GOES-R projection info 
Req  = 6378137.0;
Rpol = 6356752.31414;
L0   = -75;
Hsat = 35786023.0;
H    = Req + Hsat;

%Data info
lat_rad_1d = Goes.x;
lon_rad_1d = Goes.y;

%reate meshgrid filled with radian angles
[lat_rad,lon_rad] = meshgrid(lat_rad_1d,lon_rad_1d);

%lat/lon calc routine from satellite radian angle vectors
lambda_0 = (L0.*pi)/180.0;

a_var = sin(lat_rad).^2 + (cos(lat_rad).^2).*(cos(lon_rad).^2) ...
                        +(((Req.*Req)./(Rpol.*Rpol))*sin(lon_rad).^2);
b_var = -2.0.*H.*cos(lat_rad).*cos(lon_rad);
c_var = (H.^2)-(Req.^2);

r_s = (-1.0.*b_var - sqrt((b_var.^2)-(4.0.*a_var.*c_var)))./(2.0.*a_var);

s_x = r_s.*cos(lat_rad).*cos(lon_rad);
s_y = - r_s.*sin(lat_rad);
s_z = r_s.*cos(lat_rad).*sin(lon_rad);

lat = (180.0/pi)*(atan(((Req.*Req)/(Rpol.*Rpol)).*((s_z./sqrt(((H-s_x).*(H-s_x))+(s_y.*s_y))))));
lon = (lambda_0 - atan(s_y./(H-s_x))).*(180.0./pi);


%remove bad values (outside the Earth's disk)
IsComplex = find(imag(lat) ~= 0); lat(IsComplex) = NaN;
IsComplex = find(imag(lon) ~= 0); lon(IsComplex) = NaN;

%store
Goes.Lat = lat;
Goes.Lon = lon;

end