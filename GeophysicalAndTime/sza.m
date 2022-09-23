function SZA = sza(Lon,Lat,Time)

%find solar zenith angle of a given lat,lon, and time

%Number of days and fraction with respect to 2000-01-01T12:00Z
D = Time-datenum(2000,1,1,12,0,0);

%Geocentric apparent ecliptic longitude [rad]
g = (357.529 + 0.98560028 .* D) .* pi ./ 180;
q = 280.459 + 0.98564736 .* D;
L = (q + 1.915 .* sin(g) + 0.020 .* sin(2 * g)) .* pi ./ 180;

%Mean obliquity of the ecliptic [rad]
e = (23.439 - 0.00000036 * D) .* pi ./ 180;

% Declination [rad]
dec = asin(sin(e) .* sin(L));

% Right ascension [rad]...
ra = atan2(cos(e) .* sin(L), cos(L));
 
% Greenwich Mean Sidereal Time [h]...
GMST = 18.697374558 + 24.06570982441908 .* D;


%Local Sidereal Time [h]... 
LST = GMST + Lon ./ 15;

% Hour angle [rad]... 
h = LST ./ 12 .* pi - ra;
 
% Convert latitude... 
Lat = Lat .* pi ./ 180;

%hence, solar zenith angle [deg].
SZA = acos(sin(Lat) .* sin(dec) + cos(Lat) .* cos(dec) .* cos(h)) .* 180 / pi;
