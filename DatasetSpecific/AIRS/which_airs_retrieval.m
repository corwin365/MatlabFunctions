function Mask = which_airs_retrieval(Lon,Lat,Time)

%identified whether AIRS-3D is using the day or night retrieval, using the same
%logic as Lars' retrieval does itself (shown in C at bottom of function)
%
%Corwin Wright, c.wright@bath.ac.uk, 14/AUG/2018
%
%input is time in Matlab format, latitude, and longitude (array-safe, must all
%be same shape)
%output is a binary mask, with 1 for daytime and 0 for nighttime, of the same
%size as the inputs

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
SZAd = acos(sin(Lat) .* sin(dec) + cos(Lat) .* cos(dec) .* cos(h)) .* 180 / pi;

%create mask
Mask = SZAd .*NaN;
Mask(SZAd <  96) = 1; %day
Mask(SZAd >= 96) = 0; %night



end

%
%original C version:
%
% double sza(
%   double sec,
%   double lon,
%   double lat) {
% 
%   double D, dec, e, g, GMST, h, L, LST, q, ra;
% 
%   /* Number of days and fraction with respect to 2000-01-01T12:00Z... */
%   D = sec / 86400 - 0.5;
% 
%   /* Geocentric apparent ecliptic longitude [rad]... */
%   g = (357.529 + 0.98560028 * D) * M_PI / 180;
%   q = 280.459 + 0.98564736 * D;
%   L = (q + 1.915 * sin(g) + 0.020 * sin(2 * g)) * M_PI / 180;
% 
%   /* Mean obliquity of the ecliptic [rad]... */
%   e = (23.439 - 0.00000036 * D) * M_PI / 180;
% 
%   /* Declination [rad]... */
%   dec = asin(sin(e) * sin(L));
% 
%   /* Right ascension [rad]... */
%   ra = atan2(cos(e) * sin(L), cos(L));
% 
%   /* Greenwich Mean Sidereal Time [h]... */
%   GMST = 18.697374558 + 24.06570982441908 * D;
% 
%   /* Local Sidereal Time [h]... */
%   LST = GMST + lon / 15;
% 
%   /* Hour angle [rad]... */
%   h = LST / 12 * M_PI - ra;
% 
%   /* Convert latitude... */
%   lat *= M_PI / 180;
% 
%   /* Return solar zenith angle [deg]... */
%   return acos(sin(lat) * sin(dec) +
%           cos(lat) * cos(dec) * cos(h)) * 180 / M_PI;
% }
% 
% /*****************************************************************************/
% 
% daytime data (15 micron channels only): SZA < 96°
% 
% nighttime data (4 + 15 micron channels used): SZA >= 96°
