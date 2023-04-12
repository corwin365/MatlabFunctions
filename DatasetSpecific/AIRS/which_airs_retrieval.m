function Mask = which_airs_retrieval(Lon,Lat,Time,RowFlag)

%identified whether AIRS-3D is using the day or night retrieval, using the same
%logic as Lars' retrieval does itself (shown in C at bottom of function)
%
%Corwin Wright, c.wright@bath.ac.uk, 14/AUG/2018
%
%modified 2023/03/01 to add the option to give the option of switching at the
%row level for data, rather than the specific point. Use this by either setting
%RowFlag to specific dimension number (to set all elements along other dimensions
%to the majority choice in this dimension) or to -1 to do this in the 
%first-ordered dimension of length 90 (which may be easier in some cases)
%
%
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



%now, handle majoritarian rowflagging
if exist('RowFlag')

  %get size of mask
  sz = size(Mask);
  
  if RowFlag == -1;
    %set dimension to first dimension of length 90
    [delta,RowFlag] = min(abs(sz-90));
    if delta ~= 0;
      disp('no dimensions of length 90, majoritarian rowflagging not applied')
      return
    end
  end
 
  if RowFlag < 1 | RowFlag > numel(sz);
      disp('invalid dimension, majoritarian rowflagging not applied')
      return
  end

  %ok, we have a valid dimension. reshape so we can find the majority in a single pass
  order = unique([RowFlag,1:1:numel(sz)],'stable');
  M = reshape(permute(Mask,order),[sz(order(1)),prod(sz(order(2:end)))]);
  Sigma = sum(M,1);
  One  = find(Sigma >= size(M,1)./2);
  Zero = find(Sigma <  size(M,1)./2);  
  M(:,One) = 1;
  M(:,Zero) = 0;
  
  Mask = reshape(M,[sz(order(1)),sz(order(2:end))]);
  NewOrder = 1:1:numel(sz);
  NewOrder = [NewOrder(NewOrder < RowFlag)+1,1,NewOrder(NewOrder > RowFlag)];
  Mask = permute(Mask,NewOrder);

end


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
