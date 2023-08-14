function [Az h] = LunarAzEl(UTC,Lat,Lon,Alt)

% Programed by Darin C. Koblick 2/14/2009
%
% Updated on 03/04/2009 to clean up code and add quadrant check to Azimuth
% Thank you Doug W. for your help with the test code to find the quadrant check 
% error.
%
% Updated on 04/13/2009 to add Lunar perturbation offsets (this will
% increase the accuracy of calculations)
%
% Updated on 08/17/2010 to make use of the site altitude (this will affect
% the elevation angle)

% External Function Call Sequence:
% [Az El] = LunarAzEl('1991/05/19 13:00:00',50,10,0)


% Function Description:
% LunarAzEl will ingest a Universal Time, and specific site location on earth
% it will then output the lunar Azimuth and Elevation angles relative to that
% site.

% External Source References:

% Basics of Positional Astronomy and Ephemerides
% http://jgiesen.de/elevazmoon/basics/index.htm

% Computing planetary positions - a tutorial with worked examples
% http://stjarnhimlen.se/comp/tutorial.html

%Input Description:
% UTC (Coordinated Universal Time YYYY/MM/DD hh:mm:ss)
% Lat (Site Latitude in degrees -90:90 -> S(-) N(+))
% Lon (Site Longitude in degrees -180:180 W(-) E(+))
% Altitude of the site above sea level (km)

%Output Description:
%Az (Azimuth location of the moon in degrees)
%El (Elevation location of the moon in degrees)

%Verified output by comparison with the following source data:
%http://aa.usno.navy.mil/data/docs/AltAz.php

% Code Sequence
%--------------------------------------------------------------------------
%Do initial Longitude Latitude check

while Lon > 180
   Lon = Lon - 360;
end
while Lon < -180
   Lon = Lon + 360; 
end
while Lat > 90
    Lat = Lat - 360;
end
while Lat < -90
   Lat = Lat + 360; 
end


%Declare Earth Equatorial Radius Measurements in km
EarthRadEq = 6378.1370;

%Convert Universal Time to Ephemeris Time
jd = juliandate(UTC);%,'yyyy/mm/dd HH:MM:SS');

%Find the Day Number
d = jd - 2451543.5;

%Keplerian Elements of the Moon
%This will also account for the Sun's perturbation
 N = 125.1228-0.0529538083.*d; %    (Long asc. node deg)
 i = 5.1454; %                      (Inclination deg)
 w = 318.0634 + 0.1643573223.*d; %  (Arg. of perigee deg)
 a =  60.2666;%                     (Mean distance (Earth's Equitorial Radii)
 e = 0.054900;%                     (Eccentricity)
 M = mod(115.3654+13.0649929509.*d,360);%    (Mean anomaly deg)
  
 LMoon =  mod(N + w + M,360);                 %(Moon's mean longitude deg)
 FMoon =  mod(LMoon - N,360);                 %(Moon's argument of latitude)

 %Keplerian Elements of the Sun
 wSun = mod(282.9404 + 4.70935E-5.*d,360);    % (longitude of perihelion)
 MSun = mod(356.0470 + 0.9856002585.*d,360);  % (Sun mean anomaly)
 LSun = mod(wSun + MSun,360);                 % (Sun's mean longitude)
     
 DMoon =  LMoon - LSun;                     % (Moon's mean elongation)  


 %Calculate Lunar perturbations in Longitude
 LunarPLon = [ -1.274.*sin((M - 2.*DMoon).*(pi/180)); ...
     .658.*sin(2.*DMoon.*(pi/180)); ...
     -0.186.*sin(MSun.*(pi/180)); ...
     -0.059.*sin((2.*M-2.*DMoon).*(pi/180)); ...
     -0.057.*sin((M-2.*DMoon + MSun).*(pi/180)); ...
     .053.*sin((M+2.*DMoon).*(pi/180)); ...
     .046.*sin((2.*DMoon-MSun).*(pi/180)); ...
     .041.*sin((M-MSun).*(pi/180)); ...
    -0.035.*sin(DMoon.*(pi/180)); ...           
    -0.031.*sin((M+MSun).*(pi/180)); ...
    -0.015.*sin((2.*FMoon-2.*DMoon).*(pi/180)); ...
    .011.*sin((M-4.*DMoon).*(pi/180))];
 
 %Calculate Lunar perturbations in Latitude 
 LunarPLat = [ -0.173.*sin((FMoon-2.*DMoon).*(pi/180)); ...
    -0.055.*sin((M-FMoon-2.*DMoon).*(pi/180)); ...
    -0.046.*sin((M+FMoon-2.*DMoon).*(pi/180)); ...
    +0.033.*sin((FMoon+2.*DMoon).*(pi/180)); ...
    +0.017.*sin((2.*M+FMoon).*(pi/180))];

%Calculate perturbations in Distance
 LunarPDist = [ -0.58*cos((M-2.*DMoon).*(pi/180)); ...
    -0.46.*cos(2.*DMoon.*(pi/180))];

% Compute E, the eccentric anomaly

%E0 is the eccentric anomaly approximation estimate 
%(this will initially have a relativly high error)
E0 = M+(180./pi).*e.*sin(M.*(pi/180)).*(1+e.*cos(M.*(pi/180)));

%Compute E1 and set it to E0 until the E1 == E0
E1 = E0-(E0-(180/pi).*e.*sin(E0.*(pi/180))-M)./(1-e*cos(E0.*(pi/180)));
while E1-E0 > .000005
    E0 = E1;
    E1 = E0-(E0-(180/pi).*e.*sin(E0.*(pi/180))-M)./(1-e*cos(E0.*(pi/180)));    
end
E = E1;

%Compute rectangular coordinates (x,y) in the plane of the lunar orbit
x = a.*(cos(E.*(pi/180))-e);
y = a.*sqrt(1-e.*e).*sin(E.*(pi/180));

%convert this to distance and true anomaly
r = sqrt(x.*x + y.*y);
v = atan2(y.*(pi/180),x.*(pi/180)).*(180/pi);

%Compute moon's position in ecliptic coordinates
xeclip = r.*(cos(N.*(pi/180)).*cos((v+w).*(pi/180))-sin(N.*(pi/180)).*sin((v+w).*(pi/180)).*cos(i.*(pi/180)));
yeclip = r.*(sin(N.*(pi/180)).*cos((v+w).*(pi/180))+cos(N.*(pi/180))*sin(((v+w).*(pi/180)))*cos(i.*(pi/180)));
zeclip = r.*sin((v+w).*(pi/180)).*sin(i.*(pi/180));

%Add the calculated lunar perturbation terms to increase model fidelity
[eLon eLat eDist] = cart2sph(xeclip,yeclip,zeclip);
[xeclip yeclip zeclip] = sph2cart(eLon + sum(LunarPLon).*(pi/180), ...
                                  eLat + sum(LunarPLat).*(pi/180), ...
                                  eDist + sum(LunarPDist));
clear eLon eLat eDist;
                              
%convert the latitude and longitude to right ascension RA and declination
%delta
T = (jd-2451545.0)/36525.0;

%Generate a rotation matrix for ecliptic to equitorial
%RotM=rotm_coo('E',jd);
%See rotm_coo.m for obl and rotational matrix transformation
Obl = 23.439291 - 0.0130042.*T - 0.00000016.*T.*T + 0.000000504.*T.*T.*T;
Obl = Obl.*(pi/180);
RotM = [1 0 0; 0 cos(Obl) sin(Obl); 0 -sin(Obl) cos(Obl)]';

%Apply the rotational matrix to the ecliptic rectangular coordinates
%Also, convert units to km instead of earth equatorial radii
sol = RotM*[xeclip yeclip zeclip]'.*EarthRadEq;

%Find the equatorial rectangular coordinates of the location specified
[xel yel zel] = sph2cart(Lon.*(pi/180),Lat.*(pi/180),Alt+EarthRadEq);

%Find the equatorial rectangular coordinates of the location @ sea level
[xsl ysl zsl] = sph2cart(Lon.*(pi/180),Lat.*(pi/180),EarthRadEq);

%Find the Angle Between sea level coordinate vector and the moon vector
theta1 = 180 - acosd(dot([xsl ysl zsl],[sol(1)-xsl sol(2)-ysl sol(3)-zsl]) ...
        ./(sqrt(xsl.^2 + ysl.^2 + zsl.^2) ...
         .*sqrt((sol(1)-xsl).^2 + (sol(2)-ysl).^2 + (sol(3)-zsl).^2)));

%Find the Angle Between the same coordinates but at the specified elevation
theta2 = 180 - acosd(dot([xel yel zel],[sol(1)-xel sol(2)-yel sol(3)-zel]) ...
    ./(sqrt(xel.^2 + yel.^2 + zel.^2) ...
         .*sqrt((sol(1)-xel).^2 + (sol(2)-yel).^2 + (sol(3)-zel).^2)));
     
%Find the Difference Between the two angles (+|-) is important
thetaDiff = theta2 - theta1;

% equatorial to horizon coordinate transformation
 [RA,delta] = cart2sph(sol(1),sol(2),sol(3));
 delta = delta.*(180/pi);
 RA = RA.*(180/pi);
 
%Following the RA DEC to Az Alt conversion sequence explained here:
%http://www.stargazing.net/kepler/altaz.html

%Find the J2000 value
J2000 = jd - 2451545.0;
warning off
hourvec = datevec(UTC,'yyyy/mm/dd HH:MM:SS');
warning on
UTH = hourvec(4) + hourvec(5)/60 + hourvec(6)/3600;

%Calculate local siderial time
LST = mod(100.46+0.985647.*J2000+Lon+15*UTH,360);

%Replace RA with hour angle HA
HA = LST-RA;

%Find the h and AZ at the current LST
h = asin(sin(delta.*(pi/180)).*sin(Lat.*(pi/180)) + cos(delta.*(pi/180)).*cos(Lat.*(pi/180)).*cos(HA.*(pi/180))).*(180/pi);
Az = acos((sin(delta.*(pi/180)) - sin(h.*(pi/180)).*sin(Lat.*(pi/180)))./(cos(h.*(pi/180)).*cos(Lat.*(pi/180)))).*(180/pi);

%Add in the angle offset due to the specified site elevation
h = h + thetaDiff;

if sin(HA.*(pi/180)) >= 0
   Az = 360-Az; 
end

%Apply Paralax Correction if we are still on earth
if Alt < 100
    horParal = 8.794/(r*6379.14/149.59787e6);
    p = asin(cos(h.*(pi/180))*sin((horParal/3600).*(pi/180))).*(180/pi);
    h = h-p;
end

function jd = juliandate(varargin)
% This sub function is provided in case juliandate does not come with your 
% distribution of Matlab

[year month day hour min sec] = datevec(datenum(varargin{:}));

for k = length(month):-1:1
    if ( month(k) <= 2 ) % january & february
        year(k)  = year(k) - 1.0;
        month(k) = month(k) + 12.0;
    end
end

jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
    floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + ...
    (hour + min/60 + sec/3600)/24;












