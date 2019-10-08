% 
% loc1 = [63.2 -30 ; 63.2 0 ; 63.2 45];
% loc2 = [13.95 150; 13.95 180 ; 13.95 225];
% 
% 


function [km] = nph_haversine_mars(loc1,loc2)


% If your input is an array of lats/lons, it should be of the form:

% loc1 = [N_POSITIONS x LATLON],
% e.g. for 100 positions loc1 = <100x2 double>.


% % HAVERSINE     Compute distance between locations using Haversine formula
% %   KM = HAVERSINE(LOC1, LOC2) returns the distance KM in km between
% %   locations LOC1 and LOC2 using the Haversine formula.  LOC1 and LOC2 are
% %   latitude and longitude coordinates.
% %   
% %   Inputs
% %       loc1 = [LAT1 LON1]
% %       loc2 = [LAT2 LON2]
% % 
% %   Notes
% %       The Haversine formula is used to calculate the great-circle
% %       distance between two points, which is the shortest distance over
% %       the earth's surface.
% % 
% % For a cosmic occultation at h_cosmic ~ 800km with a GPS satellite, a
% % haversine distance of ~11436.398km is found. You can check with geometry
% % the following:
% % loc1 = [63.2 -30];
% % loc2 = [13.95 150];
% % 
% % EDIT - made array safe. I think.
% % 

%%

loc1 = deg2rad(loc1); loc2 = deg2rad(loc2);

R = 3390;                                 % Marsradius in km
delta_lat = loc2(:,1) - loc1(:,1);        % difference in latitude
delta_lon = loc2(:,2) - loc1(:,2);        % difference in longitude
a = sin(delta_lat./2).^2 + cos(loc1(:,1)) .* cos(loc2(:,1)) .* ...
    sin(delta_lon./2).^2;
c = 2 .* atan2(sqrt(a), sqrt(1-a));
km = R .* c;

% 
% %%
% 
% % Get azimuth too? BROKEN I THINK :(
% 
% lat1 = loc1(:,1); lon1 = loc1(:,2);
% lat2 = loc2(:,1); lon2 = loc2(:,2);
% 
% 
% az = atan2(cos(lat2) .* sin(lon2-lon1),...
%            cos(lat1) .* sin(lat2) - sin(lat1) .* cos(lat2) .* cos(lon2-lon1));
% 
% % Azimuths are undefined at the poles, so we choose a convention: zero at
% % the north pole and pi at the south pole.
% az(lat1 <= -pi/2) = 0;
% az(lat2 >=  pi/2) = 0;
% az(lat2 <= -pi/2) = pi;
% az(lat1 >=  pi/2) = pi;
% 

end
