function [E,A] = elevation( Sat, Rec )
% Function to compute the elevation angle and azimuth of a satellite at a receiver
%   [E,A] = elevation( Sat, Rec )
% Input Arguments:
%   Sat = cartesian coordinates of the satellite
%   Rec = cartesian coordinates of the receiver
% Output Arguments:
%   E   = Elevation angle in radians (from horizontal)
%   A   = Azimuth angle from north in radians
% NB Generate a polar plot using polar(A,E*180/pi);
%-----------------------------------------------------------------------

RS = Sat-Rec;
OR = Rec;

MagRS = sqrt(sum( (RS.^2).' ).');
MagOR = sqrt(sum( (OR.^2).' ).');

E = pi/2-acos(dot(RS,OR,2)./(MagRS.*MagOR));

Horiz = cross(cross(OR,RS,2),OR,2);
MagH  = sqrt(sum( (Horiz.^2).' ).');

North = cartsph(Rec); 
North = [North(:,1) North(:,2)+0.01 North(:,3)];
North = sphcart(North); North = North - Rec;
East  = cartsph(Rec); 
East  = [East(:,1) East(:,2) East(:,3)+0.01];
East  = sphcart(East); East = East - Rec;

A = atan2(dot(East,Horiz,2),dot(North,Horiz,2));

