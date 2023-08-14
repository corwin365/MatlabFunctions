function Sat = invelevation( Rec, Azim, Elev, RSat )
% Compute satellite coords given azimuth and elevation to satellite
%   Sat = invelevation( Rec, Azim, Elev, RSat )
% Input Arguments:
%   Rec  = Cartesian coordinates of the receiver (m)
%   Az   = Azimuth to satellite (radians)
%   El   = Elevation of satellite (radians)
%   RSat = Height of satellite orbit (above Earth) (m)
% Output Arguments:
%   Sat = Cartesian satellite coordinates (m)
%-----------------------------------------------------------------------


E = cross(repmat([0,0,1],size(Rec,1),1),Rec,2); % East
N = cross(Rec,E,2);                             % North
V = Rec;                                        % Vertical

MagE = sqrt(sum( (E.^2).' ).'); E = E./repmat(MagE,1,3);
MagN = sqrt(sum( (N.^2).' ).'); N = N./repmat(MagN,1,3);
MagV = sqrt(sum( (V.^2).' ).'); V = V./repmat(MagV,1,3);

Sat =  N.*repmat(cos(Azim).*cos(Elev),1,3) + ...
       E.*repmat(sin(Azim).*cos(Elev),1,3) + ...
       V.*repmat(sin(Elev),1,3);

P = -MagV.*sin(Elev) + sqrt( -(MagV.*cos(Elev)).^2 + RSat^2 );

Sat = Rec + Sat.*repmat(P,1,3);


