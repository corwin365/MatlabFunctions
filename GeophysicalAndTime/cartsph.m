function Sph = cartsph( Cart, CartV )
% Function to convert cartesian to spherical coordinates
%    Sph = cartsph( Cart, CartV )
% Input arguments:
%    Cart =Cartesian coordinates, [x,y,z;...]
% Optional Input Arguments:
%    CartV = Cartesian vector, [x,y,z;...]
% Output arguments:
%    Sph  = Spherical coordinates or vector, [rad(m),lat(rad),lon(rad);...]
%
% Notes:
%    For geographic coordinates, x points to Greenwich, z to the north pole
% Coordinates and vectors can be passed as cell arrays {x,y,z}. Passing two
% arguments returns a transformed vector
%
% See also SPHCART
% --------------------------------------------------------------------------
if iscell(Cart)
   Sph  = { sqrt(Cart{1}.^2+Cart{2}.^2+Cart{3}.^2),...
            atan2(real(Cart{3}),real(sqrt(Cart{1}.^2+Cart{2}.^2))),...
            atan2(real(Cart{2}),real(Cart{1})) };
else
   Sph  = [ sqrt(Cart(:,1).^2+Cart(:,2).^2+Cart(:,3).^2),...
            atan2(real(Cart(:,3)),real(sqrt(Cart(:,1).^2+Cart(:,2).^2))),...
            atan2(real(Cart(:,2)),real(Cart(:,1)))];
end

if nargin == 1, SphV = []; return; end

if iscell(Cart)
   Sph = {dot(CartV,{cos(Sph{3}).*cos(Sph{2}),sin(Sph{3}).*cos(Sph{2}),sin(Sph{2})}),...
          dot(CartV,{-cos(Sph{3}).*sin(Sph{2}),-sin(Sph{3}).*sin(Sph{2}),cos(Sph{2})}),...
          dot(CartV,{-sin(Sph{3}),cos(Sph{3}),zeros(size(Sph{1}))})};
else
   Sph = [sum(CartV.*[cos(Sph(:,3)).*cos(Sph(:,2)),sin(Sph(:,3)).*cos(Sph(:,2)),sin(Sph(:,2))],2),...
          sum(CartV.*[-cos(Sph(:,3)).*sin(Sph(:,2)),-sin(Sph(:,3)).*sin(Sph{2}),cos(Sph(:,2))],2),...
          sum(CartV.*[-sin(Sph(:,3)),cos(Sph(:,3)),zeros(size(Sph,1),1)],2)];
end