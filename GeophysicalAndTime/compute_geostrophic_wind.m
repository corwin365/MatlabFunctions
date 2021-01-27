function [u,v,u_lat,v_lon] = compute_geostrophic_wind(Data,NoWrap)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute geostrophic wind from GPH data
%
%implements u = -(g/f) .* (dZ./dy)
%and        v = -(g/f) .* (dZ./dx)
%
%
%Inputs:
%  Data - struct, containing:
%     GPH       [nlats x nlons x [arbitrary higher dimensions]] : geopotential HEIGHT
%     LonScale  [nlons]                                         : longitude (degrees)
%     LatScale  [nlats]                                         : latitude  (degrees)
%  NoWrap - if set to 1, longitudes are **not** assumed to wrap (useful for subglobal regions)
%         - to be clear, BY DEFAULT LONGITUDES ARE ASSUMED TO WRAP
%
%
%Outputs:
%     u     [nlats   x nlons-1 x nlevels x ntimes] : zonal wind 
%     v     [nlats-1 x nlons   x nlevels x ntimes] : meridional wind 
%     u_lat [nlats-1]                              : shifted latitude  grid for u
%     v_lon [nlons-1]                              : shifted longitude grid for v
%
%   Corwin Wright,c.wright@bath.ac.uk
%   2021/01/26
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(NoWrap); NoWrap = 0; end
if NoWarp ~= 1;    NoWrap = 0; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% physical constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Const.Omega = 2.0*pi/86400.0;
Const.g     = 9.81;
Const.Re    = 6371e3; %km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fix dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data.LonScale = squeeze(Data.LonScale)';
Data.LatScale = squeeze(Data.LatScale)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% deal with the spherical Earth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NoWrap ~=1;

  %the earth is spherical, so pad either end of the array in the zonal
  %direction by one point
  
  %lonscale is easy
  Data.LonScale = [Data.LonScale(  1)-mean(diff(Data.LonScale)), ...
                   Data.LonScale',                               ...
                   Data.LonScale(end)+mean(diff(Data.LonScale))]';
  
  %gph is trickier as we're working in an arbitrary number of higher dimensions
  sz = size(Data.GPH);
  Data.GPH = reshape(permute(Data.GPH,[2,1,3:numel(sz)]),sz(2),[]);
  Data.GPH = [Data.GPH(end,:);Data.GPH;Data.GPH(1,:)];
  Data.GPH = permute(reshape(Data.GPH,[sz(2)+2,sz(1),sz(3:end)]),[2,1,3:numel(sz)]);
  clear sz
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute dZ/dx and dZ/dy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reshape data to do in one pass
sz = size(Data.GPH);
%need at least three dimensions for the below logic to work
if numel(sz) == 2; sz = [sz,1]; end

%reshape arrays to make the dimension we're taking the derivative the only
%exposed dimension
Z_x = reshape(permute(Data.GPH,[2,1,3:numel(sz)]),sz(2),prod(sz([1,3:end])));
Z_y = reshape(        Data.GPH,                   sz(1),prod(sz([  2:end])));

%take derivative, in metres
dZdx = diff(Z_x,1,1)./ (diff(Data.LonScale) .* (2 .* pi .* Const.Re ./ 360));
dZdy = diff(Z_y,1,1)./ (diff(Data.LatScale) .* (2 .* pi .* Const.Re ./ 360));
clear Z_x Z_y

%reshape back
dZdy =         reshape(dZdy,[sz(1)-1,sz(2),sz(3:end)]);
dZdx = permute(reshape(dZdx,[sz(2)-1,sz(1),sz(3:end)]),[2,1,3:numel(sz)]);
clear sz



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% deal with the spherical Earth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NoWrap ~=1;
  
  %lonscale is again easy
  Data.LonScale = Data.LonScale(2:end-1);
  
  %dZdy and dZdx are again hard
  sz = size(dZdx);
  dZdx = reshape(permute(dZdx,[2,1,3:numel(sz)]),sz(2),[]);
  dZdx = dZdx(2:end-1,:);
  dZdx = permute(reshape(dZdx,[sz(2)-2,sz(1),sz(3:end)]),[2,1,3:numel(sz)]);
  
  sz = size(dZdy);
  dZdy = reshape(permute(dZdy,[2,1,3:numel(sz)]),sz(2),[]);
  dZdy = dZdy(2:end-1,:);
  dZdy = permute(reshape(dZdy,[sz(2)-2,sz(1),sz(3:end)]),[2,1,3:numel(sz)]);

  clear sz
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scale to wind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fa = 2.* Const.Omega .* sind(Data.LatScale(1:end-1)+diff(Data.LatScale)./2);
fb = 2.* Const.Omega .* sind(Data.LatScale);

%reshape again, exposing latitude, then compute u and v
szX = size(dZdx); szY = size(dZdy);
v = reshape(((-Const.g./fb) .* reshape(dZdx,szX(1),[])),szX);
u = reshape(((-Const.g./fa) .* reshape(dZdy,szY(1),[])),szY);
clear szX szY fa fb dZdx dZdy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fix scales, and done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_lon = Data.LonScale(1:end-1)+diff(Data.LonScale)./2;
u_lat = Data.LatScale(1:end-1)+diff(Data.LatScale)./2;

return