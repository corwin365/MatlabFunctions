function [u,v,u_lon,v_lat] = compute_geostrophic_wind(Data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute geostrophic wind from GPH data
%
%implements u = -(1/f) .* (dZ./dy)
%and        v = -(1/f) .* (dZ./dx)
%
%
%Inputs:
%  Data - struct, containing:
%     GPH       [nlats x nlons x [arbitrary higher dimensions]] : geopotential HEIGHT
%     LonScale  [nlons]                                         : longitude (degrees), full range required
%     LatScale  [nlats]                                         : latitude  (degrees), arbitrary range
%
%
%Outputs:
%     u     [nlons-1 x nlats   x nlevels x ntimes] : zonal wind 
%     v     [nlons   x nlats-1 x nlevels x ntimes] : meridional wind 
%     u_lon [nlons-1]                              : shifted longitude grid for u
%     v_lat [nlats-1]                              : shifted latitude  grid for v
%
%   Corwin Wright,c.wright@bath.ac.uk
%   2021/01/26
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%the earth is spherical, so pad either end of the array in the zonal
%direction by one point

%lon is easier
Data.LonScale = [Data.LonScale(  1)-mean(diff(Data.LonScale)), ...
                 Data.LonScale',                               ...
                 Data.LonScale(end)+mean(diff(Data.LonScale))]';
               
%gph is trickier as we're working in an arbitrary number of higher dimensions
sz = size(Data.GPH);
Data.GPH = permute(Data.GPH,[2,1,3:numel(sz)]);
Data.GPH = reshape(Data.GPH,sz(2),[]);
Data.GPH = [Data.GPH(end,:);Data.GPH;Data.GPH(1,:)];
Data.GPH = reshape(Data.GPH,[sz(2)+2,sz(1),sz(3:end)]);
Data.GPH = permute(Data.GPH,[2,1,3:numel(sz)]);
clear sz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute dZ/dx and dZ/dy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reshape data to do in one pass
sz = size(Data.GPH);
%need at least three dimensions for the below logic to work
if numel(sz) == 2; sz = [sz,1]; end

%reshape arrays to make the dimension we're taking the derivative the only
%exposed dimension
Z_x = reshape(permute(Data.GPH,[2,1,3:end]),sz(2),prod(sz([1,3:end])));
Z_y = reshape(        Data.GPH,             sz(1),prod(sz([  2:end])));

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

%lonscale is again easy
Data.LonScale = Data.LonScale(2:end-1);

%dZdy and dZdx are again hard
sz = size(dZdx);
dZdx = permute(dZdx,[2,1,3:numel(sz)]);
dZdx = reshape(dZdx,sz(2),[]);
dZdx = dZdx(2:end-1,:);
dZdx = reshape(dZdx,[sz(2)-2,sz(1),sz(3:end)]);
dZdx = permute(dZdx,[2,1,3:numel(sz)]);

sz = size(dZdy);
dZdy = permute(dZdy,[2,1,3:numel(sz)]);
dZdy = reshape(dZdy,sz(2),[]);
dZdy = dZdy(2:end-1,:);
dZdy = reshape(dZdy,[sz(2)-2,sz(1),sz(3:end)]);
dZdy = permute(dZdy,[2,1,3:numel(sz)]);

clear sz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scale to wind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fv = 2.* Const.Omega .* sind(Data.LatScale(1:end-1)+diff(Data.LatScale)./2);
fu = 2.* Const.Omega .* sind(Data.LatScale);

%reshape again, exposing latitude, then compute u and v
szX = size(dZdx); szY = size(dZdy);
u = reshape(((-Const.g./fu) .* reshape(dZdx,szX(1),[])),szX);
v = reshape(((-Const.g./fv) .* reshape(dZdy,szY(1),[])),szY);
clear szX szY fu fv dZdx dZdy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fix scales, and done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_lon = Data.LonScale(1:end-1)+diff(Data.LonScale)./2;
v_lat = Data.LatScale(1:end-1)+diff(Data.LatScale)./2;

return