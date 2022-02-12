function [regridded_image,x_in,y_in,Lat,Lon] = regridairglowimage(Image,cameralat,cameralon,h,newgridX,newgridY)

%%%% INPUTS
% Image                 = image data
% cameralat,cameralon   = the lat lon of the camera
% h                     = altitude of the airglow layer in KM
% newgridX,newgridY     = new regular output grid in KM created by ndgrid or meshgrid etc.

%%%% OUTPUTS
% regridded_image   = the image on the new regular grid (newgridX,newgridY)
% Lat,Lon           = the lats and lons of this new regular grid.

% input image
Z = double(flip(Image,1));
sz = size(Z);

% Radius of the earth
Re = 6371; % km

% radius of the airglow layer
r = Re + h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Find the zenith angle and azimuth of every location in the image

% try a different method:
[V1,V2] = meshgrid(linspace(-1,1,sz(1)),linspace(-1,1,sz(2)));
zen     = 90.*sqrt(V1.^2+V2.^2);
az      = atan2d(V1,V2);

% % (optional) remove anywhere where zenith angle > 90 (below the horizon)
% zen(zen >= 90) = 90;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. convert zenith angle values of the airglow layer at height h above
% the ground to arc length along that airglow layer (see diagram for
% method)
theta = zen;

% first find the angle alpha: (use projection to PQ axis method)
alpha = 90 - theta;

% find the angle phi
phi = acosd(Re.*sind(theta)./r);

% now find the angle beta (the arc angle):
beta = phi - alpha;
arclen = beta; % in DEGREES

% convert to arcdist and project into x,y distance:
arcdist = (2.*pi.*r) .* (arclen ./ 360);
x_in       = arcdist.*sind(az);
y_in       = arcdist.*cosd(az);

% great! we've now got the X,Y distance locations of every point on the image.
% Now let's put this on a regular distance grid centre on the camera.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Create a new regular distance grid and project onto this

D1 = newgridX; D2 = newgridY;

% interpolate on to this new grid:
F = scatteredInterpolant(x_in(:),y_in(:),Z(:),'linear','none');
regridded_image = F(D1,D2);
% there is probably a better way to do this interpolation step but this is
% simple and effective

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Find new lat lons of this new regular grid

% this is a little sticky - easiest way is to use reckon, but with an
% increased earth radius for the airglow layer
RS = referenceSphere('earth');
RS.Radius = RS.Radius + h;

[Lat,Lon] = reckon(cameralat,cameralon,km2deg(sqrt(D1.^2+D2.^2)),atan2d(D1,D2));

end







