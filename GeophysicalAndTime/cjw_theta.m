function Theta = cjw_theta(Pressure,Temperature,SurfacePressure)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute theta from pressure and temperature
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%09/JAN/2014
%
%inputs
%---------
%
%Temperature - temperature  (K)
%Pressure - pressure (hPa)
%SurfacePressure - (optional) surface pressure (hPa)
%
%outputs
%---------
%
%Theta - potential temperature
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%surface pressure supplied?
if nargin<3; SurfacePressure = 1000.; end

%scale pressure
Pressure = SurfacePressure./Pressure;

%compute theta
Theta =  Temperature .* Pressure.^(0.28571429);

%done!

