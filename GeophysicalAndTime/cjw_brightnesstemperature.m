function BT = cjw_brightnesstemperature(Wavenumber,Radiance)

%frequency in cm^-1
%radiance in mW/(m2.sr.cm-1)

%adapted from https://svn.ssec.wisc.edu/repos/bennartz_group/LIBRARY/idl/std_libs/CRTM_PLANCK/planck_temperature.pro
%converted to Matlab CJW 01/01/2014 (omits most of error checking, etc, to speed up translation)

%  For wavenumber input, the Planck Temperature is calculated using:
% ;
% ;                    C2 * wavenumber
% ;         T = -------------------------------
% ;                ( C1 * wavenumber^3      )
% ;              LN( ----------------- +  1 )
% ;                (          B             )
% ;
% ;       For Wavelength input :
% ;
% ;                                C2
% ;          T = -----------------------------------------
% ;                              (        C1            )
% ;               Wavelength * LN( ---------------- + 1 )
% ;                              ( Wavelength^5 * B     )
% ;
% ;       C1 and C2 are determined using:
% ;          C1 = 2*h*c*c  [W.m2]
% ;          C2 = h*c/k    [K.m]
% ;
% ;       The fundamental and derived constants used are defined in the
% ;       include file fundamental_constants.pro
% ;
% ;       A scaling is applied to the Planck Radiance results to return
% ;       Radiances in the units of mW/(m2.sr.cm-1) for wavenumber input
% ;       or W/(m2.sr.micron) for wavelength input.


%define constants
h = 6.625e-34; %Planck constant
k = 1.381e-23; %Boltzmann constant
c = 2.997e+08; %speed of light

%convert to base units
%frequency in cm^-1
%radiance in mW/(m2.sr.cm-1)
Wavenumber = Wavenumber .*  100; %cm^-1 -> m^-1
Radiance   = Radiance   ./ 1000; %mW -> W
Radiance   = Radiance   ./  100; %/cm^-1 -> /m^-1

%calculation
C1 = 2 .* h .* c .* c;
C2 = h .* c ./ k;

Numerator = C2 .* Wavenumber;
Denominator = 1+((C1 .* Wavenumber .* Wavenumber .* Wavenumber) ./ Radiance);

BT = Numerator ./ log(Denominator);

