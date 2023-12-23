function Results = wave_properties(A,k,l,m,varargin)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the following gravity wave properties:
%
%  potential energy                             (J.kg^-1)
%  signed momentum flux                         (mPa)    
%  intrinsic and ground-based frequency         (s^-1)  
%  "    "    "    "    "    " phase speed       (m.s^-1)
%  "    "    "    "    "    " group velocity    (m.s^-1)
%  wave horizontal phase propagation bearing    (deg c/w/N) <- not currently working, fix maths
%  wave horizontal group propagation bearing    (deg c/w/N) <- not currently working, fix maths
%
%replaces a 2017 routine of the same name, but NOT BACKWARDS COMPATIBLE
%
%some default choices I have made may not make sense for most applications, but
%instead reflect my local code defaults. If using this routine for the first
%time, I **strongly** suggest setting the 'Verbose' flag to true to see what's
%being assumed.
%
%'FA2003' in comments is Fritts and Alexander (Rev. Geophys 2003, doi:10.1029/2001RG000106)
%
%Corwin Wright, c.wright@bath.ac.uk
%2023/12/23
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs: (* indicates REQUIRED, all others are optional)
%
%     VarName   (type,  default)  description
%  -------------------------------------------------------------------------------- 
%  *  A         (numeric       )  wave amplitude (K)
%  *  k         (numeric       )  zonal      wavenumber (units depend on TwoPiData and MetreData flags)
%  *  l         (numeric       )  meridional wavenumber "     "     "     "     "     "     "     "   
%  *  m         (numeric       )  vertical   wavenumber "     "     "     "     "     "     "     " 
%
%     Verbose   (logical, false)  flag: should we tell the user about implemented choices?
%     TwoPiData (logical, false)  flag: are [k,l,m] defined as 2pi/f (true) or   1/f (false)?
%     MetreData (logical, false)  flag: are [k,l,m] in units of m^-1 (true) or km^-1 (false)?
%     ForceUp   (logical,  true)  flag: should we assume the waves are  travelling vertical upwards, i.e. -ve m?
%
%     Lat       (numeric,    60)  latitude                    (degrees)
%     g         (numeric,  9.81)  acceleration due to gravity (m.s^-1)
%     N         (numeric,  0.02)  Brunt-Vaisala frequency     (s^-1)
%     H         (numeric,  7000)  atmospheric scale height    (m)
%     T         (numeric,   250)  background temperature      (K)
%     Dens      (numeric,   NaN)  background air density      (kg.m^-3)
%     P         (numeric,     3)  background pressure         (hPa) - only used if Dens is NaN
%     U         (numeric,   NaN)  background zonal      wind  (m.s^-1)
%     V         (numeric,   NaN)  background meridional wind  (m.s^-1)
%     W         (numeric,     0)  background vertical   wind  (m.s^-1)
%
%Most inputs (all except the logical flags) can be supplied as a multidimensional array or a single value. If
%arrays are used then all values must be consistently either the same size OR a single value.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs: (all as fields of the struct 'Results', and all the same size/type as the input data)
%
%basic properties:
%   Ep         - potential energy (J.kg^-1)
%   Mz         - zonal      momentum flux (mPa)
%   Mm         - meridional momentum flux (mPa)
%   Omega_I    - intrinsic    frequency (s^-1)
%   Omega_G    - ground-based frequency (s^-1)
%
%group speeds:
%   Group_I_Zonal      - intrinsic    zonal
%   Group_I_Merid      - intrinsic    meridional
%   Group_I_Horizontal - intrinsic    horizontal (i.e. sqrt(zonal^2 + merid^2))
%   Group_I_Vertical   - intrinsic    vertical
%   Group_G_Zonal      - ground-based zonal
%   Group_G_Merid      - ground-based meridional
%   Group_G_Horizontal - ground-based horizontal (i.e. sqrt(zonal^2 + merid^2))
%   Group_G_Vertical   - ground-based vertical
%
%phase speeds:
%  Phase_I_Vertical   - intrinsic    vertical
%  Phase_I_Horizontal - intrinsic    horizontal
%  Phase_G_Vertical   - ground-based vertical
%  Phase_G_Horizontal - ground-based horizontal
%
%directions:
%  BearingGroup        - group velocity direction in horizontal plane <-- currently not outputted as need to fix a maths mistake
%  BearingPhase        - phase velocity direction in horizontal plane <-- currently not outputted as need to fix a maths mistake
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%validation of inputs
%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

%required variables:
addRequired(p,'A',@isreal);  %amplitude
addRequired(p,'k',@isreal);  %zonal wavenumber
addRequired(p,'l',@isreal);  %meridional wavenumber
addRequired(p,'m',@isreal);  %vertical wavenumber

%optional variables and default values
addParameter(p,'Lat',  60,@(x) validateattributes(x,{'numeric'},{'>=',-90,'<=',90}))
addParameter(p,'N',  0.02,@isreal) %i.e. N = 0.02 s^-1
addParameter(p,'g',  9.81,@isreal) %i.e. g = 9.81 ms^-1
addParameter(p,'H',  7000,@isreal) %i.e. H = 7000 m
addParameter(p,'U',   NaN,@isreal) 
addParameter(p,'V',   NaN,@isreal) 
addParameter(p,'W',     0,@isreal) %one can dream
addParameter(p,'Dens',NaN,@isreal) %one can dream
addParameter(p,'T',   250,@isreal) %typical temperature at ~40km, the best level of AIRS data
addParameter(p,'P',     3,@isreal) %typical pressure    at ~40km, the best level of AIRS data


%flags
addParameter(p,'TwoPiData',false,@islogical); %do we want to convert the data from 1/f space to 2pi/f space?
addParameter(p,'MetreData',false,@islogical); %is the k,l,m data in 1/metres or 1/kilometres?
addParameter(p,'ForceUp',  true, @islogical); %should we force the waves to be vertically ascending, i.e. -ve m?
addParameter(p,'Verbose', false, @islogical); %tell the user about choices we have made

%OK, we're done. Parse the input tree, then rename the variables to the format used in the rest of the routine
parse(p,A,k,l,m,varargin{:})
Settings = p.Results;
clearvars -except Settings A k l m 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preliminary computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tell the user what's happening
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Verbose == 1;
  disp('==========================')
  disp('Computing GW properties')
  disp('++++++++++++++++++++++++++')
  disp(' ')
  disp('Assumptions (set optional flags to override):')
end

%copy working vars to Results, to make easier to store/use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results.A = A;
Results.k = k;
Results.l = l;
Results.m = m;

%force vertical ascent?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.ForceUp == 1;
  Positive = find(m > 0);
  GWs.k(Positive) = -k(Positive);
  GWs.l(Positive) = -l(Positive);
  GWs.m(Positive) = -m(Positive);
  clear Positive
  if Settings.Verbose == 1; disp('--> Waves are propagating vertically upwards, i.e. m < 0'); end
end

%1/f -> 2pi/f conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%our calculations are carried out in 2pi/f space, but my usual input data is in 1/f space
%accordingly, convert the data to 1/f space unless told otherwise
if Settings.TwoPiData == 0
  k = 2*pi.*k;
  l = 2*pi.*l;
  m = 2*pi.*m;
  if Settings.Verbose == 1; disp('--> Wavenumbers are defined as 1/f, not 2pi/f'); end  
end; Settings = rmfield(Settings,'TwoPiData'); 

%km -> m conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.MetreData == 0
  k = k./1e3;
  l = l./1e3;
  m = m./1e3;  
  if Settings.Verbose == 1; disp('--> Wavenumbers are in units of 1/km, not 1/m'); end   
end; Settings = rmfield(Settings,'MetreData'); 

if Settings.Verbose == 1;
  disp(' ')
  disp('Setting parameters (set optional variables to override):')
end

%atmospheric density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isnan(Settings.Dens);
  Dens = cjw_airdensity(Settings.P,Settings.T);
  if Settings.Verbose == 1;
    disp('--> Computing background density based on P and T'); 
    if Settings.P ==   3;disp('----> P may not have been specified, using 3hPa (typical ~40km value)');  end
    if Settings.T == 250;disp('----> T may not have been specified, using 250K (typical ~40km value)');  end
  end
else Dens = Settings.Dens; end
Settings = rmfield(Settings,'Dens');

%Coriolis parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find the Earth's rotation rate and hence the Coriolis parameter
BigOmega = 2.*pi./(24.*60.*60);
f = 2.*BigOmega.*sind(Settings.Lat);
clear BigOmega Lat

%U and V status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Verbose == 1
  if isnan(Settings.U); disp('--> U unspecified, zonal      ground-based velocities will be all-NaN'); end
  if isnan(Settings.V); disp('--> V unspecified, meridional ground-based velocities will be all-NaN'); end
  if Settings.W == 0;   disp('--> W unspecified, assuming zero'); end
end

%compute 1/4H^2 and N^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we use these often, so give them simple names
Hfrac = 1./((2.*Settings.H).^2);
N2    = (Settings.N).^2;

%wave propagation direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Theta = atan2d(-l,-k);

%tell the user what's happening
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.Verbose == 1; disp(' ');  disp('Configuration complete, processing data'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% potential energy and momentum flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%computation
PotentialEnergy = 0.5 .* (Settings.g.^2 ./ N2) .* (A./Settings.T).^2;
MomentumFlux.x = PotentialEnergy .* -Dens .* k./m .*1000;
MomentumFlux.y = PotentialEnergy .* -Dens .* l./m .*1000;

%copy to output struct
Results.Ep = PotentialEnergy;
Results.Mz = MomentumFlux.x;
Results.Mm = MomentumFlux.y;

clear Dens PotentialEnergy MomentumFlux A

if Settings.Verbose == 1; disp('--> Ep and MF computed'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wave intrinsic frequency - intrinsic and ground-based
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this is equation 23 of FA2003
TopTerm1        = N2.*(k.^2 + l.^2);
TopTerm2        = (f.^2).*(m.^2 + Hfrac);
BottomTerm      = k.^2 + l.^2 + m.^2 + Hfrac;
OmegaHatSquared = (TopTerm1+TopTerm2)./BottomTerm;
OmegaHat        = sqrt(OmegaHatSquared);
clear TopTerm1 TopTerm2 BottomTerm H OmegaHatSquared

%copy to output struct
Results.Omega_I = OmegaHat;
Results.Omega_G = OmegaHat + k .* Settings.U.*cosd(Theta) + l .* Settings.V.*sind(Theta);

if Settings.Verbose == 1; disp('--> Frequency computed'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% group velocity - intrinsic 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute in absolute space and then project to directions, to reduce risk of singularities
kh = quadadd(k,l);

%use the full expression, equation 25 of FA2003
%testing this against the simpler midfrequency version for a few waves gives very similar results
TopTermKL =  kh.*(N2 -OmegaHat.^2);
TopTermM =  -m.*(OmegaHat.^2-f.^2);
BottomTerm = OmegaHat.*(kh.^2 + kh.^2 + m.^2 +Hfrac);
Intrinsiccg = TopTermKL./BottomTerm;

%and project
Theta = atan2d(l,k);
Group.kI = Intrinsiccg .* cosd(Theta);
Group.lI = Intrinsiccg .* sind(Theta);

%vertical handled separately
Group.mI = TopTermM ./ BottomTerm;

%rename for output
Results.Group_I_Zonal      = abs(Group.kI);
Results.Group_I_Merid      = abs(Group.lI);
Results.Group_I_Horizontal = quadadd(Group.kI,Group.lI);
Results.Group_I_Vertical   = abs(Group.mI);


clear TopTermM TopTermKL BottomTerm Intrinsiccg 
clear Hfrac N2 f OmegaHat %we're done with these

if Settings.Verbose == 1; disp('--> Group speeds computed'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% group velocity - ground-based
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute ground-based group speed
Group.kG = Group.kI + Settings.U.*cosd(Theta);
Group.lG = Group.lI + Settings.V.*sind(Theta);
Group.mG = Group.mI + Settings.W;

%rename for output
Results.Group_G_Zonal      = abs(Group.kG);
Results.Group_G_Merid      = abs(Group.lG);
Results.Group_G_Horizontal = quadadd(Group.kG,Group.lG);
Results.Group_G_Vertical   = abs(Group.mG);

clear Group

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% phase speed - intrinsic and ground-based
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%phase speed is defined as frequency / wavenumber
%it is a SPEED, not a vector, and is defined in the direction [k,l,m]

%vertical
Results.Phase_I_Vertical = Results.Omega_I./m;
Results.Phase_G_Vertical = Results.Omega_G./m;

%horizontal
Results.Phase_I_Horizontal = Results.Omega_I./kh;
Results.Phase_G_Horizontal = Results.Omega_G./kh;

if Settings.Verbose == 1; disp('--> Phase speeds computed'); end

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %% wave direction in horizontal plane
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % Results.BearingGroup = wrapTo180(  Theta+90);
% % % Results.BearingPhase = wrapTo180(-(Theta+90));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% done!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%finished!
if Settings.Verbose == 1;
  disp(' ')
  disp('--------------------------')
  disp('GW properties computed')
  disp('==========================')
end


return