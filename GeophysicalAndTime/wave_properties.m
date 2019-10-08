function [Intrinsic,True,f] = wave_properties(N,k,l,m,Latitude,U,V)
  %calculate a range of GW properties, as intrinsic (i.e. wind-relative)
  %and true (i.e. ground-relative) values

  %latitude must be specified in *degrees*. 
  %k,l,m are in units of 1/f, *not* 2pi/f
  %Everything else in base units.
  
  %U and V are optional, but true values won't be computed without them
  
  %modified 2017/SEP/03 to fix incorrect sign in ground-based freq calculation
  %modified 2018/APR/07 to tidy up phase speed calcn

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %setup
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %rescale k,l,m from a 1/f to a 2pi/f space
  %this is because the ST spits out 1/f and I use it everywhere else
  %this function has to be in the right units for BigOmega to be correct
  k = k.*2.*pi;
  l = l.*2.*pi;
  m = m.*2.*pi;
  %kh is computed in the routine from k and l

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %intrinsic frequency
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %find the Earth's rotation rate
  BigOmega = 2.*pi./(24.*60.*60);
  %and hence the Coriolis parameter
  f = 2.*BigOmega.*sind(Latitude);
  clear BigOmega Latitude
  
  %set the scale height to 7km
  H = 7000; Hfrac = 1/(4*(H^2));
  
  %hence, find OmegaHat.^2
  TopTerm1 = (N.^2).*(k.^2 + l.^2);
  TopTerm2 = (f.^2).*(m.^2 + Hfrac);
  BottomTerm = k.^2 + l.^2 + m.^2 + Hfrac;
  OmegaHatSquared = (TopTerm1+TopTerm2)./BottomTerm;
  clear TopTerm1 TopTerm2 BottomTerm H
  
  %and OmegaHat
  OmegaHat = sqrt(OmegaHatSquared);
  Intrinsic.OmegaHat = OmegaHat;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %(intrinsic) group velocity
  %dW/dk
  %equation 25 of Fritts + Alexander 2003
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  TopTermK =  k.*(N.^2 -OmegaHat.^2);
  TopTermL =  l.*(N.^2 -OmegaHat.^2);
  TopTermM = -m.*(OmegaHat.^2-f.^2);
  BottomTerm = OmegaHat.*(k.^2 + l.^2 + m.^2 +Hfrac);
  
  Group.k = TopTermK./BottomTerm;
  Group.l = TopTermL./BottomTerm;
  Group.m = TopTermM./BottomTerm;
  clear TopTermK TopTermL TopTermM BottomTerm N
  
  %horizontal absolute and direction
  Group.mag   = sqrt(Group.k.^2 + Group.l.^2);
  Group.theta = atan2(Group.l,Group.k);

  Intrinsic.Group = Group; 
  clear Group
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %(intrinsic) phase velocity
  %this is very sensitive to spikes caused by waves aligned 
  %near the cardinal directions. so, compute in r,theta terms
  %for horizontal (vertical is fine)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  %work out leading term
  LeadTerm = OmegaHat./(k.^2 + l.^2 + m.^2);
  
  %and hence the individual projections
  Phase.k = LeadTerm.*k;
  Phase.l = LeadTerm.*l;
  Phase.m = LeadTerm.*m;
  
  %and magnitude and horizontal direction
  Phase.mag   = sqrt(Phase.k .^2 + Phase.l .^2);
  
  %direction isn't defined for phase speed.
  
  %first, vertical
  Phase.m = OmegaHat./m;
    
  %work out horizontal wavelength and angle
  kh = sqrt(k.^2 + l.^2);
  theta = atan2(l,k);
  
  %hence, work out phase *speed*
  Phase.kh = OmegaHat./kh;
  Phase.theta = theta;
  
  %and projections
  Phase.k = Phase.kh.*cos(theta);
  Phase.l = Phase.kh.*sin(theta);

  Intrinsic.Phase = Phase;
  clear Phase kh OmegaHat Hfrac OmegaHatSquared theta 

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %OK. now we need to convert everything from intrinsic values
  %to true ground-based values
  %using wind speed from elsewhere (tested with ECMWF)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  if nargin > 5
    
    %project wind into direction of wave, for kh values
    WindGroupDir = U./cos(Intrinsic.Group.theta) + V./sin(Intrinsic.Group.theta);
    WindPhaseDir = U./cos(Intrinsic.Phase.theta) + V./sin(Intrinsic.Phase.theta);
    %remove outliers arising due to projecting to large angles
    WindGroupDir(abs(WindGroupDir) > 200) = NaN;
    WindPhaseDir(abs(WindPhaseDir) > 200) = NaN;
    
    %omegahat
    True.Omega    = Intrinsic.OmegaHat + k.*U + l.*V;
    
    %group velocity
    True.Group.k  = Intrinsic.Group.k + U;
    True.Group.l  = Intrinsic.Group.l + V;
    True.Group.m  = Intrinsic.Group.m; %assumes W=0, usually fair abov tsphere
    True.Group.kh = abs(Intrinsic.Group.mag+WindGroupDir);
    True.Group.th = atan2(True.Group.l,True.Group.l);
    
    %phase velocity
    True.Phase.k  = Intrinsic.Phase.k + U;
    True.Phase.l  = Intrinsic.Phase.l + V;
    True.Phase.m  = Intrinsic.Phase.m; %viz
    True.Phase.kh = abs(Intrinsic.Phase.mag+WindPhaseDir);
    True.Phase.th = atan2(True.Phase.l,True.Phase.l);
  else
    True = [];
  end

