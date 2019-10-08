function [OmegaHat,GroupVelocity,PhaseVelocity,f] = intrinsic_frequency(N,k,l,m,Latitude)
  %calculates GW intrinsic frequency using eqn 23 of Fritts + Alexander 2003
  
  %also calculates "intrinsic" group velocity, since we're here anyway
  %this needs (u,v,w) adding to it to make a true group velocity
  
  %ok then. intrinsic phase speed too, if you insist.

  %latitude must be specified in *degrees*. 
  %k,l,m are in units of 1/f, *not* 2pi/f
  %Everything else in base units.
  
  
  %here's some testing data
% %   clear all
% %   Latitude = -55;
% %   N        = 0.02;
  
% %   k = -1./(435e3);
% %   l = 1./(20000e3);
% %   m = 1./(22e3);
  
% %   k = 1./(200e3);
% %   l = 1./(200e3);
% %   m = 1./(10e3);
  


  %rescale k,l,m from a pi/f to a 2pi/f space
  %this is because the ST spits out 1/f and I use it everywhere else
  %this function has to be in the right units for BigOmega to be correct
  k = k*2*pi;
  l = l*2*pi;
  m = m*2*pi;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %intrinsic frequency
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %find the Earth's rotation rate
  BigOmega = 2.*pi./(24.*60.*60);
  %and hence the Coriolis parameter
  f = abs(2.*BigOmega.*sind(Latitude));
  clear BigOmega Latitude
  
  %set the scale height to 7km
  H = 7000; Hfrac = 1/(4*(H^2));
  
  %hence, find OmegaHat.^2
  TopTerm1 = (N.^2).*(k.^2 + l.^2);
  TopTerm2 = (f.^2).*(m.^2 + Hfrac);
  BottomTerm = k.^2 + l.^2 + m.^2 + Hfrac;
  OmegaHatSquared = (TopTerm1+TopTerm2)./BottomTerm;
  clear TopTerm1 TopTerm2 BottomTerm
  
  %and OmegaHat
  OmegaHat = sqrt(OmegaHatSquared);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %(intrinsic) group velocity
  %dW/dk
  %equation 25 of Fritts + Alexander 2003
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  TopTermK =  k.*(N.^2 -OmegaHat.^2);
  TopTermL =  l.*(N.^2 -OmegaHat.^2);
  TopTermM = -m.*(OmegaHat.^2-f.^2);
  BottomTerm = OmegaHat.*(k.^2 + l.^2 + m.^2 +Hfrac);
  
  GroupVelocity.k = TopTermK./BottomTerm;
  GroupVelocity.l = TopTermL./BottomTerm;
  GroupVelocity.m = TopTermM./BottomTerm;
  clear TopTermK TopTermL TopTermM BottomTerm

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %(intrinsic) phase velocity
  %W/k
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  PhaseVelocity.k = OmegaHat./k;
  PhaseVelocity.l = OmegaHat./l;
  PhaseVelocity.m = OmegaHat./m;
  
  %we have an occasional problem for waves aligned very close to
  %EW or NS where 1/k or 1/l are very large and hence the phase 
  %velocity becomes excessively large. identify these and set 
  %them to NaN
  
  %two tests:
    %1. large wavevector magnitude
    %2. large angle between k and l
  %must pass both of these to be NaNned 
  
  Tests = [0,0];
  
  %test 1. magnitude. this also identifies which one is the scary one
  CutOffLambda = 2500e3; %2500km is quite big, compare AIRS swath ~1800km.
  if abs(2*pi/k) > CutOffLambda; Tests(1) = Tests(1)+1; end
  if abs(2*pi/l) > CutOffLambda; Tests(2) = Tests(2)+1; end
  clear CutOffLambda
  
  
  if sum(Tests) > 0
    %test 2: angle
    CutOffRange = [85,95]; % degrees between vectors
    
    Angle = atan2d(k,l);
    if abs(Angle) > CutOffRange(1) & abs(Angle) < CutOffRange(2);
      %test 2 failed. set the wavelength which failed test 1 to NaN
      if Tests(1) == 1; PhaseVelocity.k = NaN; end
      if Tests(2) == 1; PhaseVelocity.l = NaN; end      
    end
    clear CutOffRange
  end

  clearvars -except PhaseVelocity GroupVelocity OmegaHat f
  return


