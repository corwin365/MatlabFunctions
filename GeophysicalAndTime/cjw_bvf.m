function [BVF] = cjw_bvf(Height,Theta,SmoothLevels)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute BVF from temp, altitude and theta
%
%reimplementation of equivalent IDL routine bvf.pro, which was
%itself adapted from an original piece of code by Tanya Peevey
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%09/JAN/2014
%
%inputs
%---------
%
%Height - altitude profile
%Theta - potential temperature profile
%Smooth - number of levels to boxcar-smooth by (default: no smoothing)
%
%outputs
%---------
%
%BVF - brunt-vaisala profile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BVF = zeros(numel(Height),1);
for i=1:1:numel(Height)-1
%   Grad_Freq=(9.81/Theta[*,i])*[(Theta[*,i+1]-Theta[*,i])/((Alti[*,i+1]-Alti[*,i]))]
  BVF(i) = sqrt(abs((9.81/Theta(i)).* (Theta(i+1)-Theta(i))./(Height(i+1)-Height(i))));
end

%apply smoothing, if requested
if nargin == 3; BVF = smooth(BVF,SmoothLevels);end
