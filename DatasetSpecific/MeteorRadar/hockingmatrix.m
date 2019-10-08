function HockingMatrix = hockingmatrix(theta,phi,ModRadialVelocity)
                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute Hocking (2005)'s left-hand matrix for a given set of parameters
%see page 2434 of Hocking et al (2005)
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%15/JAN/2014
%
%derived from original code sourced from Robin Davis
%
%inputs
%---------
%
%theta             -  zenith of the meteor (angle from vertical, 0 straight up)
%phi               - azimuth of the meteor (angle anticlockwise from due east)
%ModRadialVelocity - radial PERTURBATION velocity
%
%
%outputs
%---------
%
%Hocking matrix - results in order specified by matrix x of Hocking et al (2005)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tidying - some NaNs to be removed from not enough meteors in the windfitting:

% the matrix from Hocking 2005 - AM: Checked Twice - Most Recent Check
% 12-07-2014
st = sin(theta);
ct = cos(theta);
sp = sin(phi);
cp = cos(phi);
st2 = (sin(theta)).^2;
ct2 = (cos(theta)).^2;
sp2 = (sin(phi)).^2;
cp2 = (cos(phi)).^2;
st3 = (sin(theta)).^3;
ct3 = (cos(theta)).^3;
sp3 = (sin(phi)).^3;
cp3 = (cos(phi)).^3;
st4 = (sin(theta)).^4;
ct4 = (cos(theta)).^4;
sp4 = (sin(phi)).^4;
cp4 = (cos(phi)).^4;
v2  = ModRadialVelocity.^2;
L1 = [sum(st4.*cp4) sum(st4.*cp2.*sp2) sum(st2.*ct2.*cp2) sum(2*st4.*cp3.*sp) sum(2*st3.*ct.*cp3) sum(2*st3.*ct.*cp2.*sp)];
L2 = [sum(st4.*cp2.*sp2) sum(st4.*sp4) sum(st2.*ct2.*sp2) sum(2*st4.*sp3.*cp) sum(2*st3.*ct.*sp2.*cp) sum(2*st3.*ct.*sp3)];
L3 = [sum(st2.*ct2.*cp2) sum(st2.*ct2.*sp2) sum(ct4) sum(2*st2.*ct2.*cp.*sp) sum(2*ct3.*st.*cp) sum(2*ct3.*st.*sp)];
L4 = [sum(2*st4.*cp3.*sp) sum(2*st4.*sp3.*cp) sum(2*st2.*ct2.*cp.*sp) sum(4*st4.*cp2.*sp2) sum(4*st3.*ct.*cp2.*sp) sum(4*st3.*ct.*sp2.*cp)];
L5 = [sum(2*st3.*ct.*cp3) sum(2*st3.*ct.*sp2.*cp) sum(2*ct3.*st.*cp) sum(4*st3.*ct.*cp2.*sp) sum(4*st2.*ct2.*cp2) sum(4*st2.*ct2.*cp.*sp)];
L6 = [sum(2*st3.*ct.*cp2.*sp) sum(2*st3.*ct.*sp3) sum(2*ct3.*st.*sp) sum(4*st3.*ct.*sp2.*cp) sum(4*st2.*ct2.*cp.*sp) sum(4*st2.*ct2.*sp2)];

%gluing the matrix back together - Hocking's A-matrix (the big one on p2434)
A = [L1;L2;L3;L4;L5;L6]; % A == A'

%right hand side of Hocking's page 2434 equation
b = [sum(v2.*st2.*cp2);sum(v2.*st2.*sp2);sum(v2.*ct2);sum(2*v2.*st2.*cp.*sp);sum(2*v2.*st.*ct.*cp);sum(2*v2.*st.*ct.*sp)];

%result we want! (and some other stuff...)
% warning off % matrix may be singular, due to sparse data, so temporarily disable warnings
HockingMatrix = A\b; %\ -> matlab's matrix division using Gaussian elimination (?): b matrix divided by A
                       % AM 12-07-2014 - Equivalent to b/A (b would need to
                       % be changed to a row from column for operation to
                       % work here).
% warning on
end