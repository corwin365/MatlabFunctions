function POut = h2p(alt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert altitudes to pressures, with the U.S. Standard Atmosphere, 1976
%
%inputs:
%   alt - geometric altitude, in km (any size of array is fine)
%outputs:
%   POut - atmospheric pressure, in hPa, same size, double. NaN outside -5 to 149 km, or for NaN input.
%
%Up to 86 km: the standard's defining equations, eqs (18), (23), (33a), (33b) (pp. 8-12); geometric altitude is
%converted to geopotential, and the pressure follows from the seven layers of molecular-scale temperature
%(constants in ussa76_constants.m). Above 86 km, where the standard's pressure comes from a diffusion model rather
%than a formula: linear interpolation in log(pressure) between the values of its Table I (p. 68), every 0.5 km to
%100 km and every 1 km to 149 km. Agrees with the printed Table I to its 5 digits.
%
%Source: U.S. Standard Atmosphere, 1976, NOAA-S/T 76-1562 / NASA-TM-X-74335, Washington, D.C., 1976.
%Python equivalent: atmosgw.curate.atmosphere.height_to_pressure.
%
%Corwin Wright, c.wright@bath.ac.uk, 2025/04/02 (full rewrite of older 2023 version to properly vectorise,
%which was itself based on much older code by Neil)
%
%change log:
%  2026/11/25 [Claude] rewritten to implement the standard properly (was: isothermal layers between 7 pressures,
%                      with geopotential base heights used as geometric; NaN below 0 km)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = ussa76_constants();
z = double(alt);
POut = NaN(size(z));

%geopotential altitude, eq (18)
H = S.r0.*z./(S.r0 + z);

%up to 86 km (84.852 km'): equations
Low = isfinite(H) & z >= S.MinAlt & H <= S.Hb(end);
h = H(Low); h = h(:);
b = max(sum(h >= S.Hb(1:end-1),2),1);   %layer of each point
Pb = S.Pb(b)'; Tb = S.Tb(b)'; L = S.L(b)'; dH = h - S.Hb(b)';
P = Pb.*(Tb./(Tb + L.*dH)).^(S.k./L);  %eq (33a)
Iso = L == 0;
P(Iso) = Pb(Iso).*exp(-S.k.*dH(Iso)./Tb(Iso));  %eq (33b)
POut(Low) = P;

%86 to 149 km: Table I, log(pressure) linear in height. The equations end at 84.852 km' = 85.99995 km geometric
%and the table starts at 86 km: in the 5 cm between, the table's 86 km value is used
High = isfinite(z) & ~Low & H > S.Hb(end) & z <= S.TableZ(end);
POut(High) = exp(interp1(S.TableZ,log(S.TableP),max(z(High),S.TableZ(1)),'linear'));

end
