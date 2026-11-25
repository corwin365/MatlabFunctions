function z = p2h(pres)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert pressures to altitudes, with the U.S. Standard Atmosphere, 1976 (the inverse of h2p)
%
%inputs:
%   pres - atmospheric pressure, in hPa (any size of array is fine)
%outputs:
%   z - geometric altitude, in km, same size, double. NaN outside 1777.6 to 4.74e-6 hPa (-5 to 149 km), or for
%       NaN or non-positive input.
%
%Up to 86 km: the inverses of the standard's eqs (33a) and (33b) for geopotential height in each layer, then
%eq (19) for geometric height (constants in ussa76_constants.m). Above 86 km: height linear in log(pressure)
%between the values of its Table I (p. 68). See h2p for details.
%
%Source: U.S. Standard Atmosphere, 1976, NOAA-S/T 76-1562 / NASA-TM-X-74335, Washington, D.C., 1976.
%Python equivalent: atmosgw.curate.atmosphere.pressure_to_height.
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/08/14 (originally a tidied-up pres2alt_complex() by Neil Hindley)
%
%change log:
%  2026/11/25 [Claude] rewritten to implement the standard properly (was: isothermal layers between 7 pressures,
%                      with geopotential base heights used as geometric; pressures above 1013.25 hPa reused the
%                      previous element's layer; single input computed in single precision)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = ussa76_constants();
p = double(pres);
z = NaN(size(p));

PBottom = h2p(S.MinAlt);
PTop86  = S.Pb(end).*(S.Tb(end)./(S.Tb(end) + S.L(end).*(S.Hb(end)-S.Hb(end-1)))).^(S.k./S.L(end));  %at 84.852 km'

%up to 86 km: inverted equations
Low = isfinite(p) & p > 0 & p <= PBottom & p >= PTop86;
q = p(Low); q = q(:);
b = max(sum(q <= S.Pb,2),1);   %layer: highest base whose pressure is >= q
Pb = S.Pb(b)'; Tb = S.Tb(b)'; L = S.L(b)';
dH = Tb./L.*((q./Pb).^(-L./S.k) - 1);   %inverse of eq (33a)
Iso = L == 0;
dH(Iso) = -Tb(Iso)./S.k.*log(q(Iso)./Pb(Iso));   %inverse of eq (33b)
H = S.Hb(b)' + dH;
z(Low) = S.r0.*H./(S.r0 - H);   %eq (19)

%86 to 149 km: height linear in log(pressure) between Table I values. Pressures between the table's 86 km value
%and the equations' (1.2e-6 relative apart) are set to 86 km
High = isfinite(p) & p < PTop86 & p >= S.TableP(end);
z(High) = interp1(log(S.TableP),S.TableZ,min(log(p(High)),log(S.TableP(1))),'linear');

end
