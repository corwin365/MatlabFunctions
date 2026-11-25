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
%eq (19) for geometric height. Above 86 km: height linear in log(pressure)
%between the values of its Table I (p. 68). See h2p for details.
%
%Self-contained: the constants below are duplicated in h2p.m; keep the two copies identical.
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
%  2026/11/25 [Claude] constants included here (were in ussa76_constants.m), so the function needs no other file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = constants();
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

function S = constants()
%constants of the U.S. Standard Atmosphere, 1976 (identical copy in h2p.m and p2h.m)
%   g0, Rs, M0, r0, T0, P0 - g0' (m2 s-2 m'-1), R* (N m kmol-1 K-1), M0 (kg kmol-1), r0 (km), T0 (K), P0 (hPa);
%                            table 2 (p. 2) and p. 9
%   Hb, L                  - layer base geopotential heights (km', 8 values) and molecular-scale temperature
%                            gradients (K/km', 7 values); table 4 (p. 3)
%   Tb, Pb                 - base temperatures (K) and pressures (hPa) of the 7 layers, from eqs (23), (33a), (33b)
%   k                      - g0' M0 / R*, in K per km'
%   MinAlt                 - lowest valid geometric altitude, km (start of Table I, p. 51)
%   TableZ, TableP         - Table I (geometric altitude, p. 68): Z (km) and P (hPa), 86-149 km. Each value was
%                            checked against the table's independent P/P0 column.

S.g0 = 9.80665;      %g0', m2 s-2 m'-1
S.Rs = 8.31432e3;    %R*, N m kmol-1 K-1
S.M0 = 28.9644;      %M0, kg kmol-1
S.r0 = 6356.766;     %r0, km
S.T0 = 288.15;       %K
S.P0 = 1013.25;      %hPa

S.Hb = [0, 11, 20, 32, 47, 51, 71, 84.8520];      %km'
S.L  = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0];    %K/km'
S.k  = S.g0.*S.M0./S.Rs.*1e3;                     %K/km'

%base temperatures and pressures of each layer, eqs (23), (33a), (33b)
S.Tb = S.T0; S.Pb = S.P0;
for b=1:1:numel(S.L)-1
  dH = S.Hb(b+1) - S.Hb(b);
  S.Tb(b+1) = S.Tb(b) + S.L(b).*dH;
  if S.L(b) == 0; S.Pb(b+1) = S.Pb(b).*exp(-S.k.*dH./S.Tb(b));
  else;           S.Pb(b+1) = S.Pb(b).*(S.Tb(b)./(S.Tb(b) + S.L(b).*dH)).^(S.k./S.L(b));
  end
end

S.MinAlt = -5;   %km geometric

S.TableZ = [86:0.5:99.5, 100:1:149];
S.TableP = [3.7338e-3, 3.4163e-3, 3.1259e-3, 2.8602e-3, 2.6173e-3, 2.3951e-3, ...
            2.1919e-3, 2.0060e-3, 1.8359e-3, 1.6804e-3, 1.5381e-3, 1.4078e-3, ...
            1.2887e-3, 1.1798e-3, 1.0801e-3, 9.8896e-4, 9.0560e-4, 8.2937e-4, ...
            7.5966e-4, 6.9592e-4, 6.3765e-4, 5.8439e-4, 5.3571e-4, 4.9122e-4, ...
            4.5057e-4, 4.1342e-4, 3.7948e-4, 3.4846e-4, ...   %86-99.5 km, 0.5 km steps
            3.2011e-4, 2.7192e-4, 2.3144e-4, 1.9742e-4, 1.6882e-4, 1.4477e-4, ...
            1.2454e-4, 1.0751e-4, 9.3188e-5, 8.1142e-5, 7.1042e-5, 6.2614e-5, ...
            5.5547e-5, 4.9570e-5, 4.4473e-5, 4.0096e-5, 3.6312e-5, 3.3022e-5, ...
            3.0144e-5, 2.7615e-5, 2.5382e-5, 2.3401e-5, 2.1635e-5, 2.0055e-5, ...
            1.8635e-5, 1.7354e-5, 1.6194e-5, 1.5141e-5, 1.4183e-5, 1.3307e-5, ...
            1.2505e-5, 1.1769e-5, 1.1092e-5, 1.0468e-5, 9.8907e-6, 9.3568e-6, ...
            8.8617e-6, 8.4018e-6, 7.9739e-6, 7.5751e-6, 7.2028e-6, 6.8550e-6, ...
            6.5294e-6, 6.2243e-6, 5.9380e-6, 5.6691e-6, 5.4162e-6, 5.1781e-6, ...
            4.9538e-6, 4.7421e-6];                                  %100-149 km, 1 km steps
end
