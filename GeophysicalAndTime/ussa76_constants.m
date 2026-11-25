function S = ussa76_constants()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%constants of the U.S. Standard Atmosphere, 1976, shared by h2p and p2h
%
%Source: U.S. Standard Atmosphere, 1976, NOAA-S/T 76-1562 / NASA-TM-X-74335, Washington, D.C., October 1976.
%Page numbers are the document's printed pages. All values below were checked against it on 2026/11/25.
%
%outputs (struct S):
%   g0, Rs, M0, r0, T0, P0 - g0' (m2 s-2 m'-1), R* (N m kmol-1 K-1), M0 (kg kmol-1), r0 (km), T0 (K), P0 (hPa);
%                            table 2 (p. 2) and p. 9
%   Hb, L                  - layer base geopotential heights (km', 8 values) and molecular-scale temperature
%                            gradients (K/km', 7 values); table 4 (p. 3)
%   Tb, Pb                 - base temperatures (K) and pressures (hPa) of the 7 layers, from eqs (23), (33a), (33b)
%   k                      - g0' M0 / R*, in K per km'
%   MinAlt                 - lowest valid geometric altitude, km (start of Table I, p. 51)
%   TableZ, TableP         - Table I (geometric altitude, p. 68): Z (km) and P (hPa), 86-149 km. Each value was
%                            checked against the table's independent P/P0 column.
%
%Corwin Wright, c.wright@bath.ac.uk, and Claude, 2026/11/25
%
%change log:
%  2026/11/25 [Claude] created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.g0 = 9.80665;      %g0', m2 s-2 m'-1
S.Rs = 8.31432e3;    %R*, N m kmol-1 K-1
S.M0 = 28.9644;      %M0, kg kmol-1
S.r0 = 6356.766;     %r0, km
S.T0 = 288.15;       %K
S.P0 = 1013.25;      %hPa

S.Hb = [0, 11, 20, 32, 47, 51, 71, 84.8520];      %km'
S.L  = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0];    %K/km'
S.k  = S.g0.*S.M0./S.Rs.*1e3;                     %K/km'

%base temperatures and pressures of each layer
S.Tb = S.T0; S.Pb = S.P0;
for b=1:1:numel(S.L)-1
  dH = S.Hb(b+1) - S.Hb(b);
  S.Tb(b+1) = S.Tb(b) + S.L(b).*dH;
  S.Pb(b+1) = ussa76_layer_pressure(S.Pb(b),S.Tb(b),S.L(b),dH,S.k);
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

function P = ussa76_layer_pressure(Pb,Tb,L,dH,k)
%pressure dH (km') above a layer base: eq (33a), or eq (33b) if the gradient L is zero
if L == 0; P = Pb.*exp(-k.*dH./Tb);
else;      P = Pb.*(Tb./(Tb + L.*dH)).^(k./L);
end
end
