function Out = zp(Values,Direction)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to convert between pressure and height based on the assumptions
%implicit in the US Standard Atmosphere. 
%
%cleaned up version of alt2pres_complex and pres2alt_complex, originally 
%written by Neil
%
%inputs: 
%  Values:     array of values to be converted
%
%  Direction: 'z' to go from pressure to height
%             'p' to go from height to pressure
%
%
%In both directions, pressure is in HECTOPASCALS and height is in KILOMETRES
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/03/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check that the direction is specifically either 'z' or 'p'
if ~exist('Direction','var')
  disp('Error in zp.m - direction not set. Returning NaNs.')
  Out = Values.*NaN;
  return
end
if ~strcmp(Direction,'z') && ~strcmp(Direction,'p');
  disp('Error in zp.m - invalid direction, should be ''z'' to produce heights or ''p'' to produce pressures. Returning NaNs.' )
  Out = Values.*NaN;
  return
end


%% set up the assumptions needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%layer boundaries in pressure (hPa)
Basis.P = [1013.25, 226.3210, 54.7489, 8.6802, 1.1091, 0.6694, 0.0396];

%layer boundaries in height (km)
Basis.Z = [0, 11, 20, 32, 47, 51, 71];

%compute scale heights for each layer
H = -(Basis.Z(2:end)-Basis.Z(1:end-1)) ./ log(Basis.P(2:end)./Basis.P(1:end-1));
H(end+1) = H(end);


%% find the values in the lookup array needed to do the calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assume any pressures >1013.25 or altitudes < 0 are at ground level. 
if     strcmp(Direction,'z');  Values(Values > max(Basis.P)) = max(Basis.P);
elseif strcmp(Direction,'p');  Values(Values < min(Basis.Z)) = min(Basis.Z);
end

%choose a direction
if     strcmp(Direction,'z');  Bs = Basis.P;
elseif strcmp(Direction,'p');  Bs = Basis.Z;
end

%make the values and basis array the same size and shape
A = repmat(Values,[ones(ndims(Values),1)',size(Bs)]);
B = repmat(Bs',[1,size(Values)]);
B = permute(B,[2:ndims(A),1]);

%hence find the closest value in the Basis array  BELOW the current level (as we're working in layers)
[~,idx] = min(abs(A-B),[],ndims(A));
delta = Values-Bs(idx);
idx(delta < 0 & idx > 1) = idx(delta < 0 & idx > 1) - 1;

%finally, do the calculation
if     strcmp(Direction,'z');  Out = (-H(idx) .* ln(Values./Basis.P(idx))) + Basis.Z(idx);
elseif strcmp(Direction,'p');  Out = Basis.P(idx) .* exp(-(Values - Basis.Z(idx))./H(idx));
end


%done!
