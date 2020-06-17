function Y = smoothn2(X,sz,filt,std)

if nargin == 2;
  filt = 'b';
elseif nargin == 3;
  std = 0.65;
end

Bad = find(isnan(X));
X = inpaint_nans(X);

if nargin ==2;X = smoothn(X,sz,filt); 
else          X = smoothn(X,sz,filt,std); 
end

X(Bad) = NaN;
Y = X;