function y = reverse(y,Dim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reverse values in Matlab array
%
%
%inputs:
%  y                           - array to be reversed
%  Dim (optional, default = 1) - dimension to revers ealong
%
%Set Dim to zero to reverse all values regardless of dimension.
%
%outputs: 
%  y - reversed array
%
%Corwin Wright, c.wright@bath.ac.uk
%original age unknown as  it was uncommented - bad Corwin.
%Updated 2023/08/14 to allow operation along an arbitrary dimension and commented
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('Dim','var'); Dim = 1; end

if Dim == 0; 
  %reverse all values regardless of dimension
  %make sure you know what you want going in...
  sz = size(y);
  y = y(:);
  y = y(end:-1:1);
  y = reshape(y,sz);
else
  %expose the dimension we want as dimension 1, all others merged
  [y,a,b] = expose_dim(y,Dim);
  %reverse dimension 1
  y = y(end:-1:1,:);
  %put back into original order and shape
  y = permute(reshape(y,[size(y,1),a(2:end)]),b);
end
