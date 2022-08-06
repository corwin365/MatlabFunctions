function yi= interp_1d_ndims(x,y,xi,Dim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do 1d interpolation along a chosen dimension of an ND array
%Corwin Wright, c.wright@bath.ac.uk, 2022/08/06
%
%inputs: as interp1, except with the dimension specified as the last input
%outputs: array in original form except for the specified dimension
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%get size of array
sz = size(y);

%permute desired dimension to front
DimOrder = unique([Dim,1:1:numel(sz)],'stable');

%reshape to make all other dimensions lines
y = reshape(y,[sz(Dim),prod(sz(DimOrder(2:end)))]);

%interpolate
yi = interp1(x,y,xi);

%reshape back
yi = reshape(yi,[numel(xi),sz(DimOrder(2:end))]);

%and permute back
NewOrder = 1:1:numel(sz);
NewOrder = [NewOrder(NewOrder < Dim)+1,1,NewOrder(NewOrder > Dim)];
yi = permute(yi,NewOrder);


return yi

end
