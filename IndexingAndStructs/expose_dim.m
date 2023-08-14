function [Matrix,DimSize,DimOrder] = expose_dim(Array,Dim)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Take an n-dimensional Matlab array and reshape such that the array
%becomes 2D with the chosen dimension as the first and all other 
%dimensions merged in the second
%
%To restore the original dimensions assuming no changes to size, use this syntax:
%    Restored = permute(reshape(Matrix,DimSize),DimOrder);
%
%If changes have been made to the exposed dimension, e.g. interpolation or selection, use this syntax:
%    Restored = permute(reshape(Matrix,[size(Matrix,1),DimSize(2:end)]),DimOrder);
%
%inputs:
%  Array: n-dimensional input array
%  Dim:   dimension to bring to the front
%
%outputs:
%  Matrix:   2D output array, as described above
%  DimSize:  size of dimensions needed to put the array back using reshape()
%  DimOrder: order of dimensions needed to put the array back using permute()
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/08/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the reshaping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get size of original array
sz = size(Array);

%permute desired dimension to front
DimOrder = unique([Dim,1:1:numel(sz)],'stable');

%reshape to make all other dimensions lines
Matrix = reshape(permute(Array,DimOrder),[sz(Dim),prod(sz(DimOrder(2:end)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%work out the information we need to put it back
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DimSize is used with reshape() to put all the dimensions back to the right size
DimSize  = [size(Matrix,1),sz(DimOrder(2:end))];

%NewOrder is used with permute() to put all the dimensions back in the right order
DimOrder = 1:1:numel(sz);
DimOrder = [DimOrder(DimOrder < Dim)+1,1,DimOrder(DimOrder > Dim)];

