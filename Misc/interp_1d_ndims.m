function yi= interp_1d_ndims(x,y,xi,Dim,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do 1d interpolation along a chosen dimension of an ND array
%Corwin Wright, c.wright@bath.ac.uk, 2022/08/06
%
%inputs: as interp1, except with the dimension specified as the last non-varargin input (set to [] to skip)
%outputs: array interpolated along the chosen dimension
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input handling
if nargin < 4; Dim = 1; end          %assume first dimension if not specified, and set varargin to blank if not set
if numel(Dim) == 0; Dim = 1; end     %if set to blank, then set first dimension
if nargin < 5; varargin = {}; end    %make sure we have a varargin for the call to interp1 below

%get size of array
sz = size(y);

%permute desired dimension to front
DimOrder = unique([Dim,1:1:numel(sz)],'stable');

%reshape to make all other dimensions lines
y = reshape(permute(y,DimOrder),[sz(Dim),prod(sz(DimOrder(2:end)))]);

%interpolate
yi = interp1(x,y,xi,varargin{:});

%reshape back
yi = reshape(yi,[numel(xi),sz(DimOrder(2:end))]);

%and permute back
NewOrder = 1:1:numel(sz);
NewOrder = [NewOrder(NewOrder < Dim)+1,1,NewOrder(NewOrder > Dim)];
yi = permute(yi,NewOrder);


return
