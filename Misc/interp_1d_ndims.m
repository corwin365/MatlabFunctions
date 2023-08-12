function yi= interp_1d_ndims(x,y,xi,Dim,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do 1d interpolation along a chosen dimension of an ND array
%Corwin Wright, c.wright@bath.ac.uk, 2022/08/06
%
%updated 2023/08/12 to move reshaping and reordering to a separate function expose_dim() for more general use
%
%inputs: as interp1, except with the dimension specified as the last non-varargin input (set to [] to skip)
%outputs: array interpolated along the chosen dimension
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input handling
if nargin < 4; Dim = 1; end          %assume first dimension if not specified, and set varargin to blank if not set
if numel(Dim) == 0; Dim = 1; end     %if set to blank, then set first dimension
if nargin < 5; varargin = {}; end    %make sure we have a varargin for the call to interp1 below

%expose the desired dimension
[y,a,b] = expose_dim(y,Dim);

%interpolate
yi = interp1(x,y,xi,varargin{:});

%put the dimensions back in order 
yi = permute(reshape(yi,[size(yi,1),a(2:end)]),b);

return
