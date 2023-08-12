function Array = index_dim(Array,Indices,Dim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select elements of an array along a selected dimension
%Corwin Wright, c.wright@bath.ac.uk, 2023/08/12 
%
%inputs:
%  Array   - the array we want to operate on
%  Dim     - the dimension we want to operate on
%  Indices - the indices we want to extract along dimension Dim
%
%outputs:
%  Array   - the selected array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input handling
if nargin      < 3; Dim = 1; end          %assume first dimension if not specified, and set varargin to blank if not set
if numel(Dim) == 0; Dim = 1; end     %if set to blank, then set first dimension

%expose the desired dimension
[y,a,b] = expose_dim(Array,Dim);

%select
y = y(Indices,:);

%put the dimensions back in order 
Array = permute(reshape(y,[size(y,1),a(2:end)]),b);

return
