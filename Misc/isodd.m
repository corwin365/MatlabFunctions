function tf = isodd(x)
% isodd - Returns true if the input is an odd integer
% Syntax: tf = isodd(x)
% x - numeric (real) input of any number of dimensions and sizes
% tf - true for each element that is an odd integer
%   Example
%  tf = isodd([1,3,4]);
%
%   See also: iseven

% AUTHOR    : Dan Kominsky
% Copyright 2012  Prime Photonics, LC.
%%
  if ~isreal(x)
    error('iseven:badinput','iseven requires real inputs');
  else
    tf = mod(x,2)==1;
  end
end