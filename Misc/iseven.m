function tf = iseven(x)
  % iseven - Returns true if the input is an even whole number
  % Syntax: tf = iseven(x)
  % x - numeric (real) input of any number of dimensions and sizes
  % tf - true for each element that is an even whole number
  %   Example
  %  tf = iseven([1,3,4]);
  %
  %   See also: isodd
  
  % AUTHOR    : Dan Kominsky
  % Copyright 2012  Prime Photonics, LC.
  if ~isreal(x)
    error('iseven:badinput','iseven requires real inputs');
  else
    tf = mod(x,2)==0;
  end
end