function Arr = make_odd(Arr,Down)

  %which way?
  if exist('Down','var'); Shift = -1; else Shift = 1; end

  %integer?
  for iEl=1:1:numel(Arr);  if ~isinteger(Arr(iEl)); Arr(iEl) = round(Arr(iEl)); end 

  %odd number?
  NotOdd = find(mod(Arr,2) == 0);

  Arr(NotOdd) = Arr(NotOdd) + Shift;

end

