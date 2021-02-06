function Arr = make_odd(Arr,Down)

  if exist('Down','var'); Shift = -1; else Shift = 1; end

  NotOdd = find(mod(Arr,2) == 0);

  Arr(NotOdd) = Arr(NotOdd) + Shift;

end

