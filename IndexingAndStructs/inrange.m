function InRange = inrange(Array,MinMax,NoEnds)

  if nargin < 3; NoEnds = 0; end
  
  if NoEnds ~= 1;
    InRange = find(Array >= min(MinMax) & Array <= max(MinMax));
  else
    InRange = find(Array >  min(MinMax) & Array <  max(MinMax));
  end
  return

end

