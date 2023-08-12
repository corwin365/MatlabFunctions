function MinMax = minmax(Array)

  TheMin = nanmin(Array(:));
  TheMax = nanmax(Array(:));
  
  MinMax = [TheMin,TheMax];
  return

end

