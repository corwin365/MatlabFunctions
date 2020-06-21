
%wrapper for pres2alt_complex
function z = p2h(pres)


  %rearrange data into a line
  sz = size(pres);
  pres = pres(:);
  
  
  %process
  z = pres2alt_complex(pres);

  
  %and convert back
  z = reshape(z,sz);
  


end 





