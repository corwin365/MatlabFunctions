
%wrapper for alt2pres_complex

function p = h2p(alt)


  %rearrange data into a line
  sz = size(alt);
  alt = alt(:);
  
  p = alt2pres_complex(alt);
  
  %and convert back
  p = reshape(p,sz);
    
  
end % end cycling through all of alt's element(s)...



