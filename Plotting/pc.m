function pc(a,b,c)

%shorthand version of pcolor, as i type this exact sequence of commands far too often

if nargin == 1;
  pcolor(squeeze(a)); shading flat; colorbar
else
  pcolor(squeeze(a),squeeze(b),squeeze(c)); shading flat; colorbar
end


end

