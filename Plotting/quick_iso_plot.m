function quick_iso_plot(x,y,z,v,Surface,Colour,EdgeColour)

if ~exist('EdgeColour'); EdgeColour = 'none'; end


fv = isosurface(x,y,z,v,Surface);
ThePatch = patch(fv);
ThePatch.FaceColor = Colour;
ThePatch.EdgeColor = EdgeColour;
ThePatch.FaceAlpha  = 0.33;




end

