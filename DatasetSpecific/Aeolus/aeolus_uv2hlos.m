function HLOS = aeolus_uv2hlos(u,v,theta)

  HLOS = -u.*sind(theta) - v.*cosd(theta);


end
