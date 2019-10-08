function AirDensity = cjw_airdensity_mars(Height)


AirDensity =  1.89e-2 .* exp(-Height./11.1);

return

