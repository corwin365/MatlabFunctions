function lon = wrapTo90(lon)

positiveInput = (lon > 0);
lon = mod(lon, 90);
lon((lon == 0) & positiveInput) = 90;
