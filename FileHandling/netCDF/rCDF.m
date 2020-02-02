function FileContents = rCDF(FileName,OldFormat)

%alias for cjw_readnetCDF to reduce typing - no internal functionality.

if ~exist('OldFormat'); OldFormat = 0; end

FileContents = cjw_readnetCDF(FileName,OldFormat);


end