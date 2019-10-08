function [Error,Data,FilePath] = get_airs_granule(DataDir,DateNum,GranuleId)


  [y,~,~] =datevec(DateNum);
  dn = date2doy(DateNum);

  FilePath = [DataDir,'/',sprintf('%04d',y),'/',sprintf('%03d',dn),'/', ...
              'airs_',sprintf('%04d',y),'_',sprintf('%03d',dn),'_',sprintf('%03d',GranuleId),'.nc'];
  
  if ~exist(FilePath,'file'); Error = 1; Data = struct(); return;
  else
    Data = cjw_readnetCDF(FilePath);
    Error = 0;
  end
    
              

end

