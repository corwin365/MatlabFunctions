function [fList, pList] = resolve_all(Name,OutDir,NoClobber)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %find all functions needed by a routine
  %feed routine name to function as a string, returns a cell array of required files
  %
  %just remaps an internal Matlab capability to an easier-to-remember function name
  %Corwin Wright, c.wright@bath.ac.uk, 15/MAR/2016
  %
  %modified 2020/04/20: with a directory as a second argument, copies all
  %required files not present in that directory
  %
  %will overwrite existing files unless noclobber is set to 1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %get the list and return it (just info, no action taken)
  [fList, pList] = matlab.codetools.requiredFilesAndProducts(Name);
  if ~exist('OutDir'); return; end
  
  
  %if an outdir is specified, copy the files needed to outdir
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ~exist('NoClobber'); NoClobber = 0; end %assume clobber
    
  if NoClobber == 1;
    %get info on what is already in the output directory
    Info = dir(OutDir);
    Already.FileNames = {};
    Already.Dates     = {};
    for iFile=1:1:numel(Info)
      Already.FileNames{iFile} = Info(iFile).name;
      Already.Dates{    iFile} = Info(iFile).datenum;
    end; clear iFile Info
  end
  
  
  %loop over files. Could probably do this in an array-based way, but the
  %number of iterations is very unlikely to become large and this is more
  %readable
  for iFile = 1:1:numel(fList)
    
    %get the filename
    Slashes = [strfind(fList{iFile},'/'),strfind(fList{iFile},'\')]; %covers both unix-based and windows
    LastSlash = max(Slashes); %last slash is before start of filename
    FileName = fList{iFile}; FileName = FileName(LastSlash+1: end);
    
    %work out paths for copy
    InFile  = fList{iFile};
    OutFile = [OutDir,'/',FileName];
    
    %copy file, unless noclobber is set
    if NoClobber == 0;
      %copy it
      copyfile(InFile,OutFile)
      disp(['Copied ',FileName])
    else
      %noclobber mode
      %if the file doesn't exist, copy
      if ~any(strcmp(FileName,Already.FileNames))
        copyfile(InFile,OutFile)
        disp(['Copied ',FileName,' (no file present)'])
      else
        disp(['Skipped ',FileName,' (noclobber set)'])
      end
     
      %done! tidy up and loop
      
      clear InFile OutFile LastSlash Slashes
    end

end

