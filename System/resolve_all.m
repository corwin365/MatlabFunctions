function FunctionList = resolve_all_v2(MainFunction,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find all functions needed by a routine and either:
% 1. just list them in a cell array
% 2. output them all to a specified directory
% 3. merge them into a single file
%
%This is a non-backwards-compatible rewrite of an older function I wrote in 2016.
%
%Corwin Wright, c.wright@bath.ac.uk, 21/DEC/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create input parser and testing functions
p = inputParser;

%inputs - required
addRequired(p,'MainFunction', @ischar                                             );  %just a type check

%inputs - optional
addParameter(p,'Action',                         1, @(x) validateattributes(x,{'numeric'},{'integer','<=',3}));
addParameter(p,'OutPath',            './RA_outdir', @ischar);
addParameter(p,'OutFile', './RA_merged_function.m', @ischar);
addParameter(p,'Clobber',                        1, @islogical);

%parse inputs and tidy up
parse(p,MainFunction,varargin{:});
Settings = p.Results;
clear p varargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get list of all child functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the list
FunctionList = matlab.codetools.requiredFilesAndProducts(MainFunction);

%if that's all we wanted, then return it
if Settings.Action == 1; return; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% put files in directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.Action == 2;

  %first, check if the directory exists, and create it if not
  if ~exist(Settings.OutPath,'dir'); mkdir(Settings.OutPath); end

  %if Clobber is set to zero, get a list of files already in that directory
  if Settings.Clobber == 0;
    Info = dir(Settings.OutPath);
    Already.FileNames = {};
    Already.Dates     = {};
    for iFile=1:1:numel(Info)
      Already.FileNames{iFile} = Info(iFile).name;
      Already.Dates{    iFile} = Info(iFile).datenum;
    end; clear iFile Info
  end

  %now copy over the files
  for iFile = 1:1:numel(FunctionList)

    %get the filename
    Slashes = [strfind(FunctionList{iFile},'/'),strfind(FunctionList{iFile},'\')]; %covers both unix-based and windows
    LastSlash = max(Slashes); %last slash is before start of filename
    FileName = FunctionList{iFile}; FileName = FileName(LastSlash+1: end);

    %work out paths for copy
    InFile  = FunctionList{iFile};
    OutFile = [Settings.OutPath,'/',FileName];

    %copy file, unless noclobber is set
    if Settings.Clobber == 1;
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




  return;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.Action == 3;

  %check we aren't clobbering
  if Settings.Clobber ~= 1 & exist(Settings.OutFile,'file');
    disp('Output file exists and clobber is not set, skipping')
  else

    %first, pull out the parent function from the list, then remove it from the list
    Parent = FunctionList(contains(FunctionList,MainFunction));
    FunctionList(contains(FunctionList,MainFunction)) = [];

    %next, load the parent function into memory
    fid = fopen(Parent{1});
    Merged = textscan(fid,'%s','Delimiter','\n');
    Merged = Merged{1};
    fclose(fid);
    disp(['Loaded ',MainFunction])


    %now, add the child functions one-by-one
    for iFile=1:1:numel(FunctionList)

      %load up the file
      fid = fopen(FunctionList{iFile});
      Text =  textscan(fid,'%s','Delimiter','\n');
      Text = Text{1};
      fclose(fid);

      %add a few lines padding
      for iX=1:1:10; Merged{end+1} = ' '; end

      %then add the text we just loaded
      Merged = [Merged;Text];
      disp(['-->Merged ',FunctionList{iFile}])

    end

    %finally, write the new file out
    disp(['All child functions of ',MainFunction,' merged'])
    writelines(Merged,Settings.OutFile)
    warning('Some merged functions may need ''end'' adding manually. These should be highlighted by the error finder in the Matlab editor')


  end


  return;

end

