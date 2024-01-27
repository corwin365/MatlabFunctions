function FileContents = rCDF(FilePath,Method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read netCDF file into Matlab
%Corwin wright, pre-2015
%COMPATABILITY-BREAKING update 2024/01/27 - see below
%
%inputs:
%  FilePath: path to file
%  Method:   0 - [default] call Neil's nph_getnet() but switch order of MetaData
%                and real data, and reverse variable dimension order to match 
%                internal netCDF dimension ordering
%            1 - simple loop over list of fields (simple, but robust)
%            2 - as 0, but with the nph_getnet() variable order 
%
%IMPORTANT WARNING RE: BACKWARDS COMPATIBILITY:
%  Method [1] was the default before ~2017
%  Method [2] was the default before 2024/01/27
%  [0] is NOT BACKWARDS-COMPATIBLE for any variable >1D
%  old code must be updated to add a method flag of '1' or '2' before use.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set default method
if ~exist('Method'); Method = 0; end



if Method == 0 | Method == 2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %based around Neil Hindley's nph_getnet()
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %run getnet
  Data = nph_getnet(FilePath);

  %move data up to top level
  FileContents = Data.Data;

  %move metadata to subsidiary level
  MetaFields = {'Filename','Name','Dimensions','Variables','Attributes','Groups','Format'};
  for iField=1:1:numel(MetaFields); FileContents.MetaData.(MetaFields{iField}) = Data.(MetaFields{iField}); end


  if Method == 0;
    %reorder dimensions so that they're in netCDF order (i.e. reverse them)
    Fields = fieldnames(FileContents);
    for iField=1:1:numel(Fields)
      if strcmp(Fields{iField},'MetaData'); continue;
      else;     
        sz = size(FileContents.(Fields{iField}));
        if numel(sz) > 2 | (numel(sz) == 2 & (sz) == 1)
          FileContents.(Fields{iField}) = permute(FileContents.(Fields{iField}),numel(sz):-1:1);
        end
      end
    end; clear iField
  end



elseif Method == 1
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %my old method - doesn't do a lot of useful things, but a little more robust 
  %in some special cases (not many, it's very old...)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %open the file  
  NetCdfID = netcdf.open(FilePath,'NOWRITE'); %open file, read-only access
  
  %find list of variables in the file 
  [~,NumVars,~,~] = netcdf.inq(NetCdfID);
  
  for iVar=1:1:NumVars;
    [VarName,xtype,dimids,natts] = netcdf.inqVar(NetCdfID,iVar-1);
    
    
    %remove characters which should be in the variable names
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %check the varname for spaces.
    %replace them with _s
    spaces = isspace(VarName);
    VarName(spaces) = '_';
    
    %check the varname for brackets.
    %replace them with _s
    brackets = strfind(VarName,'(');
    VarName(brackets) = '_';
    brackets = strfind(VarName,')');
    VarName(brackets) = '_';
    
    %ok, go
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Command = strcat(['FileContents.',VarName,' = netcdf.getVar(',num2str(NetCdfID),',',num2str(iVar-1),');']);
    eval(Command);

  end
  
  FileContents.MetaData = ncinfo(FilePath);  
  netcdf.close(NetCdfID); % done
  
end