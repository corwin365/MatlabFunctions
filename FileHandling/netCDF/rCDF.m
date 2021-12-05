function FileContents = rCDF(FileName,OldFormat)

if ~exist('OldFormat'); OldFormat = 0; end

if OldFormat == 1
  
  %use my old method - doesn't do a lot of useful things, but more robust for
  %some uses (not many, it's very old...)
  
  
  %open the file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  NetCdfID = netcdf.open(FileName,'NOWRITE'); %open file, read-only access
  
  %find list of variables in the file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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
  
  FileContents.MetaData = ncinfo(FileName);
  
  Result = 1; %finished
  
  netcdf.close(NetCdfID); % done
  
  
  
else %use neil's method
  
  %uses Neil's getnet() routine as it's better, but modified to provide the same
  %outputs as my old function
  
  Data = nph_getnet(FileName);
  
  MetaFields = {'Filename','Name','Dimensions','Variables','Attributes','Groups','Format'};
  FileContents = Data.Data;
  
  for iField=1:1:numel(MetaFields)
    FileContents.MetaData.(MetaFields{iField}) = Data.(MetaFields{iField});
  end


end