function IAGOS = prep_iagos(FilePath,ConstantSampling)


%load an IAGOS netCDF file and add:
%1. a matlab time index
%2. fill any bad data with NaNs
%3. interpolate to constant time sampling (units are in matlabtime, i.e. days)
%4. work out spatial distance between each pair of points


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IAGOS = rCDF(FilePath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute and add timestamps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get the metadata variable that gives the time basis for this file
Variables = {IAGOS.MetaData.Variables(:).Name};
VarIdx = find(strcmp(Variables,'UTC_time') ~= 0);

%go inside this variable and find the units attribute
Variables = {IAGOS.MetaData.Variables(VarIdx).Attributes(:).Name};
VarIdx2 = find(strcmp(Variables,'units') ~= 0);
Units = IAGOS.MetaData.Variables(VarIdx).Attributes(VarIdx2).Value;
clear Variables VarIdx VarIdx2

%hence, convert the dates
IAGOS.Time = datenum(Units(15:end),'yyyy-mm-dd HH:MM:SS') + IAGOS.UTC_time./60./60./24;
clear Units

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NaNify the bad data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iVar=1:1:numel(IAGOS.MetaData.Variables)
  
  %find the fill value
  Attrib= {IAGOS.MetaData.Variables(iVar).Attributes(:).Name};
  idx = find(strcmp(Attrib,'missing_value') ~= 0);
  if numel(idx) == 0; continue; end
  FillVal = IAGOS.MetaData.Variables(iVar).Attributes(idx).Value;
  
  %convert fill value to NaN
  Var = IAGOS.(IAGOS.MetaData.Variables(iVar).Name);
  Var(Var == FillVal) = NaN;
  IAGOS.(IAGOS.MetaData.Variables(iVar).Name) = Var;
  
end
clear Attrib idx FillVal iVar Var

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate to constant time sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('ConstantSampling') 
  if numel(ConstantSampling) ~= 0
    
    OldTime = IAGOS.Time;
    NewTime = min(IAGOS.Time): ConstantSampling : max(IAGOS.Time);
   
    Fields = fieldnames(IAGOS);
    for iField=1:1:numel(Fields)
      if strcmp(Fields{iField},'MetaData'); continue; end
      
      IAGOS.(Fields{iField}) = interp1(OldTime,IAGOS.(Fields{iField}),NewTime);
      
    end
    
    IAGOS.OriginalTime = OldTime;
    
  end
end

clear OldTime NewTime iField Fields


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% travel distance (in km)
%% using lat/lon, but cross-checked with speed/time for a few hundred flights and look pretty similar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = [IAGOS.lat;IAGOS.lon];
b = circshift(a,1,2);

IAGOS.dx = nph_haversine(a',b');
IAGOS.dx([1,end])  = 0;
clear a b



end

