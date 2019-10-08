function [SaberData,OldFile] = extract_saber_data(MatlabDay,DataDir,OldFile,NoFlatten)

if ~exist('NoFlatten'); NoFlatten = 0; end
if ~exist('OldFile'); OldFile.Name = ''; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract SABER data for a given day
%only runs if calling a new file to previous iteration
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%14/FEB/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%errors:
%0. no error - successful
%1. no file found
%2. no data on this day
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the file (if it's not the one from last time) and load it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(OldFile,'Name'); OldFile.Name = ' '; end; %always load file afresh if old file not specified

%identify file for this month
[y,m,~]   = datevec(MatlabDay);
FileString = strcat(cjw_monthname(m),num2str(y),'_v2.0.mat');
if strcmp(FileString,OldFile.Name) == 0;
  %new file - load it up
  
  FileName = wildcardsearch(DataDir,FileString);
  
  if numel(FileName) == 0;  %no file
    SaberData.Error = 1;
    return;
  end;
  
  %file exists! load data
  AllSaberData = load(FileName{1}); AllSaberData = AllSaberData.SaberData;
  AllSaberData = cleandata_saber(AllSaberData);
  clear FileName;
  
  DayScale = floor(nanmin(AllSaberData.MatlabTime,[],1)); %strictly, day of profile start
  
  %convert density units
  AllSaberData.Dens = 1.225.*(AllSaberData.Dens/0.02504e21); %kg.m^-3
  
  OldFile.Name = FileString;
  
else
  %same as last call - don't reload
  AllSaberData = OldFile.Data;
  DayScale = floor(nanmin(AllSaberData.MatlabTime,[],1)); %strictly, day of profile start
end
clear m y DataDir FileString

OldFile.Data = AllSaberData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract the day we actually want
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OnThisDay = find(DayScale == MatlabDay);
if numel(OnThisDay) ==0;
  SaberData.Error = 2;
  return;
end
clear MatlabDay

%subset needed data
SaberData.Temp    = AllSaberData.Temp(   :,OnThisDay);
SaberData.Prs     = nanmean(AllSaberData.Prs,2);
SaberData.Lat     = AllSaberData.Lat(    :,OnThisDay);
SaberData.Lon     = AllSaberData.Lon(    :,OnThisDay);
SaberData.Height  = AllSaberData.Height( :,OnThisDay);
SaberData.Dens    = AllSaberData.Dens(   :,OnThisDay);
SaberData.Time    = nanmin(AllSaberData.MatlabTime,[],1);

if NoFlatten ==0;
  %produce lat and lon scale
  SaberData.Lat = nanmean(SaberData.Lat,1);
  SaberData.Lon = nanmean(SaberData.Lon,1);
end

SaberData.Error = 0; %it worked!

return
end
