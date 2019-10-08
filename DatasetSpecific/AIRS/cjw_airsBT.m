function [ReturnBT,ReturnLat,ReturnLon,ReturnTime] = cjw_airsBT(AirsFile,NetCDF,Channels,Detrend)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute detrended brightness temperature for an AIRS granule
%reads in AIRIBRAD level 1 granules
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%05/JAN/2014
%
%inputs
%---------
%
%AirsFile - (cell   1 n) filename (including path) of granule data
%NetCDF   = data is stored in a NetCDf file rather than a HDF file (1 yes 0 no, optional; default 0)
%Channels - (double 1 n) channels to average  (optional: default from Hoffman et al 2010)
%Detrend  - (double 1  ) detrend x-track with 4th-order polynomial (optional: 1 yes, 0 no, default 1)
%
%outputs
%---------
%
%ReturnBT   - (double x y) Brightness temperature for (x,y) location in granule
%ReturnLat  - (double x y) Latitude               for (x,y) location in granule
%ReturnLon  - (double x y) Longitude              for (x,y) location in granule
%ReturnTime - (double x y) Timestamp              for (x,y) location in granule
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%input handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4;Detrend = 1; end; %detrend by default
if nargin<3; %use default channels (Hoffman et al 2010)
  Channels = 2040;
  for iChannel=2041:1:2065; Channels(end+1) = iChannel;end;clear iChannel;
  for iChannel=2072:1:2087; Channels(end+1) = iChannel;end;clear iChannel;
end;
if nargin<2; %assume data in HDF
    NetCDF = 0;
end

%read in data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if   NetCDF == 0 ;
  %HDF file
  try;
    AirsSwath = read_airs_swath(AirsFile,3);
  catch;
    ReturnBT   = NaN;    ReturnLat  = NaN;    ReturnLon  = NaN;    ReturnTime = NaN;    return;
  end;
else
  %netCDF file
  AirsSwath = cjw_readnetCDF(AirsFile);
  %correct for small difference in file generation script to HDF file
  AirsSwath.radiances = AirsSwath.radiance;
  AirsSwath = rmfield(AirsSwath,'radiance');
end

clear AirsFile;

%check quality flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%from AIRIBRAD readme file:

% The AIRS  L1B product  contains a per-scan  field named “CalFlag”.  Users
% should  avoid using any channel for any scan in which the "offset problem"
% or "gain problem", or "pop detected" bits are set (bits 6, 5, and 4 
% respectively where bit 0 is LSB). Bit 0, “cold scene noise”, and bit 1, 
% “telemetry out of limit condition”, indicates conditions that can
%potentially impact data quality. Users who require pristine data should
%discard any data in which either of these bits is set.

for iFlag=4:1:6;
  BadData = find(AirsSwath.CalFlag == iFlag);

  for iX=1:1:90;
    
    %final check on data presence in file
    try
      Temporary = squeeze(AirsSwath.radiances(:,iX,:));
    catch
      ReturnBT   = NaN;    ReturnLat  = NaN;    ReturnLon  = NaN;    ReturnTime = NaN;    return;
    end
    
    Temporary(BadData) = NaN;
    AirsSwath.radiances(:,iX,:) = Temporary;
  end; clear Temporary BadData iX;

end; clear iFlag;

% Before using  any AIRS  L1B radiance, check the value of the
% corresponding  “state” to ensure that it is equal to zero. There
% is one “state” value per  field-of-view (FOV), and it is valid for all 
% all channels in that FOV. The “state” valids and their meaning are: 
%
% Process    0  normal data
% Special    1  instrument in special calibration mode when these data were
%               taken (e.g., staring at nadir)
% Erroneous  2  data known bad (e.g., instrument in safe mode)
% Missing    3  data are missing

NChannels = numel(AirsSwath.CalFlag(:,1));

BadData = find(AirsSwath.state ~= 0); %remove anything outside 'process' state
for iChannel=1:1:NChannels; %NeN is just an arbitrary property with NChannels elements
  Temporary = squeeze(AirsSwath.radiances(iChannel,:,:));
  Temporary(BadData) = NaN;
  AirsSwath.radiances(iChannel,:,:) = Temporary;
end; clear Temporary BadData iChannel;

% Individual channel readings ("radiances") must be checked for the flag
% bad value  of  –9999.0. A channel  reading  is  set  to  this  value  only
% when  no  radiance  can  be  calculated;  questionable  or suspect values
% are indicated only by QA fields

BadData = find(AirsSwath.radiances == -9999.0); 
AirsSwath.radiances(BadData) = NaN;
clear BadData;
  
  
%extract data we want
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if the file contains a variable called 'Channels', take account of that
%and select the channels we want from that
%otherwise, assume channels are numbered sequentially from 1
if isfield(AirsSwath,'Channels');
  NChannelsWanted = numel(Channels);
  ChannelsWanted  = zeros(NChannelsWanted,1);
  for i=1:1:NChannelsWanted;
      MatchingChannels = find(intersect(AirsSwath.Channels,Channels) == Channels(i));
      NFound = numel(MatchingChannels);
      if NFound == 1; ChannelsWanted(i) = MatchingChannels;
      else
          disp(char(strcat('Error - channel ',num2str(Channels(i)),' not found in dataset. Continuing without...')));
          ChannelsWanted(i) = NaN;
      end
  end
  Channels = ChannelsWanted(isfinite(ChannelsWanted));
end;


%reduce down to just the channels we want
AirsData.Latitude    = AirsSwath.Latitude;
AirsData.Longitude   = AirsSwath.Longitude;
AirsData.Radiances   = AirsSwath.radiances(Channels,:,:);
AirsData.ChannelFreq = AirsSwath.nominal_freq(Channels);
AirsData.Time        = AirsSwath.Time;
clear AirsSwath Channels;

%compute brightness temperature for each location and channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AirsData.BT = NaN(size(AirsData.Radiances));

for iChannel = 1:1:size(AirsData.Radiances,1);
  for iCrossTrack=1:1:size(AirsData.Radiances,2);
    for iAlongTrack=1:1:size(AirsData.Radiances,3);
      Radiance   = AirsData.Radiances(iChannel,iCrossTrack,iAlongTrack);
      Wavenumber = AirsData.ChannelFreq(iChannel);
      AirsData.BT(iChannel,iCrossTrack,iAlongTrack) = cjw_brightnesstemperature(Wavenumber,Radiance);
    end; clear iAlongTrack Radiance Wavenumber;
  end; clear iCrossTrack;
end; clear iChannel

%average channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AirsData.BT = squeeze(nanmean(AirsData.BT)); %mean is over first non-singleton dimension, i.e. channels


%detrend with fourth-order polynomial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Detrend == 1;
  
  AcrossTrackElements=1:1:size(AirsData.Latitude,1); %index array, purely to clarify below loop
  
  for iStripe=1:1:size(AirsData.Latitude,2);

    %select cross-track row
    TheStripe              = squeeze(AirsData.BT(:,iStripe))';
    %compute and remove the polynomial, leaving BT perturbations
    ThePolynomial  = polyval(polyfit(AcrossTrackElements,TheStripe,4),AcrossTrackElements);
    %remove polynomial from data
    AirsData.BT(:,iStripe) = TheStripe-ThePolynomial;

   end; clear iStripe TheStripe ThePolynomial;
end; clear Detrend AcrossTrackElements; %if Detrend == 1

%return variables and clear memory
ReturnBT   = AirsData.BT;
ReturnLat  = AirsData.Latitude;
ReturnLon  = AirsData.Longitude;
ReturnTime = AirsData.Time;
clear AirsData

%done!
