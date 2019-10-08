function [Results] = airs3d_t_v2_mat(Airs,Swath,Altitude,Detrend,ComputeAtten)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute appropriately-weighted 3D AIRS brightness temperatures
%uses daily mat files rather than raw HDF files
%Corwin Wright, corwin.wright@trinity.oxon.org
%12/NOV/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select the swath we actually want
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AirsSwath.radiances    = squeeze(Airs.radiances(   Swath,:,:,:));
AirsSwath.Latitude     = squeeze(Airs.Latitude(    Swath,:,:));
AirsSwath.Longitude    = squeeze(Airs.Longitude(   Swath,:,:));
AirsSwath.Time         = squeeze(Airs.Time(        Swath,:,:));
AirsSwath.CalFlag      = squeeze(Airs.CalFlag(     Swath,:,:));
AirsSwath.state        = squeeze(Airs.state(       Swath,:,:));
AirsSwath.nominal_freq = squeeze(Airs.nominal_freq(Swath,:));
AirsSwath.Channels     = Airs.Channels;
clear AirsMat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specify Gong et al channels, which we use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gong.Pressures = [2,2.5,3,4,7,10,20,30,40,60,80,100];
Gong.Levels    = zeros(numel(Gong.Pressures),14);
Gong.Levels(1,1)     = 74; %2.0 hPa
Gong.Levels(2,1)     = 75; %2.5 hPa
Gong.Levels(3,1)     = 76; %3.0 hPa
Gong.Levels(4,1)     = 77; %4.0 hPa
Gong.Levels(5,1)     = 78; %7.0 hPa
Gong.Levels(6,1)     = 79; %10.0 hPa
Gong.Levels(7,1:2)   = [81,82]; %20hPa
Gong.Levels(8,1:6)   = [102, 108, 114, 120, 125, 126 ] ;%30 hPa
Gong.Levels(9,1:7)   = [64, 88, 90, 94, 100, 106, 118]; %40 hPa
Gong.Levels(10,1:9)  = [66, 68, 70, 86, 87, 91, 93, 97, 130 ]; %60hPa
Gong.Levels(11,1:14) = [92, 98, 104, 105, 110, 111, 116, 117, 122, 123, 128, 129, 134, 140]; %80hPa
Gong.Levels(12,1:6)  = [132, 133, 138, 139, 149, 152 ]; %100 hPa

Gong.MinDetectable = sqrt(1e-3.*[3.78,3.72,3.63,3.66,3.88,4.62,2.14,0.98,0.83,0.66,0.50,0.67]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check quality flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
  if numel(BadData) ~= 0; %takes up significant time and often does nothing!

    for iX=1:1:numel(AirsSwath.Latitude(:,1));
      
      %final check on data presence in file
      try
        Temporary = squeeze(AirsSwath.radiances(:,iX,:));
      catch
%         ReturnBT   = NaN;    ReturnLat  = NaN;    ReturnLon  = NaN;    ReturnTime = NaN;    
        return;
      end
      
      Temporary(BadData) = NaN;
      AirsSwath.radiances(:,iX,:) = Temporary;
    end; clear Temporary BadData iX;
  end

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

if numel(BadData) ~= 0;
  
  for iChannel=1:1:NChannels; %NeN is just an arbitrary property with NChannels elements
    Temporary = squeeze(AirsSwath.radiances(iChannel,:,:));
    Temporary(BadData) = NaN;
    AirsSwath.radiances(iChannel,:,:) = Temporary;
  end; 
end
  
clear Temporary BadData iChannel NChannels;

% Individual channel readings ("radiances") must be checked for the flag
% bad value  of  –9999.0. A channel  reading  is  set  to  this  value  only
% when  no  radiance  can  be  calculated;  questionable  or suspect values
% are indicated only by QA fields

AirsSwath.radiances(AirsSwath.radiances == -9999.0) = NaN;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop over height levels, and compute the brightness temperature for each one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iLevel=1:1:numel(Altitude);
  
  %which channels do we want to process?
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %use nearest level from Gong et al (2012)
  %check we're in the specified height range (slight margins either side)
  PressureLevel = alt2pres_complex(Altitude(iLevel)./1000);
  
  if     PressureLevel > 110;   continue; %too far away
  elseif PressureLevel <   1.5; continue; %too far away
  else
    %work out the channels we need
    [~,idx2] = min(abs(Gong.Pressures - PressureLevel));
    Channels = Gong.Levels(idx2,:);
    
    MinDetectable = Gong.MinDetectable(idx2);
    
  end
  
  Channels = Channels(Channels > 0);
  if numel(Channels) == 0; continue; end
  
  %map the channels from the original format to the storage format
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Channels2 = Channels.*NaN;
  for iCh = 1:1:numel(Channels);
    Channels2(iCh) = find(Channels(iCh) == AirsSwath.Channels);
  end; clear iCh
  Channels = Channels2; clear Channels2

  %extract radiances and frequencies for channels
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Radiances      = AirsSwath.radiances(Channels,:,:);
  ChannelFreqs   = AirsSwath.nominal_freq(Channels);
  ChannelWeights = ones(numel(Channels),1)./numel(Channels); %i.e. all equal


  %compute brightness temperature for each location and channel
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  BT = NaN(size(Radiances));
  
  for iChannel = 1:1:size(Radiances,1);

    %work out brightness temperature for each point
    R = squeeze(Radiances(iChannel,:,:));
    W = squeeze(ChannelFreqs(iChannel));
    BT(iChannel,:,:) = cjw_brightnesstemperature(W,R);
    
    %weight
    BT(iChannel,:,:) = BT(iChannel,:,:) .* ChannelWeights(iChannel);
    
  end; clear iChannel

  
  %sum weighted channels over channel dimension, to get final results
  BT = squeeze(nansum(BT,1)); %sum is over first non-singleton dimension, i.e. channels, which are already weighted
  
  %store
  Results.Swath(iLevel,:,:) = BT;
  
  %tidy up
  clear BT Radiances ChannelFreqs ChannelWeights
  
end


if nansum(Results.Swath) == 0;
  Results = 'no data found';
  return
else
  Results.Lat  = AirsSwath.Latitude;
  Results.Lon  = AirsSwath.Longitude;
  Results.Time = AirsSwath.Time;
end



Results.OriginalT = squeeze(Results.Swath); %T as opposed to T'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%work out attenuation responses for the channels used
%only do if necessary - big bottleneck
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('ComputeAtten') == 0; ComputeAtten = 1; end %compute by default
if ComputeAtten ~= 0 ;

  for iChannel = 1:1:numel(Channels);
    %get attenuation response for the channel(s) used
    [Attenuation] = compute_airs_attenuation(Channels(iChannel));
    
    if iChannel == 1;
      AllAtt.Attenuation  = Attenuation.Attenuation;
      AllAtt.LambdaZScale = Attenuation.LambdaScale;
    else
      AllAtt.Attenuation( iChannel,:) = Attenuation.Attenuation;
      AllAtt.LambdaZScale(iChannel,:) = Attenuation.LambdaScale;
    end
    
    AllAtt.FWHM = nanmean(Attenuation.FWHM(:));
  end
  
  Results.Attenuation = AllAtt; clear AllAtt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%detrend?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Detrend ~=0;
  %detrend according to Alexander and Grimsdell (2013)
  
  for iLevel=1:1:numel(Altitude);
    AcrossTrackElements=1:1:size(Results.Swath,2); %index array, purely to clarify below loop
  
    for iStripe=1:1:size(Results.Swath,3);

      %select cross-track row
      TheStripe              = squeeze(Results.Swath(iLevel,:,iStripe));
      %compute and remove the polynomial, leaving BT perturbations
      ThePolynomial  = polyval(polyfit(AcrossTrackElements,TheStripe,4),AcrossTrackElements);
      %remove polynomial from data
      Results.Swath(iLevel,:,iStripe) = TheStripe-ThePolynomial;

    end; clear iStripe TheStripe ThePolynomial;
  end
  
end

Results.MinDetectable = MinDetectable;


return
end
