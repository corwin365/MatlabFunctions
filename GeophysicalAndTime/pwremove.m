function [BGTemp,Perturbations,PWs] = pwremove(SatelliteData,PWLat,NPWs,LatRange,Verbose)

%this function has now been replaced by the improved and faster 'pwfilter.m'
%old code is kept commented out below, but this now just reformats the inputs and passes them through to the new function

%Corwin Wright, 2022/09/23, c.wright@bath.ac.uk


%call new function
[VarOut,PWStore] = pwfilter(NPWs,0, ...
                            SatelliteData.Lon,SatelliteData.Lat,SatelliteData.Temp, ...
                            -180:20:180,-90:PWLat:90, ...
                            SatelliteData.Height,SatelliteData.Height(:,1));

%reformat outputs
Perturbations = VarOut;
PWs.PWs = squeeze(nanmean(PWStore,2));

%done



% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %remove PWs, by fitting sine waves to lat bands
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % 
% % % % % %create array to store PW amplitudes in, with geolocation
% % % % % PWs.z  = SatelliteData.Height(:,1);
% % % % % PWs.Lat = -90:PWLat:90; PWs.Lat =PWs.Lat+PWLat./2;
% % % % % PWs.Modes = 0:1:NPWs;
% % % % % PWs.PWs = NaN(numel(PWs.Lat),numel(PWs.z),numel(PWs.Modes));
% % % % % 
% % % % % %copy of raw T, for background
% % % % % SatelliteData.BGTemp = SatelliteData.Temp;
% % % % % 
% % % % % %loop over height
% % % % % % if Verbose == 1; textprogressbar('Computing planetary waves     '); end
% % % % % for iLevel=1:1:numel(SatelliteData.Height(:,1));
% % % % % %   if Verbose == 1; textprogressbar(iLevel./numel(SatelliteData.Height(:,1)).*100); end
% % % % %   
% % % % %   PWCalc.Lat  = SatelliteData.Lat(       :);
% % % % %   PWCalc.Lon  = SatelliteData.Lon(       :);
% % % % %   PWCalc.Temp = SatelliteData.Temp(      iLevel,:)';
% % % % %   
% % % % %   
% % % % %   %planetary wave removal
% % % % %   for iLat=1:1:ceil(180./PWLat);
% % % % %     
% % % % %     %find boundaries of lat range
% % % % %     LatBinStart = -90+(iLat-1).*PWLat;
% % % % %     LatBinEnd   = -90+(iLat  ).*PWLat;
% % % % %     
% % % % %     %skip if we don't need this bin
% % % % %     if LatBinStart > max(LatRange); continue; end
% % % % %     if LatBinEnd   < min(LatRange); continue; end    
% % % % %     
% % % % %     %find profiles in bin
% % % % %     InLatBin = find(PWCalc.Lat >= LatBinStart & PWCalc.Lat < LatBinEnd);
% % % % %     if numel(InLatBin) == 0; clear InLatBin LatBinStart LatBinEnd; continue; end;
% % % % %     
% % % % %     TheLons  = PWCalc.Lon( InLatBin);
% % % % %     TheTemps = PWCalc.Temp(InLatBin);
% % % % %     
% % % % %     %remove NaNs
% % % % %     Valid = find(~isnan(TheLons + TheTemps));
% % % % %     TheLons  = TheLons( Valid);
% % % % %     TheTemps = TheTemps(Valid);
% % % % % 
% % % % %     
% % % % %     if numel(Valid) == 0; continue; end
% % % % %     
% % % % %     %sort by lon
% % % % %     [~,idx] = sort(TheLons,'ascend');
% % % % %     TheLons  = TheLons(idx);
% % % % %     TheTemps = TheTemps(idx);
% % % % %     
% % % % %     %check we have enough points for a meaningful answer
% % % % %     %twice the number needed to nyquist-sample the smallest requested wave
% % % % %     if numel(TheLons) < 2.*2.*NPWs; 
% % % % %       %remvoe the satellite data - it's not meaningful
% % % % %       SatelliteData.Temp(iLevel,(InLatBin(Valid(idx)))) = NaN;
% % % % %       %and loop
% % % % %       continue;
% % % % %     end
% % % % % 
% % % % %     %find and store zonal mean
% % % % %     PWs.PWs(iLat,iLevel,1) = nanmean(TheTemps);
% % % % %     
% % % % %     %fit planetary wave modes
% % % % %     for iPW=1:1:NPWs;
% % % % %       %identify frequency of signal we're looking for
% % % % %       WaveFreq = 1./(360./iPW); %per degree
% % % % %       
% % % % %       %compute the wave using sinefit (IEEE-1057 standard fitting)
% % % % %       %flags: verbose,plot,iterate to improve
% % % % %       
% % % % %       [Params,Wave] = sinefit(TheTemps,TheLons,WaveFreq,0,0,0);
% % % % %       %remove PW from signal
% % % % %       TheTemps = TheTemps - Wave; 
% % % % %       %store the amplitude - might be useful...
% % % % %       PWs.PWs(iLat,iLevel,iPW+1) = Params(2);
% % % % % 
% % % % %       clear Wave WaveFreq
% % % % %     end; clear iPW
% % % % %     
% % % % %     %put the temperatures back where they came from
% % % % %     SatelliteData.Temp(iLevel,(InLatBin(Valid(idx)))) = TheTemps;
% % % % %     
% % % % %     clear LatBinStart LatBinEnd TheLons TheLats TheTemps idx Valid
% % % % %   end; clear iLat
% % % % %   
% % % % %   clear PWCalc
% % % % % end; clear iLevel
% % % % % % if Verbose == 1; textprogressbar(' '); end
% % % % % 
% % % % % %remove any remaining outlier T (e.g. only one point in the latbin will not fit successfully)
% % % % % SatelliteData.Temp(abs(SatelliteData.Temp) > 100) = NaN;
% % % % % 
% % % % % %compute the background
% % % % % BGTemp = SatelliteData.BGTemp - SatelliteData.Temp;
% % % % % Perturbations =  SatelliteData.Temp;
% % % % % 
% % % % % end
