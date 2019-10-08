function [Results] = get_didbase_box(LonRange,LatRange,TimeRange,Variables,Verbose,DMUF)

% % %   %testing parameters
% % %   clear all
% % %   LatRange  = [-1,1].*5+42.60;
% % %   LonRange  = [-1,1].*5+288.50-360;
% % %   TimeRange = [datenum(2016,6,1),datenum(2016,6,7)];
% % %   Verbose   = 1;
% % %   DMUF = 3000;
% % %   Variables = {'foF2','foE'};

  %first, find all stations in the box
  [Stations] = get_didbase_geo(LonRange,LatRange,TimeRange);
  
  %organise variables
  StartTime = min(TimeRange);
  EndTime   = max(TimeRange);
  
  %now, loop over the stations and get their data
  for iStation=1:1:size(Stations,1);
    Site = Stations{iStation,1};
    
    %get the data, handling optional inputs
    switch nargin
      case 3; [~,Data,Rows] = get_didbase_data(Site,StartTime,EndTime);
      case 4; [~,Data,Rows] = get_didbase_data(Site,StartTime,EndTime,Variables);
      case 5; [~,Data,Rows] = get_didbase_data(Site,StartTime,EndTime,Variables,Verbose);
      case 6; [~,Data,Rows] = get_didbase_data(Site,StartTime,EndTime,Variables,Verbose,DMUF);
      otherwise; disp('Geolocation not specified correctly');
    end
    
    if numel(Data) ==0; %no data from station
      disp(['No data available from ',Site])
    else
      
      if ~exist('Results');  end
      
      %store the data
      Station = ones(size(Data,1),1).*iStation;
      Data = cat(2,Station,Data);
      
      if ~exist('AllData');
        AllData = Data;
        Results.ColHeads = [{'Stations'} Rows{1:end} ];
      else AllData = cat(1,AllData,Data);
      end
    end
    
  end
  
  
  %format output
  Results.Data = AllData;
  Results.Stations = Stations;
  
  return