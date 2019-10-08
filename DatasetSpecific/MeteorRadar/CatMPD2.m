function [Success] = CatMPD2(DataDir,OutputDir,Site,Period)

% %testing parameters
% DataDir   = 'C:\Data\WorkData\MeteorRadar\';
% OutputDir = './';%'C:\Data\WorkData\MeteorRadar\Merged\';
% Site      = 'rothera-sk';
% 
% Period.StartTime   = datenum(2005,3,01,00,00,00);
% Period.EndTime     = datenum(2005,3,31,23,59,59);
% Period.CatDuration = 30.;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%merge meteor radar data into monthly files
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%16/JAN/2014
%
%derived from original code sourced from Robin Davis
%upgraded to handle arbitrary time periods rather than just calendar months
%does not check for previous file existence - automatically clobbers
%
%inputs
%---------
%
%DataDir   - directory containing the daily files to read
%OutputDir - directory to write the monthly files to 
%Site      - meteor radar site ID (e.g. 'rothera-sk' for rothera)
%Period    - structure containing (all in matlab time format)
%            - Period.StartTime
%            - Period.EndTime
%            - Period.CatDuration (optional) - time per chunk (if not specified, uses calendar months)
%outputs
%---------
%
%a merged file for each analysed motnh, stored in the output directory
%format:
%           column 1:  Digital time
%           column 2:  hour
%           column 3:  Height
%           column 4:  Zenith angle (theta)
%           column 5:  Azimuth (phi)
%           column 6:  Horizontal drift velocity (V Horiz)
%           column 7:  Peak amplitude of meteor signal (amax)
%           column 8:  Multiple meteor entry (0 or 1)
%           column 9:  Range
%           column 10: Ambig
%           column 11: Tau
%           column 12: snrdb (Signal to noise ratio in db)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%work out how many periods we want to analyse, and how long they each are individually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('GeneralFunctions') % add path where meteor radar functions are (AM)

disp('=====================================================');
disp('==============MERGING METEOR RADAR DATA==============');
disp('=====================================================');


if ~isfield(Period,'CatDuration');

  %set up to use calendar months.
  %backwards-compatible with Robin's method - will generate same filename format
  %includes the whole of the first and last motnh, regardless of days etc specified in input
  
  %work out the number of months between the start and end date
  [Period.Start.Year, Period.Start.Month,Period.Start.Day,~,~,~] = datevec(Period.StartTime);
  [Period.End.Year  , Period.End.Month  ,Period.End.Day  ,~,~,~] = datevec(Period.EndTime);  
  
  NPeriods             = 1                                         ...
                       + 12.*(Period.End.Year -Period.Start.Year ) ...
                       +     (Period.End.Month-Period.Start.Month);
  SubperiodDefinitions = zeros(2,NPeriods); %[start,end] of each defined period
  
  for iPeriod=1:1:NPeriods;
    
    CalendarMonth = 1+mod(iPeriod+Period.Start.Month-2,12);
    CalendarYear  = Period.Start.Year + ceil((1+iPeriod+Period.Start.Month-2)/12)-1;
    DaysInMonth   = cjw_nmonthdays(CalendarMonth,CalendarYear);
   
    SubperiodDefinitions(1,iPeriod) = datenum(CalendarYear,CalendarMonth,         01,00,00,00);
    SubperiodDefinitions(2,iPeriod) = datenum(CalendarYear,CalendarMonth,DaysInMonth,23,59,59);

  end
  
  Period = rmfield(Period,'Start');
  Period = rmfield(Period,'End');

else

  %evenly-spaced chunks of time as specified. much easier! 
  %filename format based on first time in each file
  
  TotalTime = Period.EndTime - Period.StartTime;
  
  NPeriods = ceil(TotalTime/Period.CatDuration);
  SubperiodDefinitions = zeros(2,NPeriods); %[start,end] of each defined period
  for iPeriod=1:1:NPeriods;
    SubperiodDefinitions(1,iPeriod) = Period.StartTime+(iPeriod-1).*Period.CatDuration;
    SubperiodDefinitions(2,iPeriod) = Period.StartTime+(iPeriod  ).*Period.CatDuration-1.15740112960339e-05; %this is one second in matlab time
    
    %make sure we haven't gone over the end of the period we want to analyse
    if SubperiodDefinitions(2,iPeriod) > Period.EndTime ; SubperiodDefinitions(2,iPeriod) = Period.EndTime;end;
  end
end
  
%tidy up
clear iPeriod CalendarMonth CalendarYear DaysInMonth TotalTime;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%process the data into the now-defined periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iPeriod=1:1:NPeriods; %(months by default)
  counter(iPeriod,NPeriods,'Processing periods ',12);
  
  %create storage array for this period
  Data = [];
  
  %work out days between start and end of period, and loop over them individually
  NDays = ceil(SubperiodDefinitions(2,iPeriod)-SubperiodDefinitions(1,iPeriod));
  for iDay=1:1:NDays;
    
    %identify the calendar day
    DayOfInterest = SubperiodDefinitions(1,iPeriod)+iDay-1;
    
    %locate file containing this day
    FileString = ['/mp', ...
                  num2str(datestr(DayOfInterest,'yyyy')), ...
                  num2str(datestr(DayOfInterest,'mm')),   ...
                  num2str(datestr(DayOfInterest,'dd')),   ...
                  '.',Site,'.mpd'];
    FilePath = [DataDir,Site,FileString];
    
    if exist(FilePath,'var') ~= 0;
      %we have a file for this day: extract the data
      FileData = ReadMPD(FilePath);
      

      if length(FileData) >= 1;
        %the file has contents! convert units
        [Year,Month,Day,~,~,~] = datevec(DayOfInterest);
        FileData(:,1) = FileData(:,1)/24+datenum(Year,Month,Day);
        
        %append to data-so-far
        Data = ([Data;FileData]);
        
      end
      
      %clean up
      clear FileData Year Month Day FilePath FileString DayOfInterest
   
    end
  end; clear iDay;
  
 %save data for the period (if we have any)
 if numel(Data) ~= 0;
   
   StartTime = num2str(SubperiodDefinitions(1,iPeriod));
   EndTime   = num2str(SubperiodDefinitions(2,iPeriod));
   
   if ~isfield(Period,'CatDuration');
     %monthly data output - use Robin's file naming format
     [Year,Month,~,~,~,~] = datevec(SubperiodDefinitions(1,iPeriod));
     OutFile = [OutputDir,                                     ...
       Site,                                          ...
       datestr(datenum(Year,Month,15),'mmm'),         ...
       num2str(datestr(datenum(Year,Month,15),'yyyy'))...
       '.mat'];
     clear Year Month;
   else
     %fixed-period data - generate filename based on dates stored in file
     OutFile = [OutputDir,Site,StartTime,'_',EndTime,'.mat'];
   end
   save(OutFile,'Data','StartTime','EndTime')
 end
end; clear iPeriod;


disp('=====================================================');
disp('=========================DONE========================');
disp('=====================================================');

Success = 1;