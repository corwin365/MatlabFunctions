function [Results,ResultsTable,VarList] = get_didbase_data_v2(Site,StartTime,EndTime,Variables,Verbose,DMUF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get ionosonde data from DIDBase (http://giro.uml.edu/didbase/scaled.php)
%Corwin Wright, c.wright@bath.ac.uk, 25/JAN/2017
%modified 07/MAR/2017 to account for (unidentified) changes to didbase website
%modified 14/MAY/2017 to improve speed by factor of ~100
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%required inputs:
% - Site: location ID, listed at http://giro.uml.edu/didbase/scaled.php
% - StartTime: start time, Matlab format
% - EndTime: end time, Matlab format
%
%StartTime and EndTime are precise, i.e. to get a single whole day set 
%EndTime to the beginning of the next day
%
%optional inputs: 
% - Verbose: 1 to display progress, any other value to not
% - Variables: cell array of requested variables
%            - requesting a large number of vars can be very slow
%            - defaults to all variables if not specified - slowest option
% - DMUF: maximum usable frequency for ground distance
%       - defaults to 3000km
% %
% % % %example:
% clear all
%  StartTime   = datenum(2012,7,2,12,0,0);
%  EndTime     = datenum(2012,7,3,12,0,0);
%  Site        = 'MHJ45';
%  Verbose     = 1;
%  Variables   = {'foF2','hmF2','fbEs'};
%  DMUF        = 3000;
%  [Results,ResultsTable] = get_didbase_data(Site,StartTime,EndTime,Variables,Verbose,DMUF);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%outputs:
% - Results - structure containing either all the following or a subset 
%             specified by use of the 'Variables' input
% % % Time: time, in Matlab units [[always included, no need to specify]]
% % % CS is Autoscaling Confidence Score (from 0 to 100, 999 if manual scaling, -1 if unknown)
% % % foF2 [MHz] - F2 layer critical frequency
% % % foF1 [MHz] - F1 layer critical frequency
% % % MD [] - MUF(D)/foF2
% % % MUFD [MHz] - Maximum usable frequency for ground distance D
% % % fmin [MHz] - Minimum frequency of ionogram echoes
% % % foEs [MHz] - Es layer critical frequency
% % % fminF [MHz] - Minimum frequency of F-layer echoes
% % % fminE [MHz] - Minimum frequency of E-layer echoes
% % % foE [MHz] - E layer critical frequency
% % % fxI [MHz] - Maximum frequency of F trace
% % % hF [km] - Minimum virtual height of F trace
% % % hF2 [km] - Minimum virtual height of F2 trace
% % % hE [km] - Minimum virtual height of E trace
% % % hEs [km] - Minimum virtual height of Es trace
% % % hmE [km] - Peak height of E-layer
% % % yE [km] - Half thickness of E-layer
% % % QF [km] - Average range spread of F-layer
% % % QE [km] - Average range spread of E-layer
% % % FF [MHz] - Frequence spread between fxF2 and fxI
% % % FE [MHz] - Frequence spread beyond foE
% % % hmF2 [km] - Peak height F2-layer
% % % hmF1 [km] - Peak height F1-layer
% % % zhalfNm [km] - The true height at half the maximum density in the F2-layer
% % % fminEs [MHz] - Minimum frequency of Es-layer
% % % yF2 [km] - Half thickness of F2-layer, parabolic model
% % % yF1 [km] - Half thickness of F1-layer, parabolic model
% % % TEC [10^16 m^-2] - Total electron content
% % % scaleF2 [km] - Scale hieght at the F2-peak
% % % B0 [km] - IRI thickness parameter
% % % B1 [] - IRI profile shape parameter
% % % D1 [] - IRI profile shape parameter, F1-layer
% % % foEa [MHz] - Critical frequency of auroral E-layer
% % % hEa [km] - Minimum virtual height of auroral E-layer trace
% % % foP [MHz] - Highest ordinary wave critical frequency of F region patch trace
% % % hP [km] - Minimum virtual height of the trace used to determinate foP
% % % fbEs [MHz] - Blanketing frequency of Es-layer
% % % TypeEs [] - Type Es
% % %
% % % - ResultsTable - table ntimes x (nvars+2)
% % %                   - columns are [Time,CS,Vars in order requested]
% % % 
% % % - VarList - column headings for ResultsTable, as a cross-check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% variable pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%which variables to we want?
if nargin < 4
  %if none specified, just grab all the variables DIDBASE has
  %downloading all is pretty fast nominally, but I impose a 2s wait between each
  %variable in order to avoid pumelling the server
  Variables = {'foF2','foF1','foE','foEs','fbEs','foEa','foP','fxI','MUFD', ...
               'MD','hF2','hF','hE','hEs','hEa','hP','TypeEs','hmF2','hmF1', ...
               'hmE','zhalfNm','yF2','yF1','yE','scaleF2','B0','B1','D1',...
               'TEC','FF','FE','QF','QE','fmin','fminF','fminE','fminEs'};
end
%be talkative or quiet?
if nargin < 5; Verbose = 1; end %talkative by default
%did we specify DMUF?
if nargin < 6;DMUF = 3000; end

%rename some variables, for internal tidiness
Meta.DMUF = DMUF; clear DMUF;
Meta.Site = Site; clear Site;

%% other preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify the DIDBASE URL
Meta.URL = 'http://lgdc.uml.edu/common/DIDBGetValues';

%convert the time inputs from Matlab time to the time format used by DIDBASE
%working at the whole day level here as the time format is a bit awkward
%later in the routine we'll discard unwanted parts of the first and last day
[y,m,d] = datevec(StartTime);
Meta.StartDate = [sprintf('%04d',y),'.',sprintf('%02d',m),'.',sprintf('%02d',d)];%, ...
[y,m,d] = datevec(EndTime+1); %the +1 is so we get the whole of the last day
Meta.EndDate   = [sprintf('%04d',y),'.',sprintf('%02d',m),'.',sprintf('%02d',d)];
clear y m d


%% get the data for all requested variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AllData = []; %for later use

if Verbose == 1
  disp('--------------------------------------------------------------');
  disp(['Getting data from site ',Meta.Site,' between ',Meta.StartDate,' and ',Meta.EndDate])
  disp('--------------------------------------------------------------');
end

%allow a longer timeout, needed sometimes
options = weboptions('timeout',60);

%loop over variables
  
%get the data from the website for this variable
URL = strcat(Meta.URL,'?ursiCode=',Meta.Site, ...
                      '&charName=',strjoin(Variables,','),  ...
                      '&DMUF=',num2str(Meta.DMUF), ...
                      '&fromDate=',Meta.StartDate, ...
                      '&toDate=',Meta.EndDate);
WebPage = webread(URL,options);


%% split the webpage into lines
%use a regexp on '\n' for compatibility with old Matlab versions
LineEnds = regexp(WebPage,'\n')+1; % the +1 is the first char after '\n'
LineEnds(end) = LineEnds(end)-1; %last line ends at end of file

WebPage2 = cell(numel(LineEnds),1);
for iLine=1:1:numel(LineEnds)-1
  WebPage2{iLine} = WebPage(LineEnds(iLine):LineEnds(iLine+1)-1);
end; clear iLine LineEnds
WebPage = WebPage2; clear WebPage2

%% work out which lines have data, and store them
for iLine=1:1:numel(WebPage)
  Line = WebPage{iLine};
  if numel(Line) <= 0;          continue; end %empty line
  if strcmp(Line(1),'#')        %no data on this line - check if it's a header
    if numel(Line) > 5
      if strcmp(Line(1:5),'#Time'); Headers = Line; end %header row
    end
    continue 
  end
  if strcmp(Line(1:5),'ERROR'); continue; end %problem with didbase call
  
  
  %we have a data line
  %this contains columns: [time,CS,variable,QD]
  %unsure what QD is, is NaN in every example I've looked at, so ignore it
  
  %split the line into variables, and remove spaces
  Splits = [0,regexp(Line,'[ ]+'),numel(Line)];
  Row = {};
  for iCol=1:1:numel(Splits)-1
    Row{end+1} = strtrim(Line(Splits(iCol)+1:Splits(iCol+1)));
  end
  clear iCol Splits Line
    
  %convert to numbers
  NumRow = NaN(size(Row));
  NumRow(1) = datenum(Row{1},'yyyy-mm-ddTHH:MM:SS.FFFZ');
  for iCol=2:1:numel(NumRow)
    NumRow(iCol) = str2double(Row{iCol});
  end
  
  AllData(end+1,:) = NumRow;
  clear CS Var Row NumRow
end
clear  Variable WebPage iLine Line

if numel(AllData) == 0
  if Verbose == 1; disp('No data obtained - check inputs'); end
  Results = [];
  ResultsTable = [];
  VarList = [];
  return
elseif Verbose == 1
  disp('Data obtained, post-processing');
end

%split up header row
HeadSplits = [0,regexp(Headers,'[ ]+'),numel(Headers)];
Headers2 = {};
for iCol=1:1:numel(HeadSplits)-1
  Headers2{end+1} = strtrim(Headers(HeadSplits(iCol)+1:HeadSplits(iCol+1)));
  %make it clear which QC, etc columns are which
  if strcmp(Headers2{end},'QD') == 1 
    Headers2{end} = [Headers2{end-1},'_QD'];
  end
end
Headers2{1} = 'Time'; %remove the # at start
Headers = Headers2;
clear HeadSplits Headers2 iCol options URL Variables

%% produce table of results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first, discard anything outside our exact times of interest
Good = find(AllData(:,1) >= StartTime ...
          & AllData(:,1) <= EndTime);
ResultsTable = AllData(Good,:);
clear Good StartTime EndTime AllData

if Verbose == 1; disp('Results table prepared'); end

%% finally, pull the variables into individual vars (for backwards compat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results.Time = ResultsTable(:,1);
Results.CS   = ResultsTable(:,2);

for iVar=3:1:numel(Headers)
  eval(['Results.',Headers{iVar},' = ResultsTable(:,iVar);'])
end; clear iVar 


if Verbose == 1; disp('Results struct produced'); end



%done all the work. Finally, create list of column headings for table
VarList = Headers;
  
end
