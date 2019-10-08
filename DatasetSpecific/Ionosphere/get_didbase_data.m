function [Results,ResultsTable,VarList] = get_didbase_data(Site,StartTime,EndTime,Variables,Verbose,DMUF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get ionosonde data from DIDBase (http://giro.uml.edu/didbase/scaled.php)
%Corwin Wright, c.wright@bath.ac.uk, 25/JAN/2017
%modified 07/MAR/2017 to account for (unidentified) changes to didbase website
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
%
% % %example:
% % clear all
% %  StartTime   = datenum(2012,7,2,12,0,0);
% %  EndTime     = datenum(2012,7,3,12,0,0);
% %  Site        = 'MHJ45';
% %  Verbose     = 1;
% %  Variables   = {'foF2','hmF2'};
% %  DMUF        = 3000;
% %  [Results,ResultsTable] = get_didbase_data(Site,StartTime,EndTime,Variables,Verbose,DMUF);
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
if nargin < 4;
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

if Verbose == 1;
  disp('--------------------------------------------------------------');
  disp(['Getting data from site ',Meta.Site,' between ',Meta.StartDate,' and ',Meta.EndDate])
  disp('--------------------------------------------------------------');
end

%allow a longer timeout, needed sometimes
options = weboptions('timeout',60);

%loop over variables
for iVar=1:1:numel(Variables);
  if Verbose == 1;  disp(['Querying DIDbase for ',Variables{iVar}]); end
  
% %   %to avoid hammering the server, impose a random delay on subsequent requests
% %   if iVar > 1; pause(0+rand(1)*3); end %0-3 seconds
  
  %id the variable
  Variable = Variables{iVar};
  
  %get the data from the website for this variable
  URL = strcat(Meta.URL,'?ursiCode=',Meta.Site, ...
                        '&charName=',Variable,  ...
                        '&DMUF=',num2str(Meta.DMUF), ...
                        '&fromDate=',Meta.StartDate, ...
                        '&toDate=',Meta.EndDate);
  WebPage = webread(URL,options);                
  
                  
  %% split the webpage into lines
  %use a regexp on '\n' for compatibility with old Matlab versions
  LineEnds = regexp(WebPage,'\n')+1; % the +1 is the first char after '\n'
  LineEnds(end) = LineEnds(end)-1; %last line ends at end of file
  
  WebPage2 = cell(numel(LineEnds),1);
  for iLine=1:1:numel(LineEnds)-1;
    WebPage2{iLine} = WebPage(LineEnds(iLine):LineEnds(iLine+1)-1);
  end; clear iLine LineEnds
  WebPage = WebPage2; clear WebPage2 

  %% work out which lines have data, and store them
  for iLine=1:1:numel(WebPage);
    Line = WebPage{iLine};
    if numel(Line) == 0;          continue; end; %empty line
    if strcmp(Line(1),'#');       continue; end; %no data on this line
    if strcmp(Line(1:5),'ERROR'); continue; end; %problem with didbase call
    
    
    %we have a data line
    %this contains four columns: [time,CS,variable,QD]
    %unsure what QD is, is NaN in every example I've looked at, so ignore it
    
    %split the line into variables, and remove spaces
    Splits = [0,regexp(Line,'[ ]+'),numel(Line)];
    Row = {};
    for iCol=1:1:numel(Splits)-1;
      Row{end+1} = strtrim(Line(Splits(iCol)+1:Splits(iCol+1)));
    end; 
    clear iCol Splits Line
    
    %convert to numbers 
    Time = datenum(Row{1},'yyyy-mm-ddTHH:MM:SS.FFFZ');
    CS   = str2double(Row{2});
    Var  = str2double(Row{3});
        
    AllData(:,end+1) = [iVar,Time,CS,Var]; %iVar needed later to de-dupe lines
    clear Time CS Var Row
  end;
  clear  Variable WebPage iLine Line
end; clear iVar


if Verbose == 1; 
  if numel(AllData) == 0;
    disp('No data obtained - check inputs');
    Results = [];
    ResultsTable = [];
    VarList = [];
    return
  else
    disp('Data obtained, post-processing'); 
  end
end


%% produce table of results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first, discard anything outside our exact times of interest
Good = find(AllData(2,:) >= StartTime ...
          & AllData(2,:) <= EndTime);
AllData = AllData(:,Good);
clear Good StartTime EndTime

%%ok. create a table with a unique time axis
[UniqueTimes,~,idx] = unique(AllData(2,:));
ResultsTable = NaN(numel(UniqueTimes),numel(Variables)+2);
ResultsTable(:,1) = UniqueTimes;

%and fill it
for iTime=1:1:numel(UniqueTimes);
  
  %find all lines in AllData at this time
  AtTime = AllData(:,idx == iTime);

  %CS is the same for all vars, so just store the first value
  ResultsTable(iTime,2) = AtTime(3,1);
  
  %put these data into the table
  for iVar=1:1:size(AtTime,2);
    %put data from AllData column iVarinto the correct column of ResultsTable
    %offset of +2 is to leave space for Time and CS
    ResultsTable(iTime,AtTime(1,iVar)+2) = AtTime(4,iVar); 
  end; clear iVar AtTime
  
end; clear iTime UniqueTimes


if Verbose == 1; disp('Results table formatted'); end

%% finally, pull the variables into individual vars (for backwards compat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results.Time = ResultsTable(:,1);
Results.CS   = ResultsTable(:,2);

for iVar=1:1:numel(Variables);
  eval(['Results.',Variables{iVar},' = ResultsTable(:,iVar+2);'])
end; clear iVar 


if Verbose == 1; disp('Results struct produced'); end



%done all the work. Finally, create list of column headings for table
VarList = {'Time','CS'};
for i=1:1:numel(Variables); VarList{end+1} = Variables{i}; end
return
  
end
