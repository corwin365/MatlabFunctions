function LocalDataDir = LocalDataDir()

%identifies which machine we're using, and sets the path to the data store directory accordingly
[~,TheComputerThisIsOn] = system('hostname');
TheComputerThisIsOn = upper(char(strtrim(TheComputerThisIsOn)));

%identify the user as wll, for some systems
TheUser = char(java.lang.System.getProperty('user.name'));


%set data directory path
if     strcmpi(TheComputerThisIsOn,'BETTERAVE');         LocalDataDir = 'C:\Data\'; %Corwin's red laptop
elseif strcmpi(TheComputerThisIsOn,'PASTEQUE');          LocalDataDir = 'D:\Data\'; %Corwin's legion desktop
elseif strcmpi(TheComputerThisIsOn,'neils-macbook-pro'); LocalDataDir = '/Users/neil/data/'; %Neil's Macbook
elseif isunix 
  if     strcmp(TheUser,'cw785'); LocalDataDir = '/u/f/cw785/Data/';  %Corwin on Bath system
  elseif strcmp(TheUser,'nh351'); LocalDataDir = '/u/f/nh351/Data/';  %Neil on Bath system    
  else
    warning('Unix system for whom LocalDataDir is not configured, guessing path')
    LocalDataDir = '/';    
  end
else
  warning('Non-unix system for whom LocalDataDir is not configured, guessing path')
  LocalDataDir = 'C:\Data\'; %this is a guess
end


return
end