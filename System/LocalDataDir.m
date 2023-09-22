function LocalDataDir = LocalDataDir()

%identifies which machine we're using, and sets the path to the data store directory accordingly
TheComputerThisIsOn = upper(char(java.net.InetAddress.getLocalHost.getHostName));

%identify the user as wll, for some systems
TheUser = char(java.lang.System.getProperty('user.name'));


%set data directory path
if strcmp(TheComputerThisIsOn,'BETTERAVE'); %Corwin's 2023 XPS13
  LocalDataDir = 'C:\Data\';
elseif strcmpi(TheComputerThisIsOn,'neils-macbook-pro')
  LocalDataDir = '/Users/neil/data/';
elseif isunix %assume Bath Uni filesystem
  if     strcmp(TheUser,'cw785'); LocalDataDir = '/u/f/cw785/Data/';
  elseif strcmp(TheUser,'nh351'); LocalDataDir = '/u/f/nh351/Data/';    
  else
    LocalDataDir = '/';
    warning('Unix system for whom LocalDataDir is not configured')
  end
else %assume Bath Uni filesystem
  warning('Non-unix system for whom LocalDataDir is not configured')
  LocalDataDir = 'C:\Data\';
end


return
