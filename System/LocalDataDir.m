function LocalDataDir = LocalDataDir()

%identifies which machine we're using, and sets the path to the data store directory accordingly
TheComputerThisIsOn = upper(char(java.net.InetAddress.getLocalHost.getHostName));

%set data directory path
if strcmp(TheComputerThisIsOn,'BETTERAVE'); %Corwin's 2023 XPS13
  LocalDataDir = 'C:\Data\';
elseif strcmpi(TheComputerThisIsOn,'neils-macbook-pro')
  LocalDataDir = '/Users/neil/Data/'
elseif isunix %assume Bath Uni filesystem
  LocalDataDir = '/u/f/cw785/Data/';
else %assume Bath Uni filesystem
  LocalDataDir = 'Z:\Data\';
end


return
