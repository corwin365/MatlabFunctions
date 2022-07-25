function LocalDataDir = LocalDataDir()

%identifies which machine we're using, and sets the path to the data store directory accordingly
TheComputerThisIsOn = upper(char(java.net.InetAddress.getLocalHost.getHostName));

%set data directory path
if strcmp(TheComputerThisIsOn,'DESKTOP-RNCK39U') | strcmp(TheComputerThisIsOn,'BETTERAVE') == 1;
  LocalDataDir = 'C:\Data\';
elseif strcmp(TheComputerThisIsOn,'MYRTILLE') == 1 
 LocalDataDir = 'D:\Data\';  
elseif strcmp(TheComputerThisIsOn(1:4),'ITD-') == 1 || strcmp(TheComputerThisIsOn(1:5),'NODE-') == 1
  LocalDataDir = '/home/f/cw785/scratch/Data/';
elseif isunix 
  LocalDataDir = '/u/f/cw785/Data/';
else
  LocalDataDir = 'Z:\Data\';
end


return
