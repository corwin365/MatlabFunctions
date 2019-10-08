function LocalDataDir = LocalDataDir()

%identifies which machine we're using, and sets the path to the data store directory accordingly
TheComputerThisIsOn = upper(char(java.net.InetAddress.getLocalHost.getHostName));

%set data directory path
if strcmp(TheComputerThisIsOn,'DESKTOP-RNCK39U') == 1
  %Dell laptop
  LocalDataDir = 'C:\Data\';
elseif strcmp(TheComputerThisIsOn(1:4),'MYRT') == 1 
  %Razer laptop - Myrtille
 LocalDataDir = 'D:\Data\';  
elseif strcmp(TheComputerThisIsOn(1:4),'ITD-') == 1 || strcmp(TheComputerThisIsOn(1:5),'NODE-') == 1
  %Balena
  LocalDataDir = '/home/f/cw785/scratch/Data/';
elseif isunix 
  %unix system (assumes on Bath network)
  LocalDataDir = '/u/f/cw785/Data/';
else
  %unspecified windows machine (assumes on Bath network)
  LocalDataDir = 'Z:\Data\';
end


return
