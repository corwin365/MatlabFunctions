% set function directories
SystemName = upper(char(java.net.InetAddress.getLocalHost.getHostName));
if strcmp(SystemName,'BETTARAVE') == 1
  %Dell laptop
  addpath(genpath('C:\Users\cw785-admin\Dropbox\Code\MatlabFunctions'));
elseif strcmp(SystemName,'PASTEQUE') == 1
  %Lenovo desktop
  addpath(genpath('C:\Users\cw785-admin\Dropbox\Code\MatlabFunctions'));
elseif strcmp(SystemName,'EEPC-0184') == 1
  %eepc-0184.bath.ac.uk
  addpath(genpath(['/u/f/cw785/Code/MatlabFunctions']))
  setenv('LD_LIBRARY_PATH', ['/usr/lib/x86_64-linux-gnu:',getenv('LD_LIBRARY_PATH')]);  %fixes a pathing problem on the Bath system when making system() calls
elseif isunix
  %unix system (assumes on Bath network)
  addpath(genpath(['/u/f/cw785/Code/MatlabFunctions']))
else
  %unspecified windows machine
  addpath(genpath([pwd,'\MatlabFunctions']));
end
clear SystemName


%formatting
format compact       %no gap lines between information
format longG         %number formatting
beep on              %system beep on error


%larger default fonts for plotting and text
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);

%grids on figures by default
set(0, 'defaultAxesXGrid', 'on')
set(0, 'defaultAxesYGrid', 'on')
set(0, 'defaultAxesZGrid', 'on')

%default to a white background on figures
set(0,'defaultfigurecolor',[1 1 1])

%pause on error
dbstop if error


%resume in last used directory
if ispref('StartupDirectory','LastWorkingDirectory')
  lwd = getpref('StartupDirectory','LastWorkingDirectory');
  try
    cd(lwd)
  catch
    disp('Sorry, but I could not go to your last working directory:')
    disp(lwd)
  end
  clear lwd
end

%print where we are
disp('--------------------------------------------------------------------------------')
disp(['Starting in: ',pwd])
disp('--------------------------------------------------------------------------------')
