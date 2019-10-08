
% set function directories
SystemName = upper(char(java.net.InetAddress.getLocalHost.getHostName));


if strcmp(SystemName,'DESKTOP-RNCK39U') == 1
  %Dell laptop
  addpath(genpath('C:\Users\cw785-admin\Dropbox\Code\MatlabFunctions'));
elseif strcmp(SystemName(1:4),'MYRT') == 1 
  %Myrtille (my Razer laptop)
  addpath(genpath([pwd,'\MatlabFunctions']));
elseif strcmp(SystemName(1:4),'ITD-') == 1 || strcmp(SystemName(1:4),'NODE-') == 1  
  %Balena  node
  addpath(genpath('/home/f/cw785/Matlab/matlabfunctions'));
elseif isunix 
  %unix system (assumes on Bath network)
  addpath(genpath(['/u/f/cw785/Code/matlabfunctions']))
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


%better-positioned figure controls
try
  if ~verLessThan('matlab','9.5')
    set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
    set(groot,'defaultAxesCreateFcn',  @(ax,~)set(ax.Toolbar,'Visible','off'));
  end
end

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