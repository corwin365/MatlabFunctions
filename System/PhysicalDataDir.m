function DataDir = PhysicalDataDir()

%loads data from EEPC-291's physical disks
%faster runtime, but obviously not portable
%falls back gracefully if on different machine: calls the usual location if not running on eepc-291
%will not fall back safely if on eepc-291 but the data we actually want isn't on the physical drive


persistent PhysicalDriveCalled;
[~,SystemName] = system('hostname');

if numel(SystemName) ==9 && strcmp(SystemName(1:8),'eepc-291') == 1;
  %specify the physical path
  DataDir = '/u/f/cw785/eepc291-main/';
  
  %puts a warning up the first time this is called
  if isempty(PhysicalDriveCalled);

    PhysicalDriveCalled = 1;
    
    disp(' ')
    disp('=================================================================')
    disp('*****************************************************************')
    disp('WARNING --- USING PHYSICAL DRIVE PATH, MAKE SURE DATA IS PRESENT!')
    disp(' THIS WARNING WILL NOT DISPLAY AGAIN UNTIL ALL VARS ARE CLEARED. ')
    disp('               THIS MAY BE YOUR FINAL WARNING!                   ')    
    disp('*****************************************************************')
    disp('=================================================================')
    disp(' ')
  end
  

else
  
  %just call the standard directory path function
  DataDir = LocalDataDir(); %the usual data path routine gets called
end


return
end