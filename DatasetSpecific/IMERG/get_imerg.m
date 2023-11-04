function Out = get_imerg(DateList,Path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load IMERG data in HDF5 format into a Matlab struct
%
%inputs:
%  DateList - list of Matlab dates to load. Only works at scale of >= 1 day, and will return on a 30min scale whatever you request. Non-continuous periods will include gaps.
%  Path (optional) -  path to IMERG data, in HDF5 format. Will search subdirectories of this.
%
%outputs:
%  Out - struct containing:
%    Lat     - latitude  of box centre, deg N
%    Lon     - longitude of box centre, deg E
%    Precip  - precipitation in box, mm/hr
%    Error   - random error, mm/hr
%    PLiquid - probability of precipitation being liquid
%
%Corwin Wright, c.wright@bath.ac.uk
%2023/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.TimeScale = DateList; clear DateList
if nargin < 2; 
  Settings.InDir = [LocalDataDir,'/IMERG/dl/'];
else           
  Settings.InDir = Path;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.TimeScale);


  %% get and store data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %get files for this day
  [y,m,d] = datevec(Settings.TimeScale(iDay));
  FileString = ['*IMERG.',sprintf('%04d',y),sprintf('%02d',m),sprintf('%02d',d),'-*.hdf5'];
  Files = wildcardsearch([LocalDataDir,'/IMERG/dl/'],FileString);

  if numel(Files) == 0; clear y m d FileString Files; continue; end

  textprogressbar(['Loading IMERG for ',datestr(Settings.TimeScale(iDay)),' '])

  %loop over the files
  for iFile=1:1:numel(Files)

    %load file
    Data = read_imerg(Files{iFile});


    %create storage arrays
    if ~exist('Out','var')
      Out = struct();
      Out.Time     = linspace(min(Settings.TimeScale),max(Settings.TimeScale)+1,48.*(range(Settings.TimeScale)+1)+1); Out.Time = Out.Time(1:end-1);
      Out.Lat      = Data.lat;
      Out.Lon      = Data.lon;
      Out.Precip   = NaN(numel(Out.Lat),numel(Out.Lon),48.*numel(Settings.TimeScale),'single');
      Out.Error    = Out.Precip;
      Out.PLiquid  = Out.Precip;
    end

    %store the data
    idx = closest(Data.time,Out.Time);
    Out.Precip( :,:,idx) = Data.precipitation;
    Out.Error(  :,:,idx) = Data.randomError;
    Out.PLiquid(:,:,idx) = Data.probabilityLiquidPrecipitation;

    clear Data idx

    textprogressbar(iFile./numel(Files).*100)
  end

  textprogressbar(100); textprogressbar('!')

end
