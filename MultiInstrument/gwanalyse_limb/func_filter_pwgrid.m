%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% planetary wave filter, grid-based approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data,PW] = func_filter_pwgrid(Data,Settings)


%compute a time grid to work on, and also store metadata about the PWs
PW.Time = min(Settings.TimeRange):Settings.PWTimeRes:max(Settings.TimeRange);
PW.Lon = Settings.PWLonGrid; PW.Lat = Settings.PWLatGrid; PW.Alt = Settings.PWAltGrid; PW.PWs = 0:1:Settings.NPWs;
PW.WindowSize = Settings.PWWindow; PW.MinPercent = Settings.PWMinPC;

%generate a store array for the PWs, and for the output T'
PW.PW = NaN(numel(Settings.PWLonGrid),numel(Settings.PWLatGrid),numel(Settings.PWAltGrid),Settings.NPWs+1,numel(PW.Time));
A = Data.Tp.*NaN; %working variable used internally to simplify logic

%fill it, stepping over day-by-day using a time window as specified
if Settings.Verbose == 1; textprogressbar('--> Computing planetary waves '); end

for iStep=1:1:numel(PW.Time)

  %select the data we need by finding the indices of the points in the UseWindow (points to compute from) and
  %the Output window (points to store, higher resolution)
  %***logic assumes idxO is a subset of idxU***
  idxU = inrange(Data.Time,PW.Time(iStep)+[-1,1].*(Settings.PWWindow)./2);
  idxO = inrange(Data.Time,PW.Time(iStep)+[-1,1].*(mean(diff(PW.Time)))./2);
  if numel(idxU) == 0 || numel(idxO) == 0; continue; end

  %reduce data down to just what we want to use
  PWCalcData = reduce_struct(Data,idxU,{'OriginalFiles','Note'},0);

  %compute the PWs in the use window, and store it in placeholder A
  %A will overwrite most loops - this is fine as long as idxO is a subset of idxU
  [A(idxU),b] = pwfilter(Settings.NPWs,Settings.PWMinPC,                      ...
                         PWCalcData.Lon,PWCalcData.Lat,PWCalcData.Tp,   ...
                         Settings.PWLonGrid,Settings.PWLatGrid,               ...
                         PWCalcData.Alt,Settings.PWAltGrid);
  %store the Tp and PW data
  Data.Tp(idxO) = A(idxO);
  PW.PW(:,:,:,:,iStep) = permute(b,[2,1,3,4]);

  if Settings.Verbose == 1; textprogressbar(100.*iStep./numel(PW.Time)); end

end; clear iDay OutWindow UseWindow idxU PWCalcData a b idxO iStep
if Settings.Verbose == 1; textprogressbar(100); textprogressbar('!');end

%compute background temperature
I = griddedInterpolant({PW.Lon,PW.Lat,PW.Alt,PW.Time},squeeze(nansum(PW.PW,4)));
Data.BG = I(Data.Lon,Data.Lat,Data.Alt,Data.Time);

return
end
