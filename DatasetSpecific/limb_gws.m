clear all
InData = get_limbsounders(datenum(2007,1,20),'HIRDLS','HeightScale',20:1:60);

OutData = limb_gws_f(InData)

function [OutData,PWs] = limb_gws_f(InData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%general-purpose function to find GWs from limb sounders
%Corwin Wright, c.wright@bath.ac.uk, 2023/09/04
%
%inputs:
%  required:
%    InstrumentData [struct] - data to use, in the format produced by get_limbsounders.m
%  optional:
%
%    VarName         (type,           default)  description
%    -----------------------------------------------------------------------------
%    HeightScale     (numeric,        20:1:60)  heightscale to interpolate the data onto, in km
%
%
%outputs:
%   Data: struct containing all variables, on a [profiles x height] grid
%   PWs:  array containing PWs, on a lon x lat x height x wave x days grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %% input parsing
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % %create parser object
% % % % % %%%%%%%%%%%%%%%%%%%%%%
% % % % % p = inputParser;
% % % % % 
% % % % % %validation of inputs
% % % % % %%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % %date
% % % % % CheckDates  = @(x) validateattributes(x,{'numeric'},{'>=',datenum(1979,1,1)});
% % % % % addRequired(p,'TimeRange',CheckDates); %time range
% % % % % 
% % % % % %instrument
% % % % % CheckInst  = @(x) validateStringParameter(x,fieldnames(InstInfo),mfilename,Instrument);
% % % % % addRequired(p,'Instrument',CheckInst)
% % % % % 
% % % % % %additional variables
% % % % % addParameter(p,'AdditionalVars',{},@iscell)
% % % % % addParameter(p,'OriginalZ',false,@islogical)
% % % % % 
% % % % % %parse inputs and tidy up
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % %parse inputs
% % % % % parse(p,TimeRange,Instrument,varargin{:})
% % % % % 
% % % % % %additional check on dates - must be either a single value (a day to load all of) or two values (start and end time)
% % % % % if     numel(TimeRange)  > 2;  error('TimeRange must be either one value (a day of interest) or two values (a time range)');
% % % % % elseif numel(TimeRange) == 1; TimeRange = [floor(TimeRange),floor(TimeRange)+1-1e-8]; end % strictly the latter bound is 1ms less than a day
% % % % % 
% % % % % %pull out the contents into struct "Inputs", used throughout rest of routine
% % % % % Input = p.Results;
% % % % % Input.TimeRange = TimeRange;
% % % % % clearvars -except InstInfo Input
% % % % % 
% % % % % %extract just the info for the instrument we want
% % % % % InstInfo  = InstInfo.(Input.Instrument);



Input.NPWs     = 6;
Input.MinPC    =  0.66;
Input.LonGrid  = -180:20:180;
Input.LatGrid  = -90:5:90;
Input.AltGrid  = 20:1:60;
Input.PWWindow = 3; %days
Input.TimeRange = [floor(nanmin(InData.Time,[],'all')),ceil(nanmax(InData.Time,[],'all'))];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute planetary waves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate a store array for the PWs
PWStore = NaN(numel(Input.LonGrid),numel(Input.LatGrid),numel(Input.AltGrid),Input.NPWs+1,range(Input.TimeRange)+1);
Tp = InData.Temp.*NaN;

%fill it, stepping over day-by-day using a time window as specified
textprogressbar('Computing planetary waves ')
for iDay=1:1:range(Input.TimeRange)+1
  
  %select the data we need
  TimeWindow = Input.TimeRange(1)+(iDay-1)+[-1,1].*(Input.PWWindow-1)./2;
  idx= inrange(InData.Time,TimeWindow);
  PWCalcData = reduce_struct(InData,idx,[],0);

  [Tp(idx),a] = pwfilter(Input.NPWs,Input.MinPC,                        ...
                         PWCalcData.Lon,PWCalcData.Lat,PWCalcData.Temp, ...
                         Input.LonGrid,Input.LatGrid,                   ...
                         PWCalcData.Alt,Input.AltGrid);
  PWStore(:,:,:,:,iDay) = permute(a,[2,1,3,4]);
stop
  textprogressbar(100.*iDay./(range(Input.TimeRange)+1))

end; clear iDay TimeWindow idx PWCalcData a
textprogressbar(100); textprogressbar('!')


stop
% % [InData.Tp,PWStore] = pwfilter(Input.NPWs,Input.MinPC,            ...
% %                                InData.Lon,InData.Lat,InData.Temp, ...
% %                                Input.LonGrid,Input.LatGrid,       ...
% %                                InData.Alt,Input.AltGrid);



stop


end
