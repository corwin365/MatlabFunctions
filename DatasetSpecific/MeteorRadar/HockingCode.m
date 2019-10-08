%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Important Bits to edit before running %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters.Site       = 'ascension';
Parameters.MonthlyDir = strcat(cd,'/Data/MthCat/',Parameters.Site,'_mthcat/'); %where the monthly meteor files are stored
Parameters.OutFile    = strcat('HockingOutputFiles/',Parameters.Site,'_','test.mat'); %file to write results to
delete(Parameters.OutFile) % Get rid of any previous savefile with same name


Parameters.Start.Year      = 2001; %start year
Parameters.Start.Month     =    9; %start month
Parameters.Start.Day       =    1; %start day

Parameters.End.Year        = 2002; %end year
Parameters.End.Month       =   12; %end month
Parameters.End.Day         =   31; %end day

Parameters.TimeWindow = 31;             %in days
Parameters.TimeStep = 31;                %in days

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % optional detection parameters
% Parameters.WindHeightGates = [78,83;(83:92)',(86:95)'; 95,100]; %km; height gates used for wind fitting - choose depending on radar
% Parameters.MFHeightGates   = [80,85;85,90;87.5,92.5;90,95;95,100]; %km; height gates for momentum flux binning - can be different to the above, e.g. to subtract tides
% Parameters.WindHeightGates = [78,83;83,86;86,89;89,92;92,95;95,100]; %km; height gates used for wind fitting - choose depending on radar
% Parameters.MFHeightGates   = [78,83;83,86;86,89;89,92;92,95;95,100]; %km; height gates for momentum flux binning - can be different to the above, e.g. to subtract tides
Parameters.WindHeightGates = [75: 0.25 :100; 80: 0.25 :105]';
Parameters.MFHeightGates = Parameters.WindHeightGates;

Parameters.ZenithLimits    = [15*pi/180 50*pi/180]; %radiansP
Parameters.MinMeteorsForMF = 30; % minimum number of meteors in box required to compute MF
Parameters.UnphysicalMF    = 300; %minimum cutoff to declare MF unphysical

% Removal parameters - AM: 29/05/2014

Parameters.ZenithBinWidth = 5;   % Bin Width/2 to use for removal - i.e. 10degrees  =  5 

Parameters.ZenithBins = (Parameters.ZenithLimits(1)*180/pi)+Parameters.ZenithBinWidth:...
                                                          Parameters.ZenithBinWidth*2:...
                        (Parameters.ZenithLimits(2)*180/pi)+Parameters.ZenithBinWidth;
                                 % Central points of each Zenith bin for removal 
                                
Parameters.NumSTDs = 3;          % Number of STD's to use for removing meteors radial velocity pertubations
                                 % if 0 no removal takes place. Number
                                 % doesn't need to be integer but in
                                 % general it is probably always going to
                                 % be.
%declare results arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hocking.ZonalVar = NaN(length(datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day):Parameters.TimeStep:datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day)),size(Parameters.MFHeightGates,1));
Hocking.MeridVar = NaN(length(datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day):Parameters.TimeStep:datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day)),size(Parameters.MFHeightGates,1));
Hocking.VertVar = NaN(length(datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day):Parameters.TimeStep:datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day)),size(Parameters.MFHeightGates,1));
Hocking.VertMF = NaN(length(datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day):Parameters.TimeStep:datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day)),size(Parameters.MFHeightGates,1));
Hocking.ZonalMF = NaN(length(datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day):Parameters.TimeStep:datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day)),size(Parameters.MFHeightGates,1));
Hocking.MeridMF = NaN(length(datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day):Parameters.TimeStep:datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day)),size(Parameters.MFHeightGates,1));
Hocking.Time = NaN(length(datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day):Parameters.TimeStep:datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day)),1); 
Hocking.Height = NaN(length(datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day):Parameters.TimeStep:datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day)),size(Parameters.MFHeightGates,1));
Hocking.Total = NaN(length(datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day):Parameters.TimeStep:datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day)),size(Parameters.MFHeightGates,1));

%fit winds to the meteors for this location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WindStructure.StartTime = datenum(Parameters.Start.Year, ...
                                  Parameters.Start.Month,...
                                  Parameters.Start.Day,0,0,0);
WindStructure.EndTime   = datenum(Parameters.End.Year, ...
                                  Parameters.End.Month,...
                                  Parameters.End.Day,23,59,59);                                
WindStructure.HeightGates  = Parameters.WindHeightGates;
WindStructure.DataDir = Parameters.MonthlyDir;
WindStructure.Silent = 2; % Don't display any time information for WindFitting.

WindStructure.BinSize = 2/24;

tic;[WindFit.Zonal, ...
 WindFit.Merid, ...
 WindFit.Times, ...
 WindFit.Alts] = WindCalc(Parameters.Site,WindStructure);toc;
% AM: (windcompute has been tested and agrees with Robin's code)

disp('---->Wind fitted for site and chosen wind gates.');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%core data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  StartTime = datenum(Parameters.Start.Year,Parameters.Start.Month,Parameters.Start.Day,12,00,00) + floor(Parameters.TimeStep/2);
  EndTime = datenum(Parameters.End.Year,Parameters.End.Month,Parameters.End.Day,12,00,00) - floor(Parameters.TimeStep/2);
  TimeArray = StartTime : Parameters.TimeStep : EndTime;
  
  % If using a Time Window of 31 use monthly timesteps. (N.B. This keeps
  % the timestep as set in parameters so months will overlap slightly. 
  if Parameters.TimeWindow == 31 && Parameters.TimeStep == 31; 
     TimeArray = datenum(Parameters.Start.Year,...
                         Parameters.Start.Month : (12-Parameters.Start.Month) ...
                                                   + ((Parameters.End.Year-Parameters.Start.Year)*12)...
                                                   + Parameters.End.Month-(12-Parameters.Start.Month),...
                         16);
  end
  for iTime = 1:length(TimeArray)
    tic;
      HockingResults = NaN(6,size(Parameters.MFHeightGates,1));
      
    TimeStep = TimeArray(iTime);
    disp(['Running ',datestr(TimeStep)])
    [iYear, iMonth, ~, ~, ~, ~] = datevecmx(TimeStep);

    clear MonthlyMeteorData %so it doesn't hang over from the previous month
  
    % AM 14/04/2014:
    %  This just adjusts for when you want last month of previous year  
    if iMonth == 1
        iYear = iYear-1;
        iMonth = 13;
    end

    % Best to pick middle height gate for this as NaN's for certain height
    % levels will remove whole line from data. Picking middle height gate
    % minimises this.  AM: 28/05/2014
    WindFit_forlocalwind.Zonal = WindFit.Zonal(nanmean(WindFit.Times,2) >= datenum(iYear,iMonth-1,01)...
                                             & nanmean(WindFit.Times,2) <= datenum(iYear,iMonth+1,31),:);
    WindFit_forlocalwind.Merid = WindFit.Merid(nanmean(WindFit.Times,2) >= datenum(iYear,iMonth-1,01)...
                                             & nanmean(WindFit.Times,2) <= datenum(iYear,iMonth+1,31),:);
    WindFit_forlocalwind.Times = WindFit.Times(nanmean(WindFit.Times,2) >= datenum(iYear,iMonth-1,01)...
                                             & nanmean(WindFit.Times,2) <= datenum(iYear,iMonth+1,31),:);
    WindFit_forlocalwind.Alts = WindFit.Alts(nanmean(WindFit.Times,2) >= datenum(iYear,iMonth-1,01)...
                                             & nanmean(WindFit.Times,2) <= datenum(iYear,iMonth+1,31),:);
    Parameterstemp = Parameters;
    Parameterstemp.Start.Year      = iYear; %start year
    Parameterstemp.Start.Month     = iMonth-1; %start month
    Parameterstemp.Start.Day       =    1; %start day

    Parameterstemp.End.Year        = iYear; %end year
    Parameterstemp.End.Month       = iMonth+1; %end month
    Parameterstemp.End.Day         =   31; %end day

                                                                                  
     PeriodMeteorData = MeteorWindLocal(Parameterstemp,TimeStep,WindFit_forlocalwind); 
     % AM: calculation of local wind for each meteor.
     
    if nansum(PeriodMeteorData) == 0; continue;end; %bad/no data
    
    %discard meteors we don't have wind data for    
    PeriodMeteorData = PeriodMeteorData(~isnan(PeriodMeteorData(:,13)),:);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %process individual meteor data to get MF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%AM: Pre-allocation for counting meteors and averaging Height 
TotalMeteors = NaN(size(Parameters.MFHeightGates,1),1);
MeanHeight = NaN(size(Parameters.MFHeightGates,1),1);

    for iHeight = 1:1:size(Parameters.MFHeightGates,1)       
      %select those meteors in this (height,time) box
      PeriodMeteorDataAtHeight = PeriodMeteorData(PeriodMeteorData(:,3) >= Parameters.MFHeightGates(iHeight,1) ...
        & PeriodMeteorData(:,3) <  Parameters.MFHeightGates(iHeight,2),:);
          
      %find the mean height of these meteors and the number of meteors used
      %in calculation
      TotalMeteors(iHeight,1) = size(PeriodMeteorDataAtHeight,1);
      MeanHeight(iHeight,1) = nanmean(PeriodMeteorDataAtHeight(:,3));
      
      %check enough Meteors
      if size(PeriodMeteorDataAtHeight,1)<Parameters.MinMeteorsForMF;continue;end;      
      
      %trim down data to just what we need for the rest of the analysis
      %[zenith, azimuth, horiz drift vel, radial velocity, zonal wind, merid wind] - rest discarded
      PeriodMeteorDataAtHeight = PeriodMeteorDataAtHeight(:,[4 5 6 12 13 14]);
      %extract background wind at the location of each meteor
      %assume no w wind and that all vertical perturbs are due to gws
      %need to be projected onto meteor
      BackgroundWind = [PeriodMeteorDataAtHeight(:,5),                   ... %u
                        PeriodMeteorDataAtHeight(:,6),                   ... %v
                        zeros(length(PeriodMeteorDataAtHeight(:,1)),1)];     %w
      
      %extract measured radial velocity, and split into xyz components      
      MeteorComponents =[PeriodMeteorDataAtHeight(:,4).*sin(PeriodMeteorDataAtHeight(:,1)).*cos(PeriodMeteorDataAtHeight(:,2)), ...
                         PeriodMeteorDataAtHeight(:,4).*sin(PeriodMeteorDataAtHeight(:,1)).*sin(PeriodMeteorDataAtHeight(:,2)), ...
                         PeriodMeteorDataAtHeight(:,4).*cos(PeriodMeteorDataAtHeight(:,1))];
      
                     
      %project background wind onto radial direction of meteor, and split into xyz components
      Proj = (MeteorComponents(:,1).*BackgroundWind(:,1) ...
        + MeteorComponents(:,2).*BackgroundWind(:,2) ...
        + MeteorComponents(:,3).*BackgroundWind(:,3))./PeriodMeteorDataAtHeight(:,4);
      ProjC =  [Proj Proj Proj] ...
        .* MeteorComponents ...
        ./ [PeriodMeteorDataAtHeight(:,4) PeriodMeteorDataAtHeight(:,4) PeriodMeteorDataAtHeight(:,4)];
      
      %hence, compute radial velocity PERTURBATION
      %this is: radial velocity in components minus background wind in components
      Resultc = MeteorComponents - ProjC;
      
      %also compute the magnitude of the radial perturbation velocity
      ModRadialVelocity = sqrt(Resultc(:,1).^2+Resultc(:,2).^2+Resultc(:,3).^2);
      
      %Check radial direction - is it going inwards or outwards?
      Direction = Resultc(:,3)./sqrt(Resultc(:,3).^2); % +Mve/-ve
      ModRadialVelocity = ModRadialVelocity.*Direction;
      
      %zenith of the meteor - angle from vertical (0 straight up)
      theta = PeriodMeteorDataAtHeight(:,1);
      %azimuth of the meteor - angle anticlockwise from due east
      phi = PeriodMeteorDataAtHeight(:,2);
      
      %tidying - some NaNs to be removed from not enough meteors in the windfitting:
      %necessary to avoid breaking matrix
      phi(isnan(ModRadialVelocity)) = [];
      theta(isnan(ModRadialVelocity)) = [];
      ModRadialVelocity(isnan(ModRadialVelocity)) = [];

      %AM: 28/05/2014 - std removal of extreme radial velocity 
      %                 pertubations
    if Parameters.NumSTDs ~=0  

      newModRadialVelocity = [];
      newTheta = [];
      newPhi = [];

      for bin = 1:length(Parameters.ZenithBins) 
        BinIndex = (theta*180/pi >= (Parameters.ZenithBins(bin) - Parameters.ZenithBinWidth) & ...
                    theta*180/pi < (Parameters.ZenithBins(bin)  + Parameters.ZenithBinWidth));
        tempModRadialVelocity = ModRadialVelocity(BinIndex);
        temptheta = theta(BinIndex);
        tempphi = phi(BinIndex);
        BinMean  = nanmean(tempModRadialVelocity);
        BinSTD   = nanstd(tempModRadialVelocity);
        RemIndex = (tempModRadialVelocity >= (BinMean + (Parameters.NumSTDs * BinSTD)) | ...
                          tempModRadialVelocity <= (BinMean - (Parameters.NumSTDs * BinSTD)));
        tempModRadialVelocity(RemIndex) = [];
        temptheta(RemIndex) = [];
        tempphi(RemIndex) = [];
        if isempty(tempModRadialVelocity) ~=1
        newModRadialVelocity = cat(1,newModRadialVelocity,tempModRadialVelocity);
        newTheta = cat(1,newTheta,temptheta);
        newPhi = cat(1,newPhi,tempphi);
        end
        clear tempModRadialVelocity temptheta tempphi BinMean BinSTD BinIndex RemIndex
      end

      ModRadialVelocity = newModRadialVelocity;
      theta = newTheta;
      phi = newPhi;
      clear newTheta newPhi newModRadialVelocity;

    end
      %compute the matrices from Hocking et al (2005), and store them
      %     [theta,    phi,    ModRadialVelocity]
      
      HockingResults(:,iHeight) = hockingmatrix(theta,phi,ModRadialVelocity);
      
      %tidy up
      clear BackgroundWind Direction MeteorComponents ModRadialVelocity ...
        PeriodMeteorDataAtHeight Proj ProjC RadialVelocity    ...
        Resultc phi theta
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %extract desired results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Hocking Results
    Hocking.ZonalVar(iTime,:) = HockingResults(1,:);        % Zonal Variance
    Hocking.MeridVar(iTime,:) = HockingResults(2,:);        % Merid Variance
    Hocking.VertVar(iTime,:) = HockingResults(3,:);         % Vertical Variance
    Hocking.VertMF(iTime,:) = HockingResults(4,:);          % Vertical Momentum Flux
    Hocking.ZonalMF(iTime,:) = HockingResults(5,:);         % Zonal Momentum Flux
    Hocking.MeridMF(iTime,:) = HockingResults(6,:);         % Meridional Momentum Flux
    Hocking.Time(iTime) = nanmean(PeriodMeteorData(:,1));   % Time of Calculations
    Hocking.Height(iTime,:) = MeanHeight;                   % Height of Calculations
    Hocking.Total(iTime,:) = TotalMeteors;                  % Number of Meteors used in Calculations. 

    clear iHeight PeriodMeteorData MeanHeight
   toc;
  end
    save(Parameters.OutFile,'WindFit','Parameters','Hocking');