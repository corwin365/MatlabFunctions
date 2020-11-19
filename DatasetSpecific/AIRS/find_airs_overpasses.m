function [Results,NotFound] = find_airs_overpasses(TimeRange,LonRange,LatRange,Truncate,Verbose)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find AIRS passes over a particular grid box in time and space
%return date and granule number of each relevant pass 
%any granules which could not be loaded are stored in NotFound
%
%uses Hoffmann + Alexander 3D data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off %this routine generates a lot of polynomial-duplication warnings that don't affect the results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time range
Settings.TimeScale = TimeRange(1):1:TimeRange(2); clear TimeRange

%lat/lon box
Settings.LonRange = LonRange; clear LonRange;
Settings.LatRange = LatRange; clear LatRange;

%truncate the airs granule?
%this cuts a specified number of pixels off each side of the 2D plane
%before computing the match polygon. this is useful, for example, to make 
%sure the granule is well-centred on the region rather than just clipping it
%order is [xt,at] and is symmetrical - to match just the central row of the
%granule in the AT direction but anywhere along the length, this could be
%set to [44,0];
%OPTIONAL VARIABLE
if ~exist('Truncate'); Settings.Truncate = [0,0];
else                   Settings.Truncate = Truncate;
end; clear Truncate;

%print how many we found? input anything except 0 to do so
%OPTIONAL VARIABLE
if exist('Verbose')
  if Verbose ~= 0; Verbose = 1; else Verbose = 0;
  end
else Verbose = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create polybox of the sample region
Poly.Box = polyshape(Settings.LonRange([1,2,2,1,1]),Settings.LatRange([1,1,2,2,1]));

%create results array
Results = []; %will expand out

%create list of files that we failed to locate
NotFound = []; %will expand out

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over period and find all the matching passes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.TimeScale)
  Count = 0;
  for iGranule=1:1:240;

    %load granule geolocation
    [Airs,~,Error,~] = prep_airs_3d(Settings.TimeScale(iDay),iGranule,'LoadOnly',true);
    
    %check if we loaded a granule successfully, and loop to next if not
    if Error ~=0; 
      NotFound(end+1,:) = [Settings.TimeScale(iDay),iGranule];
      clear Airs Error;
      continue; 
    end
    clear Error
    
    %truncate, if needed
    if nansum(Settings.Truncate) > 0;
      Airs.l1_lon = Airs.l1_lon(1+Settings.Truncate(1):end-Settings.Truncate(1), ...
                                1+Settings.Truncate(2):end-Settings.Truncate(2));
      Airs.l1_lat = Airs.l1_lat(1+Settings.Truncate(1):end-Settings.Truncate(1), ...
                                1+Settings.Truncate(2):end-Settings.Truncate(2));
    end
    
    %create a polybox for the AIRS granule area    
    Poly.Airs = polyshape([Airs.l1_lon(1:end,1)', ...
                           Airs.l1_lon(end,1:end), ...
                           Airs.l1_lon(end:-1:1,end)', ...
                           Airs.l1_lon(1,end:-1:1)], ...
                          [Airs.l1_lat(1:end,1)', ...
                           Airs.l1_lat(end,1:end), ...
                           Airs.l1_lat(end:-1:1,end)', ...
                           Airs.l1_lat(1,end:-1:1)]);

    %check if the box crosses the dateline, as this breaks things
    if range([min(Poly.Airs.Vertices(:,1)),max(Poly.Airs.Vertices(:,1))]) < 170;
      %normal behaviour - do the two regions overlap?
      TF = overlaps(Poly.Box,Poly.Airs);
    else
      %dateline crossing.
      
      %first, shift the data into a positive range
      Shift = find(Poly.Airs.Vertices(:,1) < 0);
      Poly.Airs.Vertices(Shift,1) = Poly.Airs.Vertices(Shift,1)+360;
      
      %then, check both the ground box and the duplicate of this in +360
      BoxA = Poly.Box;
      BoxB = Poly.Box; BoxB.Vertices(:,1) = BoxB.Vertices(:,1)+360;
      TF1 = overlaps(BoxA,Poly.Airs);
      TF2 = overlaps(BoxB,Poly.Airs);
      TF = max([TF1,TF2]); %if either are true this will be 1
      
    end
     

    
    %if we have a match, store it
    if TF == 1; 
      Results(end+1,:) = [Settings.TimeScale(iDay),iGranule]; 
      Count = Count+1;
    end
    
    clear Airs TG; Poly = rmfield(Poly,'Airs'); 
      
  end; clear iGranule
  if Verbose == 1; disp(['Processed ',datestr(Settings.TimeScale(iDay)),' - ',num2str(Count),' matches today']); end
  clear Count
end; clear iDay


warning on %switch it back on
    