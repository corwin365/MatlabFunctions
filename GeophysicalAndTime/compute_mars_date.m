function [MarsYear,SolarLongitude,CoordinatedMarsTime] = compute_mars_date(MatlabDate)
% % MatlabDate = datenum(now)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert Matlab date to Mars date
%
%Corwin Wright, c.wright@bath.ac.uk, 17/AUG/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%following the method of
%%http://jtauber.github.io/mars-clock/
%and all arbitrary-looking numbers are taken from there

JDut = juliandate(MatlabDate);

JDtt = JDut + (37+32.184)/86400;

DaysSinceJ2000Epoch = JDtt - 2451545.0;

MarsSolDate = ((DaysSinceJ2000Epoch  -4.5) ./ 1.027491252)+44796 - 0.00096;

CoordinatedMarsTime = mod(24*MarsSolDate,24);

MarsMeanAnomaly = mod(19.3870 + 0.5240275 .* DaysSinceJ2000Epoch,360);

aFMS = mod(270.3863 + 0.52403840 .* DaysSinceJ2000Epoch,360);

e = 0.09342;

vminusM = (10.691 + 3e-7 .* DaysSinceJ2000Epoch).*sind(MarsMeanAnomaly) ...
        + 0.623  .* sind(2 .* MarsMeanAnomaly) ...
        + 0.050  .* sind(3 .* MarsMeanAnomaly) ...
        + 0.005  .* sind(4 .* MarsMeanAnomaly) ...
        + 0.0005 .* sind(5 .* MarsMeanAnomaly);
      
SolarLongitude = aFMS + vminusM; 



%finally, year. lookup table of start date for years 1-40. ONLY WORKS FROM 1955 TO 2028
YearStarts = datenum(['Apr 11 1955';'Feb 26 1957';'Jan 14 1959';'Dec 01 1960';'Oct 19 1962';'Sep 05 1964';'Jul 24 1966';'Jun 10 1968';'Apr 28 1970';'Mar 15 1972';'Jan 31 1974';'Dec 19 1975';'Nov 05 1977';'Sep 23 1979';'Aug 10 1981';'Jun 28 1983';'May 15 1985';'Apr 01 1987';'Feb 16 1989';'Jan 04 1991';'Nov 21 1992';'Oct 09 1994';'Aug 26 1996';'Jul 14 1998';'May 31 2000';'Apr 18 2002';'Mar 05 2004';'Jan 21 2006';'Dec 09 2007';'Oct 26 2009';'Sep 13 2011';'Jul 31 2013';'Jun 18 2015';'May 05 2017';'Mar 23 2019';'Feb 07 2021';'Dec 26 2022';'Nov 12 2024';'Sep 30 2026';'Aug 17 2028']);
[~,MarsYear] = min(abs(MatlabDate-YearStarts));


FudgeFactor = 20; %no idea why this offset is needed, but it fixes a bug. unsatisfying. 

for iPoint=1:1:numel(MarsYear)
 if MatlabDate(iPoint) <= YearStarts(MarsYear(iPoint))+FudgeFactor
   MarsYear(iPoint) = MarsYear(iPoint)-1; 
 end
end



