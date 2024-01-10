function DayNumberList = samedays_everyyear(Years,DoYs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate a list of Matlab day numbers for the same set
%of days every specified year
%
%Corwin Wright, c.wright@bath.ac.uk
%2024/01/10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Years,DoYs] = meshgrid(Years,DoYs);

DayNumberList = datenum(Years(:),1,DoYs(:));


end