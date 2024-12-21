function FormattedString = yyyyDDD(DateNum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert a datenum() to a char of the format YYYYdDDD,
%where YYYY is the year and DDD is the day number within the year
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/11/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,~,~] = datevec(DateNum);
dn = floor(date2doy(DateNum));
FormattedString = [sprintf('%04d',y),'d',sprintf('%03d',dn)];

return
end