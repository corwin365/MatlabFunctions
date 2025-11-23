function FormattedString = yyyyDDD(Date)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert a datenum() or datetime() to a char of the format YYYYdDDD,
%where YYYY is the year and DDD is the day number within the year
%
%modified 28/01/2025 to add datetimes
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/11/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(class(Date),'datetime')
  FormattedString = [sprintf('%04d',year(Date)),'d',sprintf('%03d',day(Date(1),'dayofyear'))];


else

  [y,~,~] = datevec(Date);
  dn = floor(date2doy(Date));
  FormattedString = [sprintf('%04d',y),'d',sprintf('%03d',dn)];

end

return
end
