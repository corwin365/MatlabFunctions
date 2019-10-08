% Program to load in a month of mpd files, convert to matlab and catenate

clear

site = 'riogrande';%rothera-sk, riogrande
daycount = zeros(10,12);
data = [];
for year = 2008:2008
    monda = [31,28,31,30,31,30,31,31,30,31,30,31];
    if year/4 == round(year/4)
        monda(2) = 29;
    end
    for month = 6:12
        data = [];
        for day = 1:monda(month);
            t = datenum(year,month,day);
            if exist(['/Users/robindavis/Desktop/GravWavMomFlux/data/',site,'MPD/mp',num2str(datestr(t,'yyyy')),num2str(datestr(t,'mm')),num2str(datestr(t,'dd')),'.',site,'.mpd']) ~=0
            dat = readmpd(['/Users/robindavis/Desktop/GravWavMomFlux/data/',site,'MPD/mp',num2str(datestr(t,'yyyy')),num2str(datestr(t,'mm')),num2str(datestr(t,'dd')),'.',site,'.mpd']);
            if length(dat)>=1
            dat(:,1) = dat(:,1)/24+datenum(year,month,day);
            data = ([data;dat]);
            end
            clear dat
            end
        end
        fprintf(',')
        save(['/Users/robindavis/Desktop/GravWavMomFlux/data/MonthMeteors/',site,datestr(datenum(year,month,15),'mmm'),num2str(datestr(datenum(year,month,15),'yyyy'))],'data')
    end
end       