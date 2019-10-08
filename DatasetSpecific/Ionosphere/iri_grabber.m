function raw_iri_profiles = iri_grabber(start_time, end_time, station_lat, station_lon, thisyear, thismonth, thisday)
% Take IRI profiles from url, uses iri_2012

% url = 'http://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi';
% url = 'http://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi';
    

  url    = 'https://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi';
  params = {'model','iri_2016',...
                    'year', num2str(thisyear),...
                    'month',num2str(thismonth),...
                    'day', num2str(thisday), ...
                    'time_flag','0', ...
                    'hour','1.5', ...
                    'geo_flag', '0', ...
                    'latitude',num2str(station_lat), ...
                    'longitude',num2str(station_lon), ...
                   'height','100', ...
                   'profile','8', ...
                   'start','0', ...
                   'stop','24', ...
                   'step','0.25', ...
                  ... %'sun_n','', ...
                   ...%'ion_n','', ...
                   ...%'radio_f','', ... 
                   ...%'radio_f81','', ...
                   ...%'htec_max','', ... 
                   ...%'ne_top','', ...
                   ...%'imap','', ...
                  ... %'ffof2','fof2', ...
                   'ib0','0', ...
                   'probab','0', ...
                   ...%'ffoe',num2str(foe), ...
                   'dreg','0', ...
                   'tset','0', ...
                   'icomp','0', ...
                   'nmf2','0', ...
                   'hmf2','0', ...
                   ...%'user_nme','0', ...
                   ...%'user_hme','0', ...
                   'format','0', ...
                   'vars','00', ...
                   'vars','01', ...
                   'vars','02', ...
                   'vars','04',...
                   'vars','31', ... % hmf2 
                   'vars','33', ... % hmE 
                   'vars','44', ... % fof2
                   'vars','46', ... % foe                 
                    };
  [paramString,header] = http_paramsToString(params,1);
  [data,extras] = urlread2(url,'POST',paramString,header);


                
%%                
%find the start of the sequence where the real data start
Sequence = '<pre>';
Start = strfind(data,Sequence)+4; %gets rid of the <pre>
%and the end
Sequence = '</pre>';
End = strfind(data,Sequence)-2; %the -2gets rid of the leading < and a line break
data = data(Start:End);
clear Start End Sequence

%ASSUMING NOTHING CHANGES, this then removes all the metadata before the array
data  =data(145:end);

%finally, convert to an array of doubles
data = str2num(data);                
                
                
%%                
raw_iri_profiles = data;
end