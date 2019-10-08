%  This program reads in a .mpd file from a Skiymet meteor radar and outputs the
%  parameters stored in the file as a matrix (mpd_data).  
%
%  INPUT
%  inputfile = text string containing full path filename of input .mpd file
%
%  OUTPUT
%  The format of the output matrix (mpd_data) is as follows:
%           column 1:  Digital time
%           column 2:  hour
%           column 3:  Height
%           column 4:  Zenith angle (theta)
%           column 5:  Azimuth (phi)
%           column 6:  Horizontal drift velocity (V Horiz)
%           column 7:  Peak amplitude of meteor signal (amax)
%           column 8:  Multiple meteor entry (0 or 1)
%           column 9:  Range
%           column 10: Ambig
%           column 11: Tau
%           column 12: snrdb (Signal to noise ratio in db)
%
%  Written 19th January, 2006, by David Sandford (University of Bath).
%  1st June 2006 - Signal to noise ratio added if mpd version is 2.1 or 2.2
%
%

function [mpd_data] = readmpd(inputfile)

%  Program variables
n = 0;          %  Simple counter (number of file lines read)
mpd_data = [];   %  Output matrix
thetamax = 75*pi/180;

%  Open .mpd file
fid2 = fopen(inputfile, 'rt');

data = [];

%fprintf('     Reading .mpd file')
while ~feof(fid2)
    n=n+1;
    tline = fgetl(fid2);
    if tline ~= -1
        if n == 1
            version = sscanf(tline,'%*s %g');
            nhead = 28;
        end
        if n == 2
            station = sscanf(tline,'%*s %s');
            if version == 2.2
                nhead = 28; %30;     %  Number of header lines
                Vers = 1;
            elseif version == 2.1 && strcmp(station, 'ascension') == 1
                nhead = 26; %30;     %  Number of header lines
                Vers = 1;
            elseif version == 2
                nhead = 20;
                Vers = 2;
            elseif version == 2.1 && strcmp(station, 'esrange') == 1
                nhead = 26;
                Vers = 2;
            else
                disp('Oh Crap!!!')
            end
        end
    

            
        if n > nhead;                           %<----- Ignore file header
            if Vers == 1
                data = [data;fscanf(fid2,'%g %*c %g %*c %g %g %*c %g %*c %g %*s %g %g %g %g %g %g %g %g %g %g %g %g %g %g',[20 inf])'];
            elseif Vers == 2
                data = fscanf(fid2,'%g %*c %g %*c %g %g %*c %g %*c %g %*s %g %g %g %g %g %g %g %g %g %g %g %g %g',[19 inf])';
            end
        end
    end
end

if isempty(data) == 0
    if Vers == 1
        mpd_data = [data(:,4)+data(:,5)/60+data(:,6)/3600 data(:,4) data(:,8) (pi*data(:,11)/180) (pi*data(:,12)/180) data(:,9) data(:,17) ones(size(data,1),1) data(:,7) data(:,13) data(:,18) data(:,20)];
    elseif Vers == 2
        mpd_data = [data(:,4)+data(:,5)/60+data(:,6)/3600 data(:,4) data(:,8) (pi*data(:,11)/180) (pi*data(:,12)/180) data(:,9) data(:,17) ones(size(data,1),1) data(:,7) data(:,13) data(:,18)];
    end
    
    clear data

    for m = 1:(size(mpd_data,1)-1)
        if mpd_data(m,10) == 1 && mpd_data(m,4) < thetamax && mpd_data(m,4) >10*pi/180
            if abs(mpd_data(m+1,1)-mpd_data(m,1)) < 1/3600
                if mpd_data(m+1,7) < mpd_data(m,7)
                    mpd_data(m+1,8) = 0; % ===== Change back to 0 =====
                else
                    mpd_data(m,8) = 0; % ===== Change back to 0 =====
                end
            end
        end
    end
    mpd_data = mpd_data(find(mpd_data(:,8) == 1 & mpd_data(:,10) == 1),:);
end

%  Close input file
fclose(fid2);

%==========================================================================

%  Criteria:
%  Tau is > 0.015 noise spike
if isempty(mpd_data) == 0
    mpd_data = mpd_data(find(mpd_data(:,11) >= 0.015 & mpd_data(:,11) < 1.0),:);
    
%  Correction for the curvature of the Earth
    Re = 6378.1;
    mpd_data(:,3) = (sqrt(Re^2 + mpd_data(:,9).^2 - (2 * Re.* mpd_data(:,9).* cos(pi - mpd_data(:,4)))) - Re);
    mpd_data(:,4) = acos((mpd_data(:,9).^2 + (Re + mpd_data(:,3)).^2 - Re^2)./(2.*mpd_data(:,9).*(Re + mpd_data(:,3)))) ;
    % or mpd_data(:,4) = asin((sin(pi - mpd_data(:,4)).*Re)./(Re + mpd_data(:,3)));
    % Both give same result
    
%  Convert radial velocity to horizontal velocity  ( Horizontal velocity = Radial velocity * sin(Zenith Angle) )
    mpd_data(:,6) = mpd_data(:,6)./ sin(mpd_data(:,4));
    
% Horiz Vel is less than 200 ms-1
    mpd_data = mpd_data(find(abs(mpd_data(:,6)) < 200),:);
    
%  Knock out meteors which have a zenith of less than 10 degrees
    mpd_data = mpd_data(find(mpd_data(:,4) > ((10/180)*pi)),:);
    
%  Knock out meteors which have a zenith of greater than 75 degrees
    mpd_data = mpd_data(find(mpd_data(:,4) < ((75/180)*pi)),:);
    
end

%==========================================================================

% if size(mpd_data,1) > 0
%     %fprintf(' \n')
%     %fprintf('     mpd file succussfully read \n')
%     fprintf('     Number of meteors detected: %g \n', size(mpd_data,1))
%     fprintf(' \n')
% else
%     fprintf('     .mpd file empty \n')
%     fprintf(' \n')
% end