function [Frames,Error,Meta] = read_airglow_frames(Date,Band,Site,DataDir,Download)
%%%
%%%locates and reads in airglow frame gifs from Boston Univ format
%%%can download if needed, set Download=1 to do so.
%%%while downloading, times are also constructed from the webpage info
%%%
%%%errors:
%%% 1. file download failed
%%% 2. band unknown, no altitude value set
%%% 3. problem reading avi file
%%% 4. problem getting timestamps
%%%
%%%list of sites:
%%%  arecino,asiago,colombia,el_leoncito,germany,jicamarca,
%%%  mcdonald, mercedes, millstone, mt_john, rio_grande,
%%%  rothera, sutherland
%%%
%%%Corwin Wright, c.wright@bath.ac.uk
%%%2017/FEB/02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set up processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%defaults
if nargin < 3; Meta.Site     = 'millstone'; end
if nargin < 4; DataDir  = [LocalDataDir,'/AirglowFrames/']; end
if nargin < 5; Download = 0; end


%assume success unless proven otherwise
Error = 0;

%meta flagging
Meta.Site = Site;
Meta.Band = Band;
clear Site Band

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify file path and metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%convert date to filename string
%Aug0316
DateString = datestr(Date,'mmmddyy');

%approx altitude of band (km)
switch Meta.Band
  case 5577; Meta.Altitude = 96;
  case 6050; Meta.Altitude =  0;
  case 6444; Meta.Altitude =  0;
  case 5727; Meta.Altitude =  0;    
  case 6300; Meta.Altitude =  300;
  case 6950; Meta.Altitude =  87;
  case 5893; Meta.Altitude =  90;
  otherwise; Meta.Altitude = []; Error = 2;  %altitude unknown
end

%path
FileName = [DateString,'_',num2str(Meta.Band)];
FilePath = [DataDir,Meta.Site,'/',FileName,'.data.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if needed, download the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Download ==1;
  if exist(FilePath) ~= 2;%don't repeat!
    disp(['Don''t have frames stored for ',Meta.Site,' on ',datestr(Date),' at ',num2str(Meta.Band),'a; trying to download']);
    
   %get the file
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %get the time stamps an links for the frames, by parsing the source page
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [y,m,d] = datevec(Date);
    WebPath = ['http://sirius.bu.edu/data/?location=',Meta.Site, ...
                                          '&year=',num2str(y)...
                                          '&filt=',num2str(Meta.Band) ...
                                          '&month=',datestr(Date,'mmm'),...
                                          '&day=',num2str(d)];
    %download
    try
      Page = webread(WebPath);
    catch
      Frames = []; Error = 4;
      return;
    end
    clear WebPath
    
    
 
    %times are on lines of the form '\n</a><br/><b>hh:mm:ss UT</b>'
    %so find all chunks of text between each time these patterns occur
    %first split the file into lines
    Page = strsplit(Page,'\n');
    %then loop over the lines
    AnimLine = 0;
    Times = [];
    Frames = {};
    for iLine = 1:1:numel(Page);
      Line = Page{iLine};
      if numel(Line) < 12; continue; end
      if strcmp(Line(1:10),'<img src="') == 1;
        %animation frame: grab it
        
        %first, find path
        [startIndex,endIndex] = regexp(Line,'src\="[a-z0-9A-Z_./]+"');
        Path = Line(startIndex+5:endIndex-1); clear startIndex endIndex
        Path = ['http://sirius.bu.edu/data/',Path];
        
        %next, download
        Image = webread(Path);
        
        %and store
        if ~exist('AirGlow'); AirGlow = Image;
        else                   AirGlow = cat(3,AirGlow,Image);
        end
        
        %the time is in the filename, so pull it out
        [startIndex,endIndex] = regexp(Path,'/[A-Z0-9a-z]+_.gif');
        Time = Path(startIndex+2:endIndex-9);
        Time = datenum(Time,['HHMMSS']);
        Time = Time-floor(Time); Time = Time+Date;

        %and store
        Times(end+1) = Time;
      end
    end
    clear page AnimLine WebPath Time y d Page
    
    %save the results
    if exist('AirGlow');
      save(FilePath,'Times','AirGlow');
      clear Times AirGlow
    end

    clear y m d
  end
  
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read in times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except FilePath Error Meta

try
  load(FilePath);
  Meta.Time = Times;
  Frames = AirGlow;
catch
  Error = 3; Frames = []; Meta = [];
end


return

end

