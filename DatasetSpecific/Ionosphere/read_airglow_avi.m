function [Frames,Error,Meta] = read_airglow_avi(Date,Band,Site,DataDir,Download)
%%%
%%%locates and reads in airglow avis from Boston Univ format
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
if nargin < 4; DataDir  = [LocalDataDir,'/Airglow/']; end
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
  %from Smith et al (JGR, 2000)
  case 5893; Meta.Altitude = 90; %sodium
  case 6444; Meta.Altitude =  0; %background
  case 6950; Meta.Altitude = 87; %oxygen
  case 7774; Meta.Altitude =  1; %cloud
  %from Madrigal FPI metadata
  case 5577; Meta.Altitude = 120; %green line
  case 6300; Meta.Altitude = 230; %red line     
% % %   %more to do!!
% % %   case 6050
% % %   case 6563
  otherwise; Meta.Altitude = []; Error = 2;  %altitude unknown
end

%path
FileName = [DateString,'_',num2str(Meta.Band),'.avi'];
FilePath = [DataDir,Meta.Site,'/',FileName];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if needed, download the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Download ==1;
  if exist(FilePath) ~= 2;%don't repeat!
    
   %get the file
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %identify file path on server
    [y,m,d] = datevec(Date);
    WebPath = ['http://sirius.bu.edu/data/data/',Meta.Site,'/', ...
                num2str(y),'/',num2str(Meta.Band),'/', ...
                datestr(Date,'mmm'),'/movies/',FileName];
    
    
    %download
    try
      websave(FilePath,WebPath);
    catch
      Frames = []; Error = 1;
      return;
    end
    
    %get the time stamps for the frames, by parsing the source page
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    for iLine = 1:1:numel(Page);
      Line = Page{iLine};
      if numel(Line) < 12; continue; end
      if strcmp(Line(1:12),'</a><br/><b>') == 1;
        %the first such line is a link the the animation, so skip
        if AnimLine == 0; AnimLine = 1; continue; end
        %we have a time. reformat and store.
        Time = [sprintf('%04d',y),'-',sprintf('%02d',m),'-',sprintf('%02d',d),' ',Line(13:20)];
        Times(end+1) = datenum(Time,['yyyy-mm-dd HH:MM:SS']);
      end 
    end
    clear page AnimLine WebPath Time y d Page
    
    %save the times
    TimeFile = [FilePath,'.times.mat'];
    save(TimeFile,'Times');
    clear Times

    clear y m d
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read in file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create video objexct, and extract frames
%there are only a hundred or so frames per night at most, so just read all
try
  v = VideoReader(FilePath);
  while hasFrame(v)
    
    if ~exist('Frames'); Frames = flipud(rgb2gray(readFrame(v)));
    else Frames = cat(3,Frames,flipud(rgb2gray(readFrame(v))));
    end
    
  end
  clear v
catch; Error = 3; Frames = []; return; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read in times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([FilePath,'.times.mat']);
Meta.Time = Times;


return

end

