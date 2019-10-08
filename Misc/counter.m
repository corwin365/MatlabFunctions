function OutputString = counter(iLoop,NLoops,DescString,NPerLine)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display an inline counter for matlab loops, adding a dot for each
%loop 1-N and a new line for every Nth loop
%Corwin Wright, corwin.wright@trinity.oxon.org
%05/JAN/2014
%
%inputs
%---------
%
%iLoop       - loop number
%NLoops      - total number of process
%DescString  - (optional) text string desciring what's happening
%NPerLine    - (optional) number of dots before a new line starts
%
%outputs
%---------
%
%outputs text to string - function output is just a dummy variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%input handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4; NPerLine = 1;                    end;
if nargin<3; DescString='Processing  loop';   end;

%output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(iLoop,NPerLine) == 1;
  
  %check we don't say something like 'loops 1 - 5 of 1'.
  MaximumLoops = iLoop+NPerLine-1;
  if NLoops < MaximumLoops; MaximumLoops = NLoops;end; 
  
  %output current progress to screen
  disp([DescString,' ',num2str(iLoop),' - ',num2str(MaximumLoops),' of ',num2str(NLoops)]);
  fprintf('\b');
  
end; 

if     mod(iLoop,20) == 0; fprintf('O'); ...
elseif mod(iLoop,10) == 0; fprintf('o'); ...
elseif mod(iLoop, 5) == 0; fprintf(','); ...
else;                      fprintf('.'); ...
end;

if mod(iLoop,NPerLine) == 0;
  fprintf('\n');
end;

%dummy variable for 'output'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OutputString = ''; 

  