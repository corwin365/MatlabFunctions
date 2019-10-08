function [MatlabTime] = cjw_time_airs2matlab(AirsTime,Inverse)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert airs time to matlab time
%
%airs format is: "floating-point elapsed seconds since Jan 1, 1993"
%approx - does not account for leap seconds, etc, but good enough for our purposes!
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%22/JAN/2014
%
%inputs
%---------
%
%AirsTime - time in AIRS units
%Inverse - matlab -> airs
%
%outputs
%---------
%
%MatlabTime - time in Matlab units
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; Inverse = 0; elseif Inverse ~= 0; Inverse = 1; end;

SecondsPerDay = 24*60*60.;
AirsEpoch     = datenum('1993','yyyy');

if Inverse == 0;
  MatlabTime    =  AirsEpoch+ AirsTime/(SecondsPerDay);
elseif Inverse == 1;
  %names hre are the wrong way around, for output handling
  MatlabTime = ((AirsTime-AirsEpoch).*SecondsPerDay); 
end