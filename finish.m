
%store last working directory, for interactive sessions only: -batch runs (which may
%run several at once) neither need it nor should overwrite it, and an error here
%must never stop MATLAB from exiting
%  2026/11/25 [Claude] skipped in -batch runs; setpref errors caught (parallel -batch runs hung at exit)
disp(newline)
disp(newline)
disp('--------------------------------------------------------------------------------')
if ~batchStartupOptionUsed
  try
    setpref('StartupDirectory','LastWorkingDirectory',pwd)
    disp(['Last working directory retained: ',pwd])
  catch Err
    disp(['Could not store the last working directory: ',Err.message])
  end
end


%done!
disp('Goodbye!');
disp('--------------------------------------------------------------------------------')
