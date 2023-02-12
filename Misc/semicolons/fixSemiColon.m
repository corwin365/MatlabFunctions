function fixSemiColon(fileName)
% This function will accept the name of a MATLAB file and then replace add
% a semi-colon to all lines which do not have one.
%
% Note that currently this function cannot handle statements broken up
% into multiple lines using ellipses
%
messages = mlint(fileName,'-struct');

for ii = 1:length(messages)
    if strcmp(messages(ii).message,...
            'Terminate statement with semicolon to suppress output (within a script).')
        
        % Read and modify line that needs a semi-colon
        replaceLine = messages(ii).line;
        fileLine = readLineOfFile(fileName,replaceLine);
        fileLine = [fileLine ';'];
        
        replaceFileLine(fileName,replaceLine,fileLine)
        
    end
end
