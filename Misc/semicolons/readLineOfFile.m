function fileLine = readLineOfFile(fileName,lineNumber)
% This function will read and return a specific line of a text file
%
% For example:
% fileLine = readLineOfFile('testFile.m',5);
%
fid = fopen(fileName,'r');
for k=1:(lineNumber-1);
    fgetl(fid);
end;

fileLine = fgetl(fid);
fclose(fid);