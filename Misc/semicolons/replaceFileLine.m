function replaceFileLine(fileName,replaceLine,fileLine)
% This function will replace a specific line of a text file with another
% line
%
% Use fgetl to place the file position indicator at the line
fid = fopen(fileName,'r+');
for k=1:(replaceLine-1);
    fgetl(fid);
end;

% Call fseek between read and write operations
fseek(fid, 0, 'cof');

% Print the new values
fprintf(fid, fileLine)

% Close the file
fclose(fid);
