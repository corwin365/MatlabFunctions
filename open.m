function out = open(name)
%OPEN	 Open file in appropriate application.
%   OPEN NAME opens the file or variable specified by NAME in the
%   appropriate application. For example, variables open in the Variables
%   Editor, MATLAB code files open in the MATLAB Editor, and Simulink
%   models open in Simulink.
%
%   If NAME is a MAT-file, MATLAB returns the variables in NAME to a
%   structure.
%
%   If NAME does not include an extension, MATLAB searches for variables
%   and files according to the function precedence order.
%
%   A = OPEN(NAME) returns a structure if NAME is a MAT-file or a figure
%   handle if NAME if a figure. Otherwise, OPEN returns an empty array.
%
%   Examples:
%
%     open f2                   First looks for a variable named f2, then
%                               looks on the path for a file named f2.slx,
%                               f2.mdl, f2.mlapp, f2.mlx, or f2.m.  Error
%                               if can't find any of these.
%
%     open f2.mat               Error if f2.mat is not on path.
%
%     open d:\temp\data.mat     Error if data.mat is not in d:\temp.
%
%
%   OPEN is user-extensible.  To open a file with the extension ".XXX",
%   OPEN calls the helper function OPENXXX, that is, a function
%   named 'OPEN', with the file extension appended.
%
%   For example,
%      open('foo.log')       calls openlog('foo.log')
%      open foo.log          calls openlog('foo.log')
%
%   You can create your own OPENXXX functions to set up handlers
%   for new file types.  OPEN will call whatever OPENXXX function
%   it finds on the path.
%
%   See also EDIT, LOAD, OPENFIG, OPENVAR, WHICH, WINOPEN.
%

%   Copyright 1984-2023 The MathWorks, Inc.
arguments
    name {mustBeTextScalar, mustBeNonzeroLengthText};
end

assert(~isdeployed, message('MATLAB:open:noOpen'));

name = convertStringsToChars(name);

if nargout
    out = [];
end

% If we found a variable that matches, use that.
% Open the variable, and get out.

if isvarname(name) && evalin('caller', "exist('" + name + "', 'var')")
    evalin('caller', "openvar('" + name + "', " + name + ");");
    return;
else
    % special handling for netCDF files
    [~,~,x] = fileparts(name);
    if strcmpi(x,'.nc')
        data = rCDF(name);
        % always write to variable 'netCDF' in caller workspace
        assignin('caller','netCDF',data);
        if nargout
            out = data;
        end
        return
    end
end

fullpath = matlab.lang.internal.introspective.safeWhich(name);
if isempty(fullpath)
    if isfile(name)
        fullpath = name;
    else
        nameResolver = matlab.lang.internal.introspective.resolveName(name);
        classInfo = nameResolver.classInfo;
        if isempty(classInfo) || ~(classInfo.isClass || classInfo.isMethod)
            fullpath = nameResolver.nameLocation;
        end
    end
end

if isempty(fullpath)
    if isfolder(name)
        error(message('MATLAB:open:noOpenFolder', name));        
    else
        if ~isempty(name) && edit(name)
            return;
        end

        error(message('MATLAB:open:fileNotFound', name));
    end
end

% If no extension was specified, call which with a '.' appended, so that
% we can see if the exact match is available.
if ~hasExtension(name)
    %Get all files/dirs which have just the name
    tmpPath = which([name '.'], '-all');
    if ~isempty(tmpPath)
        for i = 1:length(tmpPath)
            %If we find a file, set the path to it, and stop.  This means
            %we find files in the same order as which -all returns them.
            if isfile(tmpPath{i})
                fullpath = tmpPath{i};
                break;
            end
        end
    end
end

if ~isfile(fullpath)
    error(message('MATLAB:open:fileNotFound', fullpath));
else
    % let finfo decide the filetype
    [~, openAction] = finfo(fullpath);
    if isempty(openAction)
        openAction = @defaultopen;
        % edit.m does not opens p files
        % check here if the .p extension was supplied by the user
        % If the user did not specify .p then which command appended the .p and
        % we need to strip it off before calling openp.
    elseif strcmp(openAction, 'openp')
        [~,~,ext] = fileparts(name);
        % g560308/g479211 is there wasn't an extension specified and a .p file
        % was found, then search for an associated .m file.
        if ~strcmp(ext, '.p')
            fullpath = fullpath(1:end-2);
            % if the .m file associated with the .p file does not exist, error out.
            if ~isfile([fullpath, '.m'])
                error(message('MATLAB:open:openFailure', [fullpath,'.p']));
            end
            
        end
        
    end
    
    try
        % if opening a mat file be sure to fetch output args
        if any(strcmp(openAction, ["openmat","openparquet"])) || nargout
            out = feval(openAction,fullpath);
        else
            feval(openAction,fullpath);
        end
    catch exception
        % we only want the message from the exception, not the stack trace
        throw(exception);
    end
end

%------------------------------------------
% Helper method that determines if filename specified has an extension.
% Returns 1 if filename does have an extension, 0 otherwise
function result = hasExtension(s)
[~,~,ext] = fileparts(s);
result = ~isempty(ext);

%------------------------------------------
function out = defaultopen(name)
% Default action to open unrecognized file types.
out = [];
edit(name);
