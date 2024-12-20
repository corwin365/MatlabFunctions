function uiopen(type,direct)
%UIOPEN Present file selection dialog with appropriate file filters.
%
%   UIOPEN presents a file selection dialog.  The user can either choose a
%   file to open or click cancel.  No further action is taken if the user
%   clicks on cancel.  Otherwise the OPEN command is evaluated in the base
%   workspace with the user specified filename.
%
%   These are the file filters presented using UIOPEN.
%
%   1st input argument   Filter List
%   <no input args>      *.m, *.fig, *.mat,
%                        *.mdl, *.slx  (if Simulink is installed),
%                        *.cdr         (if Stateflow is installed),
%                        *.rtw, *.tmf, *.tlc, *.c, *.h, *.ads, *.adb
%                                      (if Simulink Coder is installed),
%                        *.*
%   MATLAB               *.m, *.fig, *.*
%   LOAD                 *.mat, *.*
%   FIGURE               *.fig, *.*
%   SIMULINK             *.mdl, *.slx, *.*
%   EDITOR               *.m, *.mdl, *.cdr, *.rtw, *.tmf, *.tlc, *.c, *.h, *.ads, *.adb, *.*
%
%   If the first input argument is unrecognized, it is treated as a file
%   filter and passed directly to the UIGETFILE command.
%
%   If the second input argument is true, the first input argument is
%   treated as a filename.
%
%   Examples:
%       uiopen % displays the dialog with the file filter set to all MATLAB
%              %files.
%
%       uiopen('matlab') %displays the dialog with the file
%                         %filter set to all MATLAB files.
%
%       uiopen('load') %displays the dialog with the file filter set to
%                      %MAT-files (*.mat).
%
%       uiopen('figure') %displays the dialog with the file filter set to
%                        %figure files (*.fig).
%
%       uiopen('simulink') %displays the dialog with the file filter set to
%                          %model files (*.mdl,*.slx).
%
%       uiopen('editor') %displays the dialog with the file filter set to
%                        %"All MATLAB files". This filters out binary
%                        %formats: .mat, .fig, .slx.
%                        %All files are opened in the MATLAB Editor.
%
%   See also UIGETFILE, UIPUTFILE, OPEN, UIIMPORT.

%   Copyright 1984-2023 The MathWorks, Inc.

% Check for Swing capability if web UI is disabled
import matlab.internal.capability.Capability;
if ~matlab.ui.internal.dialog.FileDialogHelper.isWebUI()
    Capability.require(Capability.Swing);
end

if nargin > 0
    type = convertStringsToChars(type);
end

if nargin > 1
    direct = convertStringsToChars(direct);
end

if nargin == 0
    type = '';
end

if nargin < 2
    direct = false;
end

if direct
    fn = type;
else
    if ~matlab.ui.internal.dialog.FileDialogHelper.isWebUI()
        % Error if MATLAB is running in no JVM mode.
        warnfiguredialog('uiopen');
    end
    if isempty(type)
        % Do not provide a filter list. UIGETFILE called below will pick up all
        % the filters available in the MATLAB installation automatically.
        filterList = '';
    else
        filterList = matlab.ui.internal.dialog.UIOpenHelper.getFilterList(type);
    end
    
    % Is it a .APP or .KEY directory on the Mac?
    % If so, open it properly.
    if ismac && ~iscell(filterList) && (ischar(filterList) && isfolder(filterList))
        [~, ~, ext] = fileparts(filterList);
        if strcmpi(ext, '.app') || strcmpi(ext, '.key')
            unix(['open "' filterList '" &']);
            return;
        end
    end
    [fn,pn] = uigetfile(filterList,getString(message('MATLAB:uistring:uiopen:DialogOpen')));
    if isequal(fn,0)
        return;
    end
    fn = fullfile(pn,fn);
end

try
    % send open requests from editor back to editor
    if strcmpi(type,'editor')
        edit(fn);
    else
        % Is it a MAT-file?
        [~, ~, ext] = fileparts(fn);
        if strcmpi(ext, '.mat')
            quotedFile = ['''' strrep(fn, '''', '''''') ''''];
            evalin('caller', ['load(' quotedFile ');']);
            matlab.ui.internal.dialog.UIOpenHelper.setStatusBar(~isempty(whos('-file', fn)));
            return;
        elseif strcmpi(ext, '.nc') | strcmpi(ext, '.nc4')
          %is it a netCDF file?
          quotedFile = ['''' strrep(fn, '''', '''''') ''''];
          evalin('caller', ['netCDF = rCDF(' quotedFile ');']);
          return;
        end
        
        sans = [];
        % If open creates variables, ans will get populated in here.
        % We need to assign it in the calling workspace later
        open(fn);
        
        if ~isempty(sans)
            vars = sans;
            % Shove the resulting variables into the calling workspace
            status = true;
            if isstruct(vars)
                names = fieldnames(vars);
                status = ~isempty(names);
                for i = 1:length(names)
                    assignin('caller',names{i}, vars.(names{i}));
                end
            else
                assignin('caller','ans',vars);
            end
            matlab.ui.internal.dialog.UIOpenHelper.setStatusBar(status);
        end
    end
catch ex
    err = ex.message;
    % Strip hyperlinks since errordlg does not support them
    err = regexprep(err, '<a\s+href\s*=\s*"[^"]*"[^>]*>(.*?)</a>','$1');
    errordlg(err);
end

end
