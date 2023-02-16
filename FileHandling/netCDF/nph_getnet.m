
% Neil Hindley, University of Bath, UK. 2017.

% Read a netcdf file's variables, attributes and other into a structure

% Now called nph_getnet.m to distinguish from previous versions, which now
% just call this function.

% EDIT: 2023-02-13: Fixed a bug where all variables, where strings or
% numerics, were getting converted to numerics by single() and double()
% commands.

% EDIT: 2022-09-09 getnet.m has been working well for a while now, but I've
% encoutered an issue is the data is split into Groups. If the variables
% are in Groups, getnet doesn't find them. Fortunately, the ncinfo function
% finds and extracts Groups - it just doesn't extract the Data. I've fixed
% this now into structured data groups.

% EDIT 2018-05-24: Okay, getnet's been long overdue an overhaul, here
% I've added the functionality to read attributes from individual
% variables, but have had to reorganise the structuring so I'm afraid it's
% not retro-compatible any more. Sorry past-neil :(

% EDIT Also added the ability to scale and add offsets, if the variable's
% attributes are specified EXACTLY right.

% EDIT, okay so I think the output is going to look quite like the output
% of ncinfo, only with the variables' values actually in as well.

% Optional Arguments:
%
% 'double' | 'single' : choose whether to output the
% data fields as double or single (useful if data is massive).
%


function NC = nph_getnet(filepath,varargin)

%% Open netcdf file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nci = ncinfo(filepath);
ncid = netcdf.open(filepath,'nc_nowrite');

% Assign output:
NC = nci;

%% Attributes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NC.Attributes = struct;
NC.Attributes.Global = nci.Attributes;

% ALSO ASSIGN ATTRIBUTES FROM VARIABLES TO THE ATTRIBUTES FIELD:
for i = 1:length(nci.Variables)
    NC.Attributes.(nci.Variables(i).Name) = nci.Variables(i).Attributes;
end

%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% If all the data is in variables, just load the variables...

if any(strcmpi(varargin,'single'))
    singleflag = 1;
else
    singleflag = 0;
end

NC.Data = struct;

for i = 1:length(nci.Variables)
    
    numericflag = 1;
    if any(strcmpi(nci.Variables(i).Datatype,{'char','string'}))
        numericflag = 0;
    end
    
    try
        % load data as netcdf
        if numericflag
            if singleflag
                NC.Data.(nci.Variables(i).Name) = single(netcdf.getVar(ncid,i-1));
            else
                NC.Data.(nci.Variables(i).Name) = double(netcdf.getVar(ncid,i-1));
            end
        else
            NC.Data.(nci.Variables(i).Name) = netcdf.getVar(ncid,i-1);
        end
        %         % load data as netcdf
        %         if singleflag
        %             NC.Data.(nci.Variables(i).Name) = single(netcdf.getVar(ncid,i-1));
        %         elseif any(strcmpi(nci.Variables(i).Datatype,'string'))
        %             NC.Data.(nci.Variables(i).Name) = string(netcdf.getVar(ncid,i-1));
        %         else
        %             NC.Data.(nci.Variables(i).Name) = double(netcdf.getVar(ncid,i-1));
        %         end
        % try to deal with unsupported datatype problem...
    catch err
        warning(['Couldn''t load the variable "' nci.Variables(i).Name '", trying as HDF5...'])
        if strcmpi(err.identifier,'MATLAB:imagesci:netcdf:unrecognizedVarDatatype')
            % load as HDF5:
            NC.Data.(nci.Variables(i).Name) = h5read(filepath,['/' nci.Variables(i).Name '/']);
        else
            rethrow(err)
        end
    end
    
end

% %%%% Sometimes, all the data is in Groups first (like with VIIRS data).
% %%%% This means we have to specify the Group as well to get the data.
% if isempty(nci.Variables)
%     for i = 1:length(nci.Groups)
%         NC.Data.(nci.Groups(i).Name) = struct;
%     end
% end


%% Extract Data from within Groups
groupflag = 0;
% reding data from within a group doesn't seem easy using the basic
% netcdf.getVar() commands, so going to try using the higher level ncread()
% function. Chances are this will be a lot slower, but should work for now.
if ~isempty(nci.Groups)
    
    groupflag = 1;
    
    for i = 1:length(nci.Groups)
        
        % first, make the Group names in the Data structure:
        NC.Data.(nci.Groups(i).Name) = struct;
        
        % find the Group variables
        for j = 1:length(nci.Groups(i).Variables)
            
            numericflag = 1;
            if any(strcmpi(nci.Groups(i).Variables(j).Datatype,{'char','string'}))
                numericflag = 0;
            end
            
            % load data as netcdf
            if numericflag
            if singleflag
                NC.Data.(nci.Groups(i).Name).(nci.Groups(i).Variables(j).Name) ...
                    = single(ncread(filepath,['/' nci.Groups(i).Name '/' nci.Groups(i).Variables(j).Name ]));
            else
                NC.Data.(nci.Groups(i).Name).(nci.Groups(i).Variables(j).Name) ...
                    = double(ncread(filepath,['/' nci.Groups(i).Name '/' nci.Groups(i).Variables(j).Name ]));
            end
            else
                NC.Data.(nci.Groups(i).Name).(nci.Groups(i).Variables(j).Name) ...
                    = ncread(filepath,['/' nci.Groups(i).Name '/' nci.Groups(i).Variables(j).Name ]);
            end
            
        end
        
    end
    
end


%% Apply scale factors and offsets if they exist %%%%%%%%%%%%%%%%%%%%%%%%%%
% note these names are for the ECMWF netcdf formats, not clear if they'll
% work for other formats, but you've got the attributes anyway now so you
% can do it afterward.

% if ~groupflag % i have a problem here where the below doesn't index right if the data is in groups

for i = 1:length(nci.Variables)
    
    % First check that there are actually some attributes for this
    % variable...
    if ~isempty(nci.Variables(i).Attributes)
        
        % Fill values to NaN:
        if any(strcmpi({nci.Variables(i).Attributes.Name},'_FillValue'))
            ind = find(strcmpi({nci.Variables(i).Attributes.Name},'_FillValue'));
            fill_value = nci.Variables(i).Attributes(ind).Value;
            if isnumeric(fill_value)
                NC.Data.(nci.Variables(i).Name)(NC.Data.(nci.Variables(i).Name) == fill_value) = NaN;
            else
                % uh oh, the fill value is a char/string???????
%                 NC.Data.(nci.Variables(i).Name)(NC.Data.(nci.Variables(i).Name) == fill_value) = NaN;                
            end
        end
        % missing values to NaN:
        if any(strcmpi({nci.Variables(i).Attributes.Name},'missing_value'))
            ind = find(strcmpi({nci.Variables(i).Attributes.Name},'missing_value'));
            missing_value = nci.Variables(i).Attributes(ind).Value;
            NC.Data.(nci.Variables(i).Name)(NC.Data.(nci.Variables(i).Name) == missing_value) = NaN;
        end
        % if any of the variables have the term 'scale_factor':
        if any(strcmpi({nci.Variables(i).Attributes.Name},'scale_factor'))
            ind = find(strcmpi({nci.Variables(i).Attributes.Name},'scale_factor'));
            NC.Data.(nci.Variables(i).Name) = NC.Data.(nci.Variables(i).Name) .* double(nci.Variables(i).Attributes(ind).Value);
        end
        % if any of the variables have the term 'add_offset':
        if any(strcmpi({nci.Variables(i).Attributes.Name},'add_offset'))
            ind = find(strcmpi({nci.Variables(i).Attributes.Name},'add_offset'));
            NC.Data.(nci.Variables(i).Name) = NC.Data.(nci.Variables(i).Name) + double(nci.Variables(i).Attributes(ind).Value);
        end
        
    end
    
end

% end

%% Close netcdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

netcdf.close(ncid); clear ncid;

end







% if nargin == 1, varnums = 0:length(nci.Variables)-1; end;
% % varnums = 0:length(nci.Variables)-1;
% if nargin <= 2, attnums = 0:length(nci.Attributes)-1; end;
% % attnums = 0:length(nci.Attributes)-1;
%
% varnames = {nci.Variables.Name};
%
% % Cope with dashes in variable names:
% varnames = strrep(varnames,'-','_');
% % varnames = regexprep(varnames,'-','_');
%
%
% if ~isempty(nci.Attributes),
%     attnames = {nci.Attributes.Name};
% end
%
% ncid = netcdf.open(filepath,'nc_nowrite');
%
% for v = varnums,
%
%     nc.(varnames{v+1}) = double(netcdf.getVar(ncid,v));
%
% end
%
% for a = attnums,
%
%     nc.Attributes.(attnames{a+1}) = double(netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attnames{a+1}));
%
% end
%
%
% netcdf.close(ncid); clear ncid;
%
% % end


