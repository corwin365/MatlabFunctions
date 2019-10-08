
% Neil Hindley, University of Bath, UK. 2017.

% Read a netcdf file's variables, attributes and other into a structure

% EDIT 2018-05-24: Okay, getnet.m's been long overdue an overhaul, here
% I've added the functionality to read attrributes from individual
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


function NC = getnet(filepath,varargin)

% filepath = '/Volumes/NJMitchell-Scratch/Data/neil/e5_southgeorgia/sg_01.nc';

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

NC.Data = struct;

for i = 1:length(nci.Variables)
    if any(strcmpi(varargin,'single'))
        NC.Data.(nci.Variables(i).Name) = single(netcdf.getVar(ncid,i-1));
    else
        NC.Data.(nci.Variables(i).Name) = double(netcdf.getVar(ncid,i-1));
    end
        
end

%% Apply scale factors and offsets if they exist %%%%%%%%%%%%%%%%%%%%%%%%%%
% note these names are for the ECMWF netcdf formats, not clear if they'll
% work for other formats, but you've got the attributes anyway now so you
% can do it afterward.
for i = 1:length(nci.Variables)
    
    % Fill values to NaN:
    if any(strcmpi({nci.Variables(i).Attributes.Name},'_FillValue'))
        ind = find(strcmpi({nci.Variables(i).Attributes.Name},'_FillValue'));
        fill_value = nci.Variables(i).Attributes(ind).Value;
        NC.Data.(nci.Variables(i).Name)(NC.Data.(nci.Variables(i).Name) == fill_value) = NaN;
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
        NC.Data.(nci.Variables(i).Name) = NC.Data.(nci.Variables(i).Name) .* nci.Variables(i).Attributes(ind).Value;
    end
    % if any of the variables have the term 'add_offset':
    if any(strcmpi({nci.Variables(i).Attributes.Name},'add_offset'))
        ind = find(strcmpi({nci.Variables(i).Attributes.Name},'add_offset'));
        NC.Data.(nci.Variables(i).Name) = NC.Data.(nci.Variables(i).Name) + nci.Variables(i).Attributes(ind).Value;
    end
    
end

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


