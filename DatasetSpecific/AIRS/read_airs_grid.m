function buffer = read_airs_grid(filename,content_flag,content_list,gridname)

%   function buffer = read_airs_grid(filename,content_flag,content_list,gridname)
%
% Created by Stephen Licata (SL) on 03-17-2005.
% Updated on 12-05-2006 by SL.
%   Function name was changed to a standard form and documentation/comments were updated.
% Updated on 05-16-2007 by SL.
%   The content_flag = 4 option was added to dump the names of all data fields within a specific grid.
%
% DESCRIPTION:
% This function reads a Level 3 granule data file in the HDF-EOS grid format
% and extracts one or more data items from a grid structure within the data file.
% These data items could be actual data points or it could be characteristics
% about the data grid itself, such as grid name and coordinates.
%
% INPUT ARGUMENTS (REQUIRED)
%
% filename     - The fully qualified path to a Level 3 EOS-HDF grid format granule file.
%
% content_flag - An integer (0-4) that specifies the type of data to be extracted, as follows:
%                0: A cell array of the names of the grid(s) in that file.
%                1: The name and values of the grid's dimension parameters.
%                2: The name and values of the grid's attribute parameters.
%                3: The name and values of the grid's data field parameters.
%                4: A cell array of the names of all data fields in a specific grid.
%
% INPUT ARGUMENTS [OPTIONAL]:
%
% INPUT OPTIONS:
% content_list - An array of names for the content items to be returned. If left unspecified,
%                the function will return ALL the parameters in that content category.
%
% gridname     - A text expression for the data grid within the granule file that is to be
%                examined. This function will only process one data grid at a time. In the
%                (typical) case that there is only ONE data grid in the granule file, this
%                argument can be left unspecified.
%
% RETURN VALUES:
%
% buffer       - MATLAB data structure whose content is based on "content_flag".
%                Option 0: buffer will be a cell array of all the grid names within a file.
%                Option 4: buffer will be a cell array of all the parameter names within a specified grid.
%                Options 1-3: buffer will be a generic MATLAB data structure in which each member
%                has a "name" and a "data" component, which are actual data values.
%
% SIDE EFFECTS
%
%              Parameter names in the input file containing with a period character will be saved
%              in the output data structure (buffer) with an underscore in place of the period.
%*****************************************************************************************************

   prog_name = 'read_airs_grid';

   buffer        = [];
   field_list    = {};

   if ~exist('gridname','var')
      gridname = [];
   end

   if ~exist('content_list','var')
      content_list = {};
   end

   if ischar(content_list)
      tmp = content_list;
      content_list = {tmp};
      clear tmp
   end

% Abort the program if no data file has been provided.
   if isempty(filename)
      disp([prog_name ': ERROR - No input filename was specified.'])
      return
   end

% Abort the program if no data type has been specified.
   if isempty(content_flag)
      disp([prog_name ': ERROR - No content code (type) was provided.'])
      return
   end

% Set up a set of names (labels) for the different types of data queries.
   type_list = {'grid','dimension','attribute','field'};

% Get a file id value.
   fid   = hdfgd('open',filename,'read');

% Abort the program if the file does not exist or cannot be read.
   if fid == -1
      disp([prog_name ': ERROR - ' filename ' could not be opened.'])
      status = hdfgd('close', fid);
      return
   end

% Get the number of data grids in the file (normally just 1)
   [ngrid,gridlist] = hdfgd('inqgrid',filename);
   gridlist_org     = gridlist;

% Abort the program if no data grid(s) is/are found.
   if ngrid == 0 | isempty(gridlist)
      disp([prog_name ': ERROR - ' filename ' has no data grid.'])
      status = hdfgd('close', fid);
      return;
   end

% Convert the gridlist to an array.
   for i=1:ngrid
      [t,gridlist]  = strtok(gridlist,',');
      grid_list{i}  = t;
   end

% If only the grid list is requested, return that array and end the program.
   if content_flag == 0
      buffer   = grid_list;
      status   = hdfsw('close', fid);
      return;
   end

% Accept a user-specified data grid name or the first (and only) available grid within the file.
% Only continue processing if the data set is confined to a single grid.
% Abort the program if more than one data grid can be extracted from 'gridname'.
   if ~isempty(gridname)
      gridname = gridname;
   else
      if ngrid == 1
         gridname  = grid_list{1};
      else
         disp([prog_name ': ERROR - This file contains multiple data grids - Choose just one.'])
         disp([prog_name ': Grid List = ' gridlist_org])
         status = hdfgd('close', fid);
         return;
      end
   end

% Attach an ID to this grid.
   grid_id = hdfgd('attach', fid, gridname);
   if (grid_id == -1)
      disp([prog_name ': ERROR - Unable to open grid ' gridname])
      status = hdfgd('detach',grid_id);
      status = hdfgd('close', fid);
      return;
   end

% If only the list of parameters name within a grid is requested, provide it then end the program.
   [nfield,fieldlist,rank,ntype] = hdfgd('inqfields',grid_id);
   if nfield == 0
      disp([prog_name ': ERROR - No ' gridname ' data fields were located.'])
      status = hdfgd('detach',grid_id);
      status = hdfgd('close', fid);
      return;
   end 
   for i=1:nfield
      [t,fieldlist]  = strtok(fieldlist,',');
      field_list{i}  = t;
   end

% If only the data field list is need, return it and end the program.
   if (content_flag == 4) 
      buffer   = field_list;
      status   = hdfgd('detach',grid_id);
      status   = hdfgd('close', fid);
      return;
   end


% Select the content list (i.e., set of parameter names) from the grid itself
% unless these names are specified directly by the user.
   if isempty(content_list)
      switch content_flag
         case {1}
            [ndim,dimlist,dims] = hdfgd('inqdims',grid_id);
            if ndim == 0
               disp([prog_name ': ERROR - The ' gridname ' grid has no data dimensions.'])
               status = hdfgd('detach',grid_id);
               status = hdfgd('close', fid);
               return;
            end
            for i=1:ndim
            [t,dimlist]  = strtok(dimlist,',');
               content_list{i}  = t;
            end
         case {2}
            [nattrib,attrlist]    = hdfgd('inqattrs',grid_id);
            if nattrib == 0
               disp([prog_name ': ERROR - The ' gridname ' grid has no attributes.'])
               status = hdfgd('detach',grid_id);
               status = hdfgd('close', fid);
               return;
            end
            for i=1:nattrib
               [t,attrlist]  = strtok(attrlist,',');
               content_list{i}  = t;
            end
         case {3}
            content_list  = field_list;
         otherwise
            disp([prog_name ': ERROR - No content list (based on content flag) was generated.'])
            status = hdfgd('detach',grid_id);
            status = hdfgd('close', fid);
            return;
      end
   end

% Abort the program if the content list is just a single blank string entry.
   if length(content_list) == 1 & length(content_list{1}) == 0
      disp([prog_name ': ERROR - No set of ' type_list{content_flag} ' parameter names was specified.'])
      status = hdfgd('detach',grid_id);
      status = hdfgd('close', fid);
      return
   end

% Abort the program if the content list is still undefined.
   if isempty(content_list)
      disp([prog_name ': ERROR - No set of ' type_list{content_flag} ' parameter names was specified.'])
      status = hdfgd('detach',grid_id);
      status = hdfgd('close', fid);
      return
   end

   buffer = [];

   for i = 1:length(content_list)

      item_name = content_list{i};

% Discard any parameter names that have a '=' sign in the name.
% NOTE: This is an optional feature based on experience with these data files.
      if ~isempty(findstr(item_name,'=')) 
         continue;
      end

% Extract the data value based on the parameter name and data type.
      switch content_flag
         case {1}
            [item_val     ] = hdfgd('diminfo',grid_id,item_name);
         case {2}
            [item_val,fail] = hdfgd('readattr',grid_id,item_name);
         case {3}
            [item_val,fail] = hdfgd('readfield',grid_id,item_name,[],[],[]);
         otherwise
            disp([prog_name ': ERROR - No valid content flag was specified.'])
            status = hdfgd('detach',grid_id);
            status = hdfgd('close', fid);
            return;
      end

      if ~isempty(item_val) 

% Now replace any '.' characters in the parameter name with '_'.
         if ~isempty(findstr(item_name,'.'))
            new_name  = strrep(item_name,'.','_');
            item_name = new_name;
         end

% Now add that parameter data set as a field to the buffer structure.
         if ~isfield(buffer,item_name)
            eval(['buffer.' item_name ' = item_val;'])
         end

      end

   end

   status = hdfgd('detach',grid_id);
   status = hdfgd('close', fid);

