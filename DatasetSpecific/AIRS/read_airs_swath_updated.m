function buffer = read_airs_swath_updated(filename,content_flag,content_list,swathname)

%   function buffer = read_airs_swath(filename,content_flag,content_list,swathname)
%
% Created by Stephen Licata (SL) on 03-17-2005.
% Updated on 12-05-2006.
%   Function name was changed to a standard form and documentation/comments were updated.
% Updated on 05-18-2006 by SL.
%   The content_flag = 4 option was added to enable a dump of the swath's data fields.
%
% DESCRIPTION:
% This function reads a Level 1 or 2 granule data file in the HDF-EOS format
% and extracts one of the following information types into a data buffer:
%
% INPUT ARGUMENTS (REQUIRED)
%
% filename  - The fully qualified path to a Level 1/2 EOS-HDF swath format granule file.
%
% content_flag - An integer (0-4) that specifies the type of data to be extracted, as follows:
%                0: A cell array of the names of all swaths in the specified file.
%                1: The name and values of the swath's dimension parameters.
%                2: The name and values of the swath's attribute parameters.
%                3: The name and values of the swath's data field parameters.
%                4: A cell array of the names of all parameters in a specific swath.
%
% INPUT ARGUMENTS [OPTIONAL]:
%
% content_list - An array of names for the content items to be returned. If left unspecified,
%                the function will return ALL the parameters in that content category.
%
% swathname    - A text expression for the data swath within the granule file that is to be
%                examined. This function will only process one data swath at a time. In the
%                (typical) case that there is only ONE data swath in the granule file, this
%                argument can be left unspecified.
%
% RETURN VALUES:
%
% buffer       - MATLAB data structure whose content is based on "content_flag". 
%                Option 0: buffer will be a cell array of all the swath names within a file. 
%                Option 4: buffer will be a cell array of all the parameter names within a specified swath. 
%                Options 1-3: buffer will be a generic MATLAB data structure in which each member
%                has a "name" and a "data" component, which are actual data values.
%
% SIDE EFFECTS:
%
%              Options 3 and 4, by default, retrieve NOTH the geo-location AND science data fields.
%
%              Parameter names in the input file containing with a period character will be saved
%              in the output data structure (buffer) with an underscore in place of the period.
%*****************************************************************************************************

import matlab.io.hdfeos.*

   prog_name = 'read_airs_swath';

   buffer = [];

   if ~exist('swathname','var')
      swathname = [];
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
   type_list = {'swath','dimension','attribute','field','field'};

% Get a file id value.
   fid   = sw.attach('open',filename,'read');

% Abort the program if the file does not exist or cannot be read.
   if fid == -1
      disp([prog_name ': ERROR - ' filename ' could not be opened.'])
      status = hdfsw('close', fid);
      return
   end

% Get the number of data swaths in the file (normally just 1)
   [nswath, swathlist] = hdfsw('inqswath',filename);
   swathlist_org       = swathlist;

% Abort the program if no data swath(s) is/are found.
   if nswath == 0 | isempty(swathlist)
      disp([prog_name ': ERROR - ' filename ' has no data swath.'])
      status = hdfsw('close', fid);
      return;
   end

% Get the name of one or more available swaths.
   for i=1:nswath
      [t,swathlist]  = strtok(swathlist,',');
      swath_list{i}  = t;
   end

% If only the swath list is requested, return that information and end the program.
   if content_flag == 0
      buffer   = swath_list;
      status   = hdfsw('close', fid);
      return;
   end

% Accept a user-specified data swath name or the one and only available swath name.
% Abort the program if more than one data swath has been extracted from 'swathlist'.
   if ~isempty(swathname)
      swathname = swathname;
   else
      if nswath > 1
         disp([prog_name ': ERROR - This file contains multiple data swaths - Choose just one.'])
         disp([prog_name ': Swath List = ' swathlist_org])
         status = hdfsw('close', fid);
         return;
      else
         swathname = swath_list{1};
      end
   end

% Attach an ID to this swath.
   swath_id = hdfsw('attach', fid, swathname);
   if (swath_id == -1)
      disp([prog_name ': ERROR - Unable to open swath ' swathname])
      status = hdfsw('detach',swath_id);
      status = hdfsw('close', fid);
      return;
   end

% Get both the geo-location and science parameters in one long list.
   [ngfields,gfieldlist,rank,ntype] = hdfsw('inqgeofields',swath_id);
   [nfields,fieldlist,rank,type]    = hdfsw('inqdatafields',swath_id);

   if ((nfields + ngfields) == 0)
      disp([prog_name ': ERROR - The ' swathname ' swath has no data fields.'])
      status = hdfgd('detach',swath_id);
      status = hdfgd('close', fid);
      return;

   else
      tmp    = ' ';
      if ~isempty(gfieldlist)
         tmp      = gfieldlist;
      end
      if ~isempty(fieldlist)
         tmp      = [tmp ',' fieldlist];
      end
      tmp       = strrep(tmp,',,',',');
      tmp       = strrep(tmp,' ','');
      for i=1:(ngfields+nfields)  
         [t,tmp]       = strtok(tmp,','); 
         field_list{i} = t; 
      end
   end

% Stop here if all we need is the list of data fields.
   if (content_flag == 4)
      buffer     = field_list;
      status     = hdfsw('detach',swath_id);
      status     = hdfsw('close', fid);
      return;
   end

% Select the content list (i.e., set of parameter names) from the swath itself
% unless these names are specified directly by the user.
   if isempty(content_list)
      switch content_flag
         case {1}
            [ndim,tmp,dims] = hdfsw('inqdims',swath_id);
            if (ndim < 1) 
               disp([prog_name ': ERROR - The ' swathname ' swath has no dimension information.'])
               status = hdfsw('detach',swath_id);
               status = hdfsw('close', fid);
               return;
            end
            for i=1:ndim  
               [t,tmp]         = strtok(tmp,','); 
               content_list{i} = t; 
            end
         case {2}
            [nattrib,tmp]    = hdfsw('inqattrs', swath_id);
            if (nattrib < 1) 
               disp([prog_name ': ERROR - The ' swathname ' swath has no attribute information.'])
               status = hdfsw('detach',swath_id);
               status = hdfsw('close', fid);
               return;
            end
            for i=1:nattrib  
               [t,tmp]         = strtok(tmp,','); 
               content_list{i} = t; 
            end
         case {3}
            content_list    = field_list;
         otherwise
            disp([prog_name ': ERROR - No content list (based on content flag) was generated.'])
            status = hdfsw('detach',swath_id);
            status = hdfsw('close', fid);
            return;
      end
   end

% Abort the program if the content list is just a single blank string entry.
   if length(content_list) == 1 & length(content_list{1}) == 0
      disp([prog_name ': ERROR - No set of ' type_list{content_flag} ' parameter names was specified.'])
      status = hdfsw('detach',swath_id);
      status = hdfsw('close', fid);
      return
   end

% Abort the program if the content list is still undefined.
   if isempty(content_list)
      disp([prog_name ': ERROR - No set of ' type_list{content_flag} ' parameter names was specified.'])
      status = hdfsw('detach',swath_id);
      status = hdfsw('close', fid);
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
            [item_val     ] = hdfsw('diminfo',swath_id,item_name);
         case {2}
            [item_val,fail] = hdfsw('readattr',swath_id,item_name);
         case {3}
            [item_val,fail] = hdfsw('readfield',swath_id,item_name,[],[],[]);
         otherwise
            disp([prog_name ': ERROR - No valid content flag was specified.'])
            status = hdfsw('detach',swath_id);
            status = hdfsw('close', fid);
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

% Close out the swath and file.
   status = hdfsw('detach',swath_id);
   status = hdfsw('close', fid);

