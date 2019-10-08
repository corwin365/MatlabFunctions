function A=get_MLS(file, swath, READATT)
%
%function A=readL2GP(FILE, SWATH)
%function A=readL2GP(DIRECTORY, SWATH)
%
% eg
%FILE  = '/data/emls/l2gp/v02.20/2005/028/MLS-Aura_L2GP-O3_v02-20-c06_2005d243.he5'
%SWATH = 'O3'
%  or
%FILE  = '/data/emls/l2gp/v02.20/2005/028/MLS-Aura_L2GP-DGG_v02-020-c06_2005d243.he5'
%SWATH = 'O3-640'   %(a diagnostic product stored in the DGG file)
%  or
%DIRECTORY = '/data/emls/l2gp/v02.20/2005/028/'
%SWATH     = 'O3-640'  % (a diagnostic product stored in the DGG file)
%
%If the first argument is a directory, the program attempts to guess the
%correct file, first looking for files ending with '.he5' that have SWATH in their
%names and, if none are found, then for those that have 'DGG' in their
%names.  In either case, if multiple files are found, the last (alphabetically) is
%used.  This selects the latest "cycle" (eg 'c06', above) if level-2 has been run repeatedly.
%
%READATT = true [false] 
%  true causes attributes to be returned in a
%  substructure, "A.att.(attribute_name)"
%
%the structure field 'L2gpValue' is the product and 'L2gpPrecision' is the precision
%
%
%Note: some fields are single or int32 and may need to be caste to double
% (eg   double(A.Status) ) for compatibility with certain matlab operations
%
% Michael J. Schwartz  04-Oct-2006
%
%Copyright 2006, by the California Institute of Technology. ALL RIGHTS RESERVED. 
%United States Government Sponsorship acknowledged. 
%Any commercial use must be negotiated with the Office of Technology Transfer 
%at the California Institute of Technology.
% 
%This software is subject to U.S. export control laws and regulations and has been 
%classified as 4D994.  By accepting this software, the user agrees to comply with all 
%applicable U.S. export laws and regulations.  User has the responsibility to obtain 
%export licenses, or other export authority as may be required before exporting such 
%information to foreign countries or providing access to foreign persons.



if exist('READATT') ~=1, READATT=false; end
foo=dir(file);

if isempty(foo), error(['file ' file ' not found.']); end

if isdir(file),
    Path=fileparts([file '/']);
else
    Path=fileparts(file);
end

for i=1:length(foo),
    [Pathstr, Name, Ext]=fileparts(foo(i).name);
    if strcmp(Ext, '.met'), r(i)=0; else r(i)=1; end
end
foo=foo(find(r));  %get rid of .met extensions

if length(foo)==1, 
    file=foo.name;
else
    a=strfind({foo.name}, swath);
    for i=1:length(a), aa(i)=~isempty(a{i}); end
    aa=find(aa);
    if isempty(aa),
        a=strfind({foo.name}, 'DGG');
        for i=1:length(a), aa(i)=~isempty(a{i}); end
        aa=find(aa);
        if ~isempty(aa),
            file=foo(aa(end)).name;
        else
            error(['file ' file ' not found.']);
        end
    else
        file=foo(aa(end)).name;
    end
end;

b=file;
for i=1:3,
    [a,b]=strtok(b, '_');
end

r=find(a=='-');
A.version=a(1:(r(end)-1));  %parse 'version' from the file name 
A.l2cycle=a((r(end)+1):end);  %parse 'l2cycle' from the file name 

[a,b]=strtok(b, '_');
A.date=strtok(a, '.');



file=fullfile(Path, file);
if strcmp(file(end+(-3:0)), '.he5'),

    A.file=file;
    A.swath=swath;

    if (0)
        fi=hdf5info(file, 'ReadAttributes', false);
        for i=1:length(fi.GroupHierarchy.Groups(1).Groups(2).Groups),
            SWATH{i}=strtokn(fi.GroupHierarchy.Groups(1).Groups(2).Groups(i).Name, 'last', '/');
        end
        N=find(strcmp(SWATH, swath));
        if isempty(N), error(['Can not find swath ' swath ' in ' file '.']); end
        datafields_=    {fi.GroupHierarchy.Groups(1).Groups(2).Groups(N).Groups(1).Datasets.Name};
        geofields_=    {fi.GroupHierarchy.Groups(1).Groups(2).Groups(N).Groups(2).Datasets.Name};

        for i=1:length(datafields_),
            A.(strtokn(datafields_{i}, 'last', '/'))=hdf5read(A.file, datafields_{i});
        end
        for i=1:length(geofields_),
            A.(strtokn(geofields_{i}, 'last', '/'))=hdf5read(A.file, geofields_{i});
        end

    else
        datafields={'L2gpValue' 'L2gpPrecision' 'Quality' 'Status', 'Convergence'};
        geofields= {'ChunkNumber' 'Latitude' 'LineOfSightAngle' 'LocalSolarTime' 'Longitude' 'OrbitGeodeticAngle' 'Pressure' 'SolarZenithAngle' 'Time'};

        for i=1:length(datafields),
            try
                A.(datafields{i})=hdf5read(A.file, ['/HDFEOS/SWATHS/' swath '/Data Fields/' datafields{i}]);
            catch
                A.(datafields{i})=[];
            end
        end
        for i=1:length(geofields),
            try
                A.(geofields{i})=hdf5read(A.file, ['/HDFEOS/SWATHS/' swath '/Geolocation Fields/' geofields{i}]);
            catch
                A.(geofields{i})=[];
            end
        end
        if READATT, fi=hdf5info(file, 'ReadAttributes', false); end
    end




    if READATT,
        Att=fi.GroupHierarchy.Groups(1).Groups(1).Groups.Attributes;
        r=find(Att(1).Name=='/');
        for i=1:length(Att),
            name=strrep(Att(i).Name((r(end)+1):end), ' ', '_');
            %    A.(name)=hdf5read(file, Att(i).Name);
            %    if strcmp(class(A.(name)), 'hdf5.h5string'), A.(name)=A.(name).Data; end
            A.att.(name)=hdf5read(file, Att(i).Name);
            if strcmp(class(A.att.(name)), 'hdf5.h5string'), A.att.(name)=A.att.(name).Data; end
        end
    else
        A.att=struct([]);  %included so that structure arrays may be created including elements with and without
    end




else
    warning([file '  does not end .he5'])
end

%if isfield(A, 'PCF1'), A=rmfield(A, 'PCF1'); end
%if isfield(A, 'PCF2'), A=rmfield(A, 'PCF2'); end

