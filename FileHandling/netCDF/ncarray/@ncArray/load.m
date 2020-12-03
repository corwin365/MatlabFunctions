% [val,coord1,coord2,...] = load(self,'coord_name1',range1,'coord_name2',range2,...)
% Load a subset of a variable based on range of coordiante variables.
% The names of the coordinates (coord_name1, coord_name2,...) coorespond to the standard_name attribute.
% Only 1-dimensional coordinates are currently supported.
%
%
% Example
% [temp,lon,lat,depth,time] = load(self,'longitude',[0 10],'latitude',[0 10])

function varargout = load(self,varargin)

c = coord(self);

for i = 1:length(c)
  c(i).v = full(c(i).val);  
  % per default take all data along a dimension
  c(i).index = ':';
  c(i).sub = c(i).v;
end

% loop over all constraints
for i = 1:2:length(varargin) 
  name = varargin{i};
  
  j = find(strcmp(name,{c.standard_name}));
  if isempty(j)
    warning(['no coordinate has the standard_name ' name ...
             '. Try to use variable names.']);
  
    j = find(strcmp(name,{c.name}));
    if isempty(j)
      error(['no coordinate has the name ' name '.']);
    end  
  end
    
  range = varargin{i+1};
  
  if numel(range) == 1
    dist = abs(c(j).v - range);
    [mindist,i] = min(dist);
    %mindist
  else    
    i = find(range(1) < c(j).v & c(j).v < range(end)); 
    i = min(i):max(i);
  end
  
  c(j).index = i;
  c(j).sub = c(j).v(i);
end

idx = substruct('()',{c.index});
data = subsref (self,idx);

varargout = {data,c.sub};


% i = 1;

% mask = xr(1) <= x & x <= xr(2);
% l = find(mask);

% [ij{:}] = ind2sub(size(mask),l);

% for j=1:len
% mins(j) = min(ij{j});
% maxs(j) = max(ij{j});
 
