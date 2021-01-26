function ZG = bin2matN(NDims,varargin)

%calls the appropriate version of bin2mat for 1d, 2d or 3d data
%in future, extend to ND, but this will do for now

switch NDims
  case 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %call bin2mat1d
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x    = double(varargin{1});
    z    = double(varargin{2});
    xi   = double(varargin{3});
    
    if numel(varargin) < 4;
      ZG = bin2mat1d(x,z,xi);
    else
      ZG = bin2mat1d(x,z,xi,varargin{4:end});
    end
    
  case 2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bin2mat2d
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x    = double(varargin{1});
    y    = double(varargin{2});
    z    = double(varargin{3});
    xi   = double(varargin{4});
    yi   = double(varargin{5});
    
    if numel(varargin) < 6;
      ZG = bin2mat2d(x,y,z,xi,yi);
    else
      ZG = bin2mat2d(x,y,z,xi,yi,varargin{6:end});
    end
    
  case 3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bin2mat3d
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x    = double(varargin{1});
    y    = double(varargin{2});
    z    = double(varargin{3});
    v    = double(varargin{4});
    xi   = double(varargin{5});
    yi   = double(varargin{6});
    zi   = double(varargin{7});
    
    if numel(varargin) < 8;
      ZG = bin2mat3d(x,y,z,v,xi,yi,zi);
    else
      ZG = bin2mat3d(x,y,z,v,xi,yi,zi,varargin{8:end});
    end
  
  case 3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bin2mat4d
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x    = double(varargin{1});
    y    = double(varargin{2});
    z    = double(varargin{3});
    a    = double(varargin{4});    
    v    = double(varargin{5});
    xi   = double(varargin{6});
    yi   = double(varargin{7});
    zi   = double(varargin{8});
    ai   = double(varargin{9});    
    
    if numel(varargin) < 10;
      ZG = bin2mat3d(x,y,z,a,v,xi,yi,zi,ai);
    else
      ZG = bin2mat3d(x,y,z,a,v,xi,yi,zi,ai,varargin{10:end});
    end
    
    
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin2mat1d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ZG = bin2mat1d(x,z,xi,varargin)
% BIN2MAT - create a matrix from scattered data without interpolation


%check inputs
error(nargchk(3,inf,nargin,'struct'));

%make sure the vectors are column vectors
x = x(:);
z = z(:);

if all(any(diff(cellfun(@length,{x,z}))));
  error('Inputs x and z must be the same size');
end

%process optional input
fun=@mean;
test=1;
if ~isempty(varargin)
  fun=varargin{1};
  if ~isa(fun,'function_handle');
    fun=str2func(fun);
  end
  
  %test the function for non-scalar output
  test = feval(fun,rand(5,1),varargin{2:end});
  
end

%grid nodes
m=numel(xi);

%limit values to those within the specified grid
xmin=min(xi);
xmax=max(xi);

gind =(x>=xmin & x<=xmax);

%find the indices for each x and y in the grid
[junk,xind] = histc(x(gind),xi);


%break the data into a cell for each grid node
blc_ind=accumarray([xind],z(gind),[m 1],@(x){x},{NaN});

%evaluate the data in each grid using FUN
if numel(test)>1
  ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind,'uni',0);
else
  ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin2mat2d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function ZG = bin2mat2d(x,y,z,XI,YI,varargin)
% BIN2MAT - create a matrix from scattered data without interpolation
%
%   ZG = BIN2MAT(X,Y,Z,XI,YI) - creates a grid from the data
%   in the (usually) nonuniformily-spaced vectors (x,y,z)
%   using grid-cell averaging (no interpolation). The grid
%   dimensions are specified by the uniformily spaced vectors
%   XI and YI (as produced by meshgrid).
%
%   ZG = BIN2MAT(...,@FUN) - evaluates the function FUN for each
%   cell in the specified grid (rather than using the default
%   function, mean). If the function FUN returns non-scalar output,
%   the output ZG will be a cell array.
%
%   ZG = BIN2MAT(...,@FUN,ARG1,ARG2,...) provides aditional
%   arguments which are passed to the function FUN.
%
%   EXAMPLE
%
%   %generate some scattered data
%    [x,y,z]=peaks(150);
%    ind=(rand(size(x))>0.9);
%    xs=x(ind); ys=y(ind); zs=z(ind);
%
%   %create a grid, use lower resolution if
%   %no gaps are desired
%    xi=min(xs):0.25:max(xs);
%    yi=min(ys):0.25:max(ys);
%    [XI,YI]=meshgrid(xi,yi);
%
%   %calculate the mean and standard deviation
%   %for each grid-cell using bin2mat
%    Zm=bin2mat(xs,ys,zs,XI,YI); %mean
%    Zs=bin2mat(xs,ys,zs,XI,YI,@std); %std
%
%   %plot the results
%    figure
%    subplot(1,3,1);
%    scatter(xs,ys,10,zs,'filled')
%    axis image
%    title('Scatter Data')
%
%    subplot(1,3,2);
%    pcolor(XI,YI,Zm)
%    shading flat
%    axis image
%    title('Grid-cell Average')
%
%    subplot(1,3,3);
%    pcolor(XI,YI,Zs)
%    shading flat
%    axis image
%    title('Grid-cell Std. Dev.')
%
% SEE also RESHAPE ACCUMARRAY FEVAL

% A. Stevens 3/10/2009
% astevens@usgs.gov

%check inputs
error(nargchk(5,inf,nargin,'struct'));

%make sure the vectors are column vectors
x = x(:);
y = y(:);
z = z(:);

if all(any(diff(cellfun(@length,{x,y,z}))));
  error('Inputs x, y, and z must be the same size');
end

%process optional input
fun=@mean;
test=1;
if ~isempty(varargin)
  fun=varargin{1};
  if ~isa(fun,'function_handle');
    fun=str2func(fun);
  end
  
  %test the function for non-scalar output
  test = feval(fun,rand(5,1),varargin{2:end});
  
end

%grid nodes
xi=XI(1,:);
yi=YI(:,1);
[m,n]=size(XI);

%limit values to those within the specified grid
xmin=min(xi);
xmax=max(xi);
ymin=min(yi);
ymax=max(yi);

gind =(x>=xmin & x<=xmax & ...
  y>=ymin & y<=ymax);

%find the indices for each x and y in the grid
[junk,xind] = histc(x(gind),xi);
[junk,yind] = histc(y(gind),yi);

%break the data into a cell for each grid node
blc_ind=accumarray([yind xind],z(gind),[m n],@(x){x},{NaN});

%evaluate the data in each grid using FUN
if numel(test)>1
  ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind,'uni',0);
else
  ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind);
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin2mat3d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ZG = bin2mat3d(x,y,c,z,XI,YI,CI,varargin)

%BIN2MAT3D
% Bins data (Z) in 3D coordinates (x,y,c)  nto 3D bins specified by XI, YI,
% CI

% A. Stevens 3/10/2009
% astevens@usgs.gov

%check inputs
error(nargchk(7,inf,nargin,'struct'));

%make sure the vectors are column vectors
x = x(:);
y = y(:);
c = c(:);
z = z(:);

if all(any(diff(cellfun(@length,{x,y,z,c}))));
  error('Inputs x, y,c and z must be the same size');
end

%process optional input
fun=@mean;
test=1;
if ~isempty(varargin)
  fun=varargin{1};
  if ~isa(fun,'function_handle');
    fun=str2func(fun);
  end
  
  %test the function for non-scalar output
  % test = feval(fun,rand(5,1),varargin{2:end});
  
end

%grid nodes
xi=squeeze(XI(1,:,1));
yi=squeeze(YI(:,1,1));
ci=squeeze(CI(1,1,:));
[m,n,l]=size(XI);

%limit values to those within the specified grid
xmin=min(xi);
xmax=max(xi);
ymin=min(yi);
ymax=max(yi);
cmin=min(ci);
cmax=max(ci);
gind =(x>=xmin & x<=xmax & ...
  y>=ymin & y<=ymax & c>=cmin & c<=cmax );

%find the indices for each x and y in the grid
[junk,xind] = histc(x(gind),xi);
[junk,yind] = histc(y(gind),yi);
[junk,cind] = histc(c(gind),ci);

%break the data into a cell for each grid node
blc_ind=accumarray([yind xind cind],z(gind),[m n l],@(x){x},{NaN});

%evaluate the data in each grid using FUN
if numel(test)>1
  ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind,'uni',0);
else
  ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind);
end
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin2mat4d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ZG = bin2mat4d(x,y,c,d,z,XI,YI,CI,DI,varargin)

%BIN2MAT4D
% Bins data (Z) in 4D coordinates (x,y,c)  into 4D bins specified by XI, YI,
% CI

%check inputs
error(nargchk(9,inf,nargin,'struct'));

%make sure the vectors are column vectors
x = x(:);
y = y(:);
c = c(:);
d = d(:);
z = z(:);

if all(any(diff(cellfun(@length,{x,y,z,c,d}))));
  error('Inputs x, y, c, d and z must be the same size');
end

%process optional input
fun=@mean;
test=1;
if ~isempty(varargin)
  fun=varargin{1};
  if ~isa(fun,'function_handle');
    fun=str2func(fun);
  end
  
  %test the function for non-scalar output
  % test = feval(fun,rand(5,1),varargin{2:end});
  
end

%grid nodes
xi=squeeze(XI(1,:,1));
yi=squeeze(YI(:,1,1));
ci=squeeze(CI(1,1,:));
di=squeeze(DI(1,1,:));
[m,n,l,o]=size(XI);

%limit values to those within the specified grid
xmin=min(xi);
xmax=max(xi);
ymin=min(yi);
ymax=max(yi);
cmin=min(ci);
cmax=max(ci);
dmin=min(di);
dmax=max(di);
gind =(x>=xmin & x<=xmax & ...
       y>=ymin & y<=ymax & ...
       c>=cmin & c<=cmax & ...
       d>=dmin & d<=dmax );

%find the indices for each x and y in the grid
[junk,xind] = histc(x(gind),xi);
[junk,yind] = histc(y(gind),yi);
[junk,cind] = histc(c(gind),ci);
[junk,dind] = histc(d(gind),di);

%break the data into a cell for each grid node
blc_ind=accumarray([yind xind cind,dind],z(gind),[m n l o],@(x){x},{NaN});

%evaluate the data in each grid using FUN
if numel(test)>1
  ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind,'uni',0);
else
  ZG=cellfun(@(x)(feval(fun,x,varargin{2:end})),blc_ind);
end
end

