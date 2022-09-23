
% nph_sinefit.m
%
% Created based on a really simple sine-fitting method I saw online, simply
% using the \ operator. It's witchcraft, but it seems to work.

% ALSO: added the ability to add weightings for the fit (turns out this was
% quite hard and involved lots of matrix maths):
% [yfit,F] = nph_sinefit(x,y,period_estimate,'weights',w);
%
% EDIT: Added the ability to fit multiple periods at once. Simply set
% period_estimate to be a vector of the periods you want to fit.
% Tis results in an output matrix that has all the factors Acos() and Bsin()
% for each period estimate, i.e. Acos(t1) + Bsin(t1) + Ccos(t2) + Dsin(t2) + ... + E.
% This has meant I've had to change to output structure to make it easier to
% reconstruct the signal afterwards. The new structure looks like, for a
% 3-period fit for example:
% F = [ A     B     0
%       C     D     0
%       E     F     G]; where G is the signal mean part.
%
% For a single-period fit, the output still looks like:
% F = [ A     B     C]; where C is the signal mean part.
%
% Have fun!
%
%
%
% INPUTS:
%
% x - the coordinates of the corresponding y values. This allows for
% irregularly-spaced data!
% y - vector of signal data at the points in x to be sine-fitted.
% period_estimate - estimate of period to which to fit a sinusoid, in
% whatever units x was in.
%
% OUTPUTS:
%
% OUT - the fitted sinusoidal wave, equal to the length of IN and evaluated
% at the same points specified in IN.
%
% F - 3-element vector of coefficients [A B C], where the fit is defined
% as:               yfit = C + A*cos(2PIx/t) + B*sin(2PIx/t);
% where C is obv the mean "offset" of the signal (DC) and t is the estimated period.
%
% If you're interested, phase is computed as phi = atan2(A,B), or atan2(F(1),F(2))
%
% [yfit,F] = nph_sinefit(x,y,period_estimate);
%

function varargout = nph_sinefit(x,IN,period_estimate,varargin)

if any(strcmpi(varargin,'weights'))
    method = 'weights';
else
    %     if any(strcmpi(varargin,'weights2'))
    %     method = 'weights2';
    %     else
    method = 'backslash';
    %     end
end

% first, linearise to column vectors
origsize = size(IN);
IN = IN(:); x = x(:); period_estimate = period_estimate(:);

% then cope with NaNs by just ignoring them:
nanlocs = isnan(x) | isnan(IN);
% finiteinds = intersect(find(~isnan(x)),find(~isnan(IN)));
xfinite = x(~nanlocs); % make a list of x values that are non-nan
y = IN(~nanlocs);


switch method
    % BACKSLASH METHOD ====================================================
    case 'backslash'
        
        % % % % %         % create fit matrix:
        % % % % %         X = ones(length(y),3); % three part fit, FIT = A*cos(2*pi*t/period) + B*sin(2*pi*t/period) + C;
        % % % % %         X(:,1) = cos((2*pi*xfinite) ./ period_estimate);
        % % % % %         X(:,2) = sin((2*pi*xfinite) ./ period_estimate);
        
        % NEW: allow for simultaneous multi-component fit, like so:
        % A*cos(2*pi*t/periods(1)) + B*sin(2*pi*t/periods(1)) ...
        % C*cos(2*pi*t/periods(2)) + D*sin(2*pi*t/periods(2)) ...
        % ... + E;
        
        % why not use the outer product to assemble the sines and
        % cosines... (nice!)
        cosparts = cos((2*pi*xfinite) * (1./period_estimate)');
        sinparts = sin((2*pi*xfinite) * (1./period_estimate)');
        oneparts = ones(size(xfinite));
        
        % assemble them:
        X = ones(length(y),(2*length(period_estimate))+1);
        sincos = 1;
        for i = 1:2:(2*length(period_estimate))
            X(:,i)   = cosparts(:,sincos);
            X(:,i+1) = sinparts(:,sincos);
            sincos = sincos+1;
        end
        X(:,end) = oneparts;
        
        % collapse the fit using matrix inversion (the witchcraft)
        F = X \ y;
        
        % % % % % %         THIS IS SLOWWWW USING MATLAB'S INBUILT FITTING
        % % % % % %         % WEIGHTS METHOD ======================================================
        % % % % % %     case 'weights'
        % % % % % %
        % % % % % %         w = varargin{find(strcmpi(varargin,'weights'))+1};
        % % % % % %         w = w(:); % column vector
        % % % % % %
        % % % % % %         % Fit: 'untitled fit 1'.
        % % % % % %         [xData, yData, weights] = prepareCurveData( x, y, w );
        % % % % % %
        % % % % % %         % Set up fittype and options.
        % % % % % %         ft = fittype( 'a*cos(2*pi*x./360) + b*sin(2*pi*x./360) + c', 'independent', 'x', 'dependent', 'y' );
        % % % % % %         opts = fitoptions( 'Method', 'NonLinearLeastSquares' );
        % % % % % %         opts.Display = 'Off';
        % % % % % %         opts.StartPoint = [sqrt(2) sqrt(2) 0];
        % % % % % %         opts.Weights = weights;
        % % % % % %
        % % % % % %         % Fit model to data.
        % % % % % %         [fitresult, ~] = fit( xData, yData, ft, opts );
        % % % % % %
        % % % % % %         F = nan(3,1);
        % % % % % %         F(1) = fitresult.a;
        % % % % % %         F(2) = fitresult.b;
        % % % % % %         F(3) = fitresult.c;
        
        % WEIGHTS METHOD ======================================================
    case 'weights'
        
        w = varargin{find(strcmpi(varargin,'weights'))+1};
        w = w(:); % column vector
        
        % remove any duds found earlier
        w = w(~nanlocs);
        
        % Very slightly faster in singles:
        x = single(x); y = single(y); w = single(w);
        
        % Use the sneaky outer product again :) tiny bit faster than
        % repmat I think.
        % W = ones(size(w)) * w' .* eye(length(w));
        % ^ still rather slow :(
        % EDIT: I've removed the above but want to leave it because of its
        % simplicity. It's faster not to generate an NxN matrix of
        % weightings computationally, but when writing this up analytically I'd use it.
        
        % NEW: allow for simultaneous multi-component fit, like so:
        % A*cos(2*pi*t/periods(1)) + B*sin(2*pi*t/periods(1)) ...
        % C*cos(2*pi*t/periods(2)) + D*sin(2*pi*t/periods(2)) ...
        % ... + E;
        
        % why not use the outer product to assemble the sines and
        % cosines... (nice!)
        cosparts = cos((2*pi*xfinite) * (1./period_estimate)');
        sinparts = sin((2*pi*xfinite) * (1./period_estimate)');
        oneparts = ones(size(xfinite));
        
        % assemble them:
        X = ones(length(y),(2*length(period_estimate))+1);
        sincos = 1;
        for i = 1:2:(2*length(period_estimate))
            X(:,i)   = cosparts(:,sincos);
            X(:,i+1) = sinparts(:,sincos);
            sincos = sincos+1;
        end
        X(:,end) = oneparts;
        
        % Matrix Inversion!
        % F = inv(X' * W * X) * X' * W * y;
        % F =    (X' * W * X) \ (X' * W * y); % <- this one is about 50% faster
        F =    (X' * (repmat(w,1,size(X,2)) .* X)) \ (X' * (w .* y)); % <- slightly faster still!
%         F =    (X * (repmat(w,1,size(X,2)) .* X)') \ (X' * (w' .* y')); % <- slightly faster still!
        
        % Return to double:
        F = double(F);
        
end
% =========================================================================

% Reshape F ino something for manageable for multiple-period fits:
final_form = zeros(length(period_estimate),3);

Fr = reshape(F(1:(2*length(period_estimate))),[2 length(period_estimate)])';
final_form(1:length(period_estimate),1:2) = Fr;
final_form(length(period_estimate),3) = F(end);

F = final_form;

% Evaluate the fitted sinusoid: (use the original x-range, NaNs included)
yfit = zeros(size(x));
for i = 1:length(period_estimate)
    yfit = yfit + ...
        F(i,1).*cos((2.*pi.*x) ./ period_estimate(i)) + ...
        F(i,2).*sin((2.*pi.*x) ./ period_estimate(i)) + ...
        F(i,3);
end

% p = 0;
% yfit = F(end);
% for i = 1:2:(size(X,2)-1)
%     p = p + 1;
%     yfit = ...
%         yfit + ...
%         F(i)*cos((2*pi*x) ./ period_estimate(p)) + ...
%         F(i+1)*sin((2*pi*x) ./ period_estimate(p));
% end

% yfit = F(3) + F(1)*cos((2*pi*x) ./ period_estimate) + F(2)*sin((2*pi*x) ./ period_estimate);

yfit = reshape(yfit,origsize);

% if everything was NAN, output all NANs.
if all(isnan(IN)) || all(isnan(x)) || isempty(IN) || isempty(x)
    yfit = nan(size(yfit));
    F = nan(size(F));
end

% output(s): might add more in future
switch nargout
    case {1,0}
        varargout{1} = yfit;
    case 2
        varargout{1} = yfit;
        varargout{2} = F;
end

%% PLOTTING
if nargin > 3
    if any(strcmpi(varargin,'plot'))
        figure; hold all; grid on;
        hold on; plot(x,IN,'.k');
        hold on; plot(x,yfit,'.b');
    end
end
end







%
%
%
% return
% %
%
% %% least squares sinusoid fitting:
%
% p_est = 25;
%
%  t = (1:100)';
%  X = ones(100,3);
%  X(:,2) = cos((2*pi)/p_est*t);
%  X(:,3) = sin((2*pi)/p_est*t);
%  y = 2*cos((2*pi)/p_est*t-pi/4)+randn(size(t));
%  y = y(:);
%  beta = X\y;
%  yhat = beta(1)+beta(2)*cos((2*pi)/p_est*t)+beta(3)*sin((2*pi)/p_est*t);
%
%  figure;
%  plot(t,y,'b');
%  hold on
%  plot(t,yhat,'r','linewidth',2);
%
%
%  %% try sine fitting wind:
%
% u = reshape(squeeze(HWD.Data.u(5,:,100:110)),[1 24*11]);
%
% IN = u;
%
% IN = IN(:); % linearise to column vector
% t = (1:length(IN))'; % must all be column vectors
% % but does it actuslly have to be regularly spaced?! Investigate!
%
% % use only non-nans: (crucial step!)
% finiteinds = ~isnan(IN);
% t = t(finiteinds);
% IN = IN(finiteinds);
%
% p_est = 24; % period estimate, elements, must be regularly spaced
%
%  X = ones(length(IN),3);
%  X(:,2) = cos((2*pi)/p_est*t);
%  X(:,3) = sin((2*pi)/p_est*t);
%
% %  y = 2*cos((2*pi)/p_est*t-pi/4)+randn(size(t));
% %  y = y(:);
%
% B = X \ IN;
%
% F = B(1)+B(2)*cos((2*pi)/p_est*t)+B(3)*sin((2*pi)/p_est*t);
%
% figure; plot(t,IN);
% hold on; plot(t,F,'r','linewi',2)
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
