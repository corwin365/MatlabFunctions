function [c] = quadadd(varargin)
%Pythagorian addition of terms

Sigma = 0;
for iArg=1:1:numel(varargin)
  Sigma = Sigma + varargin{iArg}.^2;
end
c = sqrt(Sigma);



return; end

