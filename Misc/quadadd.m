function [c] = quadadd(a,b)
%Pythagorian addition of the two terms

% % c = sqrt(a.^2 + b.^2);

%turns out there is a built in function that does this with some sanity checking. use this instead!
c = hypot(a,b);
return
end

