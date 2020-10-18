function T = theta2t(theta,p)


RoverCp=0.2896;
p0=1000.;


% theta=t*(p0/p)^RoverCp

T = theta .* (p./p0).^RoverCp;


% T = theta .* (Prs./1000).^0.2896;

end

