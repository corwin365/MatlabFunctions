function u = mean_angle(phi,dimension)
	u = angle(nanmean(exp(i*pi*phi/180),dimension))*180/pi;
end
