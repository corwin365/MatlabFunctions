function Pressure = jra55_prs(LnSP,NLevs)

% % 
% % NLevs = 60;
% % LnSP = ln(1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert JRA-55 levels to pressure
%array safe
%Corwin wright, 09/MAR/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 
% % % % 
% % % % 
% % % % %get the A and B coefficients for this number of levels
% % % % switch NLevs
% % % %   case 60;
% % % %     A = [0.,0.,0.,0.,0.,0.,0.,0.,133.05101128,364.90414887,634.60271645,959.79716729,1347.68004166,1790.90739595,2294.84168995,2847.48477771,3468.87148812,4162.95646297,4891.8808325,5671.8242398,6476.71299639,7297.46989472,8122.15979125,8914.08220106,9656.1819105,10329.43617777,10912.63844424,11369.64783084,11695.37159747,11861.25308739,11855.43431635,11663.35536558,11285.40406449,10729.94940557,10014.61505351,9167.24703583,8226.2449077,7201.5689803,6088.67300853,4950.,4000.,3230.,2610.,2105.,1700.,1370.,1105.,893.,720.,581.,469.,377.,301.,237.,182.,136.,97.,65.,39.,20.,0.];
% % % %     B = [0,9.97000000e-01,9.94000000e-01,9.89000000e-01,9.82000000e-01,9.72000000e-01,9.60000000e-01,9.46000000e-01,9.26669490e-01,9.04350959e-01,8.79653973e-01,8.51402028e-01,8.19523200e-01,7.85090926e-01,7.48051583e-01,7.09525152e-01,6.68311285e-01,6.24370435e-01,5.80081192e-01,5.34281758e-01,4.88232870e-01,4.42025301e-01,3.95778402e-01,3.50859178e-01,3.07438181e-01,2.65705638e-01,2.25873616e-01,1.89303522e-01,1.55046284e-01,1.24387469e-01,9.64456568e-02,7.23664463e-02,5.21459594e-02,3.57005059e-02,2.28538495e-02,1.33275296e-02,6.73755092e-03,2.48431020e-03,1.13269915e-04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% % % % end
% % % % 
% % % % %produce array for final pressure data
% % % % InSize = size(SurfacePressure);
% % % % PHalf = NaN([InSize,NLevs+1]);
% % % % 
% % % % %use A and B coefficients to get PHalf for each level
% % % % for iLev=1:1:NLevs+1; PHalf(:,:,iLev) = A(iLev)+B(iLev).*SurfacePressure; end
% % % % 
% % % % %then get pressure for each level
% % % % Pressure = NaN([InSize,NLevs]);
% % % % for iLev=2:1:NLevs;
% % % %   Pressure(:,:,iLev) = 0.5.*(PHalf(:,:,iLev)+PHalf(:,:,iLev+1));
% % % % end
% % % % 
% % % % %scale to hPa
% % % % Pressure = Pressure ./100;


%I only care about the stratosphere for now. So:

if NLevs == 60;
  Pressure =  [99850.0, 99550.0, 99149.9, 98549.8, 97699.6, 96599.4, 95299.1, 93698.6, 91798.2, 89697.8, 87347.0, 84696.1, 81795.4, 78694.6, 75444.0, 72042.9, 68441.7, 64741.2, 60990.1, 57189.5, 53388.7, 49587.9, 45837.6, 42187.2, 38636.8, 35186.3, 31886.6, 28736.1, 25736.4, 22885.7, 20186.0, 17686.4, 15386.9, 13287.5, 11388.1, 9689.0, 8164.3, 6763.8, 5515.0, 4466.6, 3608.1, 2914.5, 2353.0, 1898.9, 1532.0, 1235.1, 997.1, 804.9, 649.3, 524.0, 422.2, 338.3, 268.4, 208.9, 158.4, 116.0, 80.5, 51.5, 29.0, 10.0]./100;
else
  disp('error: no 60 levels');
  stop
end


return