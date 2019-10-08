function [Attenuation] = compute_airs_attenuation(Channel)

%compute AIRS attenuation in a given channel via convolution of the weighting kernel and a theoretical wave
%based on Alexander and Barnet (JAS, 2007)

%load weighting functions
Weights = load([LocalDataDir,'/AIRS/airsweights.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract weighting function for the channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%extract channel 
WF  = Weights.Weights(Channel,:);

%extract z-axis for channel
Z   = Weights.Heights./1000; 
%first three points are small, but must increase monotonically 
%all values added here are basically zero relative to the heights we're interested in
Z(1:3) = 0.001.*(1:1:3); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%produce a regular height grid
HeightGrid = 0:1:80;

%find peak of WF
[~,MaxResponse] = max(WF);
HL = Z(MaxResponse);
[~,Lev] = min(abs(HL-HeightGrid));

%interpolate weighting function onto grid
WF2 = inpaint_nans(interp1(Z,WF,HeightGrid));

%produce scales to loop over
LambdaScale = 2:2:100;
Step = pi./10; PhaseScale  = 0:Step:(2.*pi)-Step; clear Step

%produce storage arrays
WCs = NaN(numel(HeightGrid),numel(LambdaScale),numel(PhaseScale));


%loop over possible wavelengths
for iLambda = 1:1:numel(LambdaScale);
  Lambda = LambdaScale(iLambda);
    %loop over possible phases within that wavelength
    for jPhase = 1:1:numel(PhaseScale);
      Phase = PhaseScale(jPhase);
      
      %generate wave of given wavelength, amplitude unity
      Wave = cos(HeightGrid./Lambda.*(2.*pi)+Phase);
      
      
      %convolve it with the weighting function
      WC = conv(WF2,Wave,'same');

      %store the result
      WCs(:,iLambda,jPhase) = WC;
      
    end; clear jPhase
  
end; clear iLambda

%find optimal phase at each height
WCs = squeeze(nanmax(WCs,[],3));

% %normalise at each height
for iLevel=1:1:numel(HeightGrid);
  WCs(iLevel,:) = WCs(iLevel,:)./nanmax(squeeze(WCs(iLevel,:)));
end

%extract height level of the WF peak
WFLev = WCs(Lev,:);



%remove short wavelengths (less than FWHM)
%function seems to slightly overestimate FWHM to my eye - 0.9 is to correct for this
%this may be due to the skew of the peaks, which extend higher vertically up than vertically down
FWHM = 0.90.*fwhm(HeightGrid,WF2);
WFLev(LambdaScale < FWHM) = NaN;
Attenuation.FWHM = FWHM.*1000;

%prep for return
Attenuation.LambdaScale = LambdaScale.*1000; %m
Attenuation.Attenuation = WFLev;

%also return the weighting function used, for inspection
Attenuation.WF.WF = WF2;
Attenuation.WF.Z  = HeightGrid;
