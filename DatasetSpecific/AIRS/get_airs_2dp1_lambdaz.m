%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D+1 st, version 1. This computes independent horizontal wavenumber fits.
%called from gwanalyse_airs_3d, probably not useful on its own
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NewFields = get_2dp1_lambdaz(Airs,Extra,Settings)


    
  %glue the extra levels we retained to the top and bottom
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Vars = {'Tp','ret_z','BG','ret_temp'};
  Dims = [3,1,3,3];
  for iVar=1:1:numel(Vars)
    Working.(Vars{iVar}) = Airs.(Vars{iVar});
    if isfield(Extra,'Bottom');
      if numel(Extra.Bottom.ret_z) > 0;
        Working.(Vars{iVar}) = cat(Dims(iVar),Extra.Bottom.(Vars{iVar}),Working.(Vars{iVar}));
      end;
    end
    if isfield(Extra,'Top');
      if numel(Extra.Top.ret_z) > 0;
        Working.(Vars{iVar}) = cat(Dims(iVar),Working.(Vars{iVar}),Extra.Top.(Vars{iVar}));
      end;
    end;
    Airs.(Vars{iVar}) = Working.(Vars{iVar});
  end; clear iVar Working Vars Dims

  
  %% 2DST every level, and store as a height-weighted sum
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %define scales
  Scales{1} = 1:1:size(Airs.l1_lat,1)/2;
  Scales{2} = 1:1:size(Airs.l1_lat,2)/2; Scales{2} = [reverse(-Scales{2}),Scales{2}];
  
  %thin scales?
  if Settings.Thin
    
    %some coarser scales that should still resolve wave structure well
    %these numbers assume a standard-sized AIRS granule
    lx = unique(floor((20.* 90)./logspace(log10(20),log10( 20.*90),30))); %UP TO 30 freqs xt
    ly = unique(floor((18.*135)./logspace(log10(18),log10(18.*135),30))); %UP TO 30 freqs at    
    
    Scales{1} = lx; Scales{2} = ly;
    clear lx ly
    
% % %     s1 = Scales{1}; Scales{1} = s1([1:1:5,6:2:end]);
% % %     s2 = Scales{2}; Scales{2} = s2([1:1:5,6:2:end]);
% % %     clear s1 s2
  end

  textprogressbar('Initial ST pass ')
  for iLevel=1:1:numel(Airs.ret_z);
    
    %st level
    ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel), ...
                           'FullST',true,           ...
                           'c',Settings.c1,         ...
                           'Scales',Scales);
    
    %compute weight
    if Settings.Weight == 1; CFac = exp(-(Airs.ret_z(iLevel)-42) ./ (2*scale_height(42))); %42 is arbitrary, but must be consistent across levels
    else                     CFac = 1;
    end
    ST.ST  = ST.ST .* CFac;
    
    %store
    if iLevel == 1;
      %first level: just store
      Store.ST   = ST.ST;
      freqs = ST.freqs; %needed later
    else
      %subsequent level: integrate into existing average
      Store.ST = ((iLevel-1).*Store.ST + ST.ST)./iLevel;
    end
    
    %and loop
    textprogressbar(iLevel./numel(Airs.ret_z).*100)
  end
  textprogressbar('!')
  clear iLevel CFac ST
 
  
  %% find spectral maxima in the two fields
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %find spectral peaks
  PeakStore = NaN([size(Airs.l1_lat),Settings.NPeaks,3]);
  
  textprogressbar('Finding peaks ')
  for iX=1:1:size(Airs.l1_lat,1);
    for iY=1:1:size(Airs.l1_lat,2);
      
      %pull out x/y location spectral space
      Frame = abs(Store.ST(  :,:,iX,iY));
      
      %peak finder setup
      if Settings.NPeaks > 1; cent = FastPeakFind(Frame,0,Settings.Filt); %find n peaks
      else                    cent = [];                                  %used in logic below
      end 
      
      if numel(cent) == 0;
        %if we either found no peaks or only wanted a single peak,
        %take the overall-maximum as the largest peak
        [~,idx] = max(Frame(:));
        [x,y] = ind2sub(size(Frame),idx);
        cent(2) = x; cent(1) = y;
        clear idx x y
      end
      

      %convert from list to [x1,y1;x2,y2;...];
      Cent(:,1) = cent(1:2:end); %position on k_(along-track)  axis
      Cent(:,2) = cent(2:2:end); %position on k_(across-track) axis
      
      %find values of these indices
      for iZ=1:1:size(Cent,1); Cent(iZ,3) = Frame(Cent(iZ,2),Cent(iZ,1)); end
      
      %remove small values
      Good = find(Cent(:,3) > Settings.Threshold);
      Cent = Cent(Good,:);
      clear Good
      
      %sort the data, select the number of peaks desired, and store
      [~,idx] = sort(Cent(:,3),'desc');
      if size(Cent,1) > Settings.NPeaks; Cent = Cent(idx(1:Settings.NPeaks),:);
      else;                              Cent = Cent(idx,:);
      end
      clear idx
      PeakStore(iX,iY,1:size(Cent,1),:) = Cent;
      clear Cent
      
    end
    textprogressbar(iX./size(Airs.l1_lat,1).*100)
  end
  clear iX iY iZ Frame iMode Working cent
  textprogressbar('!')
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% we now have the N strongest spectral signals at each geographic point,
  % found using a coarse c to provide reasonable spectral fit.
  %now put back and do a finer fit to these frequencies
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %find all the unique scales of importance in the data
  clear Scales
  s1 = unique(PeakStore(:,:,:,1,:)); s1(isnan(s1)) = []; Scales{1} = s1; clear s1;
  s2 = unique(PeakStore(:,:,:,2,:)); s2(isnan(s2)) = []; Scales{2} = s2; clear s2;
  
  %frequencies from original ST, to fix horizontal wavelengths
  f1 = freqs{1};
  f2 = freqs{2};
  
  %create storage arrays for the fine STs
  Store.STs  = complex(NaN([numel(Scales{1}),numel(Scales{2}), ...
                       size(Airs.l1_lat),numel(Airs.ret_z)], ...
                       'single'));
  
  textprogressbar('Computing fine STs ')
  for iLevel=1:1:numel(Airs.ret_z);
    
    %do ST
    ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel), ...
                           'FullST',true,          ...
                           'c',Settings.c2(1:2),   ...
                           'Scales',Scales);
    %store ST
    Store.STs(:,:,:,:,iLevel) = ST.ST;
    
    textprogressbar(iLevel ./ size(Airs.ret_z,1) .* 100);
    
  end
  clear ST iLevel
  textprogressbar('!')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% find the requested cospectra, and store the resulting GW properties
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %results arrays
  Store.A  = NaN([size(Airs.l1_lat),size(Airs.ret_z,1),numel(Settings.Steps),Settings.NPeaks]);
  Store.L  = Store.A;
  Store.F1 = Store.A;
  Store.F2 = Store.A;
  
  
  textprogressbar('Computing cospectra ')
  for iLevel=1:1:size(Airs.ret_z,1);
    for iStep = 1:1:numel(Settings.Steps)
      
      %identify the level-pair to compute a cospectrum over
      ST1 = Store.STs(:,:,:,:,iLevel); %this level is one of them...
      if iLevel+Settings.Steps(iStep) <= size(Airs.ret_z,1);      
        ST2 = Store.STs(:,:,:,:,iLevel+Settings.Steps(iStep));
      else; continue; end
      dZ = Airs.ret_z(iLevel) - Airs.ret_z(iLevel+Settings.Steps(iStep));
      
      %take the cospectrum.
      CoSpectrum = ST1 .* conj(ST2);
      
      %then pull out the indices we need
      for iX=1:1:size(Airs.l1_lat,1);
        for iY=1:1:size(Airs.l1_lat,2);
          for iZ=1:1:Settings.NPeaks;
              if ~isnan(PeakStore(iX,iY,iZ,3));
                
                %find new indices
                idx1 = closest(PeakStore(iX,iY,iZ,1),Scales{1});
                idx2 = closest(PeakStore(iX,iY,iZ,2),Scales{2});
                
                %cospectral value
                CS = CoSpectrum(idx1,idx2,iX,iY);
                
                %amplitude
                Store.A(iX,iY,iLevel,iStep,iZ) = sqrt(abs(CS));
                
                %vertical wavelength
                dPhi = angle(CS)./(2*pi);
                
                Store.L(iX,iY,iLevel,iStep,iZ) = dZ./dPhi;
                
                %horizontal wavelengths
                Store.F1(iX,iY,iLevel,iStep,iZ) = f1(PeakStore(iX,iY,iZ,2));
                Store.F2(iX,iY,iLevel,iStep,iZ) = f2(PeakStore(iX,iY,iZ,1));
                
            end
          end
        end
      end
    end
    textprogressbar(iLevel ./ size(Airs.ret_z,1) .* 100);
  end;
  clear iLevel iX iY iZ iStep CoSpectrum CS dPhi
  clear idx1 idx2 dZ ST1 ST2 f1 f2 freqs PeakStore
  Store = rmfield(Store,{'STs','ST'});
  textprogressbar('!')


  %combine and return
  NewFields.A  =    squeeze(Store.A( :,:,:,1,1));
  NewFields.m  = 1./squeeze(Store.L( :,:,:,1,1));
  NewFields.F1 =    squeeze(Store.F1(:,:,:,1,1));
  NewFields.F2 =    squeeze(Store.F2(:,:,:,1,1));

return
