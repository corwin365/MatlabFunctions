

% IN = Airs.Tp;
%
% nfreqs = 1000;
%
% c = [0.25 0.25 0.25];
% point_spacing = [xt_spacing at_spacing z_spacing];
% minwavelengths = [25 25 6];
% maxwavelengths = [10000 10000 45];


% We expect IN to be in the form XxYxZ, where the phase shift method will be
% applied to the LAST dimension, i.e. the dimension 3 here.


function OUT = nph_2dst_plus1(IN,nfreqs,point_spacing,c,varargin)

% first, get sizes
IN = squeeze(IN);
sz = size(IN);

% ERROR CHECKING:
if length(sz) ~= 3
    error('Expected 3D input data.')
end
if ~isscalar(nfreqs)
    error('Please enter a scalar value for number of 3D frequencies to be analysed.')
end
if length(point_spacing) ~= 3 || length(c) ~= 3
    error('Expected "point_spacing" and "c" to be 1x3 vectors [X Y Z].')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STEP 1: Find dominant freqs using 3d FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the 3D fft to find the dominant [NFREQS] freqs in the data.
% Then only use the horizontal freqs.
% disp('finding dominant horz frequencies...')

% distance between layers:
dz = point_spacing(3);

% find dominant freqs in 3D:
[scales3d,scalemag] = nph_finddomfreqs(IN,nfreqs);

% flip negative Z scales. Because of the ambiguity, anywhere where the sign
% of the vertical freq is negative, use the flipped horizontal freqs, ie
% 5,5,-3 => -5,-5,3. This means that we assume a positive vertical
% wavenumber in all cases before we even analyse with the 2DST. This can
% easily be flipped afterwards depending on what assumption we want to
% make.
neglocs = scales3d(3,:) < 0;
scales3d(:,neglocs) = -scales3d(:,neglocs);

% % combine with the opposing freq directions, to ensure that we get all +ve
% % and -ve phase shifts:
% scales3d = [scales3d -scales3d];
% scalemag = [scalemag  scalemag];

% Find unique combinations of scales:
[~,ia,~] = unique(scales3d(1:2,:)','rows','stable');
scales3d = scales3d(:,ia);
scalemag = scalemag(ia);
% the 'stable' is important to keep them sorted by magnitude

% sort them by magnitude again just to make sure:
[scalemag,ib] = sort(scalemag,'descend');
scales3d = scales3d(:,ib);


scales = scales3d(1:2,:);
nscales = length(scales);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STEP 2: Analyse for those freqs with the 2DST for each level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('computing 2DSTs...')

allST = struct;

for z = 1:size(IN,3)
    allST(z).ST = nph_ndst_dev(IN(:,:,z),scales,point_spacing(1:2),c(1:2),varargin{:},'full');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STEP 3: Compute cospectra and phase shifts between levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this latest approach, we compute the cospectra between all adjacent
% layers, then average these cospectra over a given vertical range
% disp('computing phase shifts...')

% First, need to specify an index matrix so that we can do some logical indexing
% later to find freqs/phases/amplitudes corresponding to the dominant
% covarying amplitudes.
indmap = permute(repmat((1:nscales)',1,size(IN,1),size(IN,2)),[2 3 1]);

% Also need to do one so we can find the co-varying horizontal wavelengths
% that these phase shifts corresponded to:
lhmap1 = permute(repmat((scales(1,:))',1,size(IN,1),size(IN,2)),[2 3 1]);
lhmap2 = permute(repmat((scales(2,:))',1,size(IN,1),size(IN,2)),[2 3 1]);
% and convert these to wavelengths:
fullwid1 = size(IN,1)*point_spacing(1);
fullwid2 = size(IN,2)*point_spacing(2);
lhmap1 = fullwid1 ./ lhmap1;
lhmap2 = fullwid2 ./ lhmap2;

% preallocate OUTPUT wavelengths:
F = struct; % sorry need this for assigning at the end
F.C = zeros(size(IN));
F.F1 = nan(size(IN));
F.F2 = nan(size(IN));
F.F3 = nan(size(IN));



% Here, we use the scaling parameter c to adjust the vertical range over
% which we average the horizontal cospectra.
% For 3D AIRS, c = 1 gives a vertical weighting FWHM of 3km, i.e. just
% one level and the next. Setting c = 0.25 gives a weighting FWHM of 12km.
zfwhm = dz ./ c(3);
% convert to STD
zstd = zfwhm ./ 2.355;
% zstd = 9; % km

% how many levels either side do we want to consider for Gaussian weighting?
% We can adjust this, it isn't the weighting but number of levels to be
% considered. I found 5 height levels either side was ok, but if you want
% you can just do 2 standard deviations:
levelseitherside = round(1*zstd); % or 1 std? keep it simple?
rng = (levelseitherside*2)+1;

%compute weightings here for the weighted average over these heights
% % % % wvec = exp(- ((1:levelseitherside)*dz ./ (sqrt(2)*zstd) ).^2);
% % % % wvec = [wvec(end:-1:1) wvec];
% since we now do phase shift over 2*dz levels, we don't need an even
% lengthed weighting vector:
wvec = exp(- ((-levelseitherside:levelseitherside)*dz ./ (sqrt(2)*zstd) ).^2);

% esnure the weightings all add up to 1:
wvec = wvec ./ (sum(wvec(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outside of the main loop, go through and compute cospectra for all
% the levels and the next. This is faster - as long as you can keep
% them all in memory.
% EDIT: Ok, I realised that the minimum vertical wavelength you can have
% with this method is 2 x dz. This means it's pointless to do the phase
% shift for each dz - all you're doing is forcing us to rely on tiny phase
% shifts which may be below the error.
% SOLUTION: For each height, do the cospectra for the heights above and
% below. It may seem counter intuitive not to use the specified height, but
% I think it makes sense practically. The given height will just be used
% in the next one.
allC = struct;
for i = 1:size(IN,3)
    switch i
        case 1 % bottom height
            allC(i).C = allST(i).ST.ST .* conj(allST(i+2).ST.ST);
        case size(IN,3) % top height
            allC(i).C = allST(i-2).ST.ST .* conj(allST(i).ST.ST);
        otherwise % all the middle ones
            allC(i).C = allST(i-1).ST.ST .* conj(allST(i+1).ST.ST);
    end
end
% % allC = struct;
% % for i = 1:size(IN,3)-1
% % %     if i == (size(IN,3)-1)
% %         allC(i).C = allST(i).ST.ST .* conj(allST(i+1).ST.ST);
% % %     else
% % %         allC(i).C = allST(i).ST.ST .* conj(allST(i+2).ST.ST);
% % %     end
% % %     % you can even give it a little 3x3 smooth:
% % %     allC(i).C = smoothn(allC(i).C,[3 3 1]);
% % end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ph = nan(size(IN));

% Now compute the weighted average phase shifts for each vertical layer.
% for each height
for z = 1:size(IN,3)
    
    % first, find the range of levels we're going to be phase shifting
    % across:
    zrng = z + (-levelseitherside:levelseitherside);
    
    % fix issues at the top and bottom:
    if min(zrng) < 1
        zrng = 1:max(zrng);
    end
    if max(zrng) > size(IN,3)
        zrng = min(zrng):size(IN,3);
    end
    
    % Take the average (this is reasonable while they are split into
    % real and imaginary components - think of a vector on the argand diagram,
    % we find the average to determine the average phase shift)
    Cavg = zeros(size(allST(1).ST.ST));
    for i = 1:length(zrng)
        
        % % % %         % compute cospectrum between each level and the next one
        % % % %         Chere = allST(zrng(i)).ST.ST .* conj(allST(zrng(i+1)).ST.ST);
        % % % %         % you can even give it a little 3x3 smooth:
        % % % %         Chere = smoothn(Chere,[3 3 1]);
        
        % or just load it from the pre-computed cospectra (faster)
        Chere = allC(zrng(i)).C;
        
        % combine using weighted average into our group for this selected
        % level z:
        Chere = (wvec(i) .* Chere);
        Cavg = Cavg + Chere;
        %         Cavg = Cavg + (wvec(i) .* Chere);
        % because the weights all add up to 1, this should work
    end
    
    % compute phase shift of this average cospectrum
    phasediff = angle(Cavg);
    
    % % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % %     %%%% SELECT ONLY POSITIVE PHASE SHIFTS
    % % % % %     % So earlier we doubled our input scales so that opposing horizontal
    % % % % %     % directions could also be considered, i.e. 5,5 and -5,-5. This means
    % % % % %     % that somewhere in our analysed STs there should be the same wave but
    % % % % %     % with +ve and -ve phase shifts. We can now choose whichever we want,
    % % % % %     % which will flip the horizontal direction as intended.
    % % % % %     Cavg(phasediff <= 0) = 0;
    % % % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
% % % %     %%%% Currently, the phase shift goes from -pi to +pi. But this will only
% % % %     %%%% give half of a full cycle in either direction, so you can never have
% % % %     %%%% more than a half wavelength.
% % % %     %%%% If you want to get a full wavelength, we need it to go from 0 to 2pi.
% % % %     %%%% Think about it. If you've got 3km spacing, and you take phase shift
% % % %     %%%% over 2 layers, so 6km vertical gap. The smallest wavelength you can
% % % %     %%%% have is 6km. But to get this, you'd need a phase shift ~2pi. By
% % % %     %%%% flipping round to -pi to +pi, you limit the minimum wavelength to
% % % %     %%%% twice this.
% % % %     %%%% I haven't worked out how to preserve the sign yet, but maybe you
% % % %     %%%% can't. Perhaps this is related to the ambiguity thing?
% % % %     neglocs = phasediff < 0;
% % % %     phasediff(neglocs) = (2*pi) + phasediff(neglocs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% APPLY MIN/MAX WAVELENGTH CONDITIONS
    % anywhere where the phase shift would give a LZ < MIN or LZ > MAX, set
    % these coefficients to zero so they don't get picked up?
    if any(strcmpi(varargin,'minwavelengths'))
        minlz = varargin{find(strcmpi(varargin,'minwavelengths'))+1}(3);
    end
    if any(strcmpi(varargin,'maxwavelengths'))
        maxlz = varargin{find(strcmpi(varargin,'maxwavelengths'))+1}(3);
    end
    
    % apply:
    kz = phasediff ./ (2 * dz);
    lz = (2*pi) ./ kz;
    
    goodlocs = abs(lz) > minlz | abs(lz) < maxlz;
    Cavg(~goodlocs) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % now find the locations of the maximum covarying amplitudes:
    [~,inds] = max(abs(Cavg),[],3);
    
    % use our index map from earlier to turn these into real indeces to
    % pluck out the values we want:
    locs = indmap == repmat(inds,1,1,nscales);
    
    % get the phase differences:
    phasediff(~locs) = 0;
    phasediff = sum(phasediff,3);
    
    % deal with artefacts from the min/max/conditions
    phasediff(phasediff == 0) = NaN;
    
    Ph(:,:,z) = phasediff;
    
    % compute LZ:
    kz = phasediff ./ (2 * dz);
    lz = (2*pi) ./ kz;
    
    % Finally, if there were literally no freqs that could be selected that
    % were within the min/max freqs, set them to the Nyquist min/max possible:
    if any(strcmpi(varargin,'minwavelengths'))
        lz(abs(lz) < minlz) = 2 * dz;
    end
    if any(strcmpi(varargin,'maxwavelengths'))
        lz(abs(lz) > maxlz) = dz * sz(3);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% SELECT ONLY POSITIVE PHASE SHIFTS / WAVELENGTHS
    % So, remember right at the start we found all negative vertical scales
    % in the dominant scales3d, and flipped them and their corresponding
    % horizontal scales so that they were positive? Well, surely this now
    % means that all phase shifts must be positive with increasing height?
    lz = abs(lz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % and subscribe!
    F.F3(:,:,z) = 1./lz;
    
    % also find the corresponding LH for these locations:
    lh1 = lhmap1; lh1(~locs) = 0; lh1 = sum(lh1,3);
    lh2 = lhmap2; lh2(~locs) = 0; lh2 = sum(lh2,3);
    F.F1(:,:,z) = 1./lh1;
    F.F2(:,:,z) = 1./lh2;
    
    
end

% assemble some other stuff:
freqs = repmat(sz(1:2)',1,size(scales,2)) .* repmat(point_spacing(1:2)',1,size(scales,2));
freqs = freqs ./ scales;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STEP 4: Finally, assemble your outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ST = struct;

ST.type                     = '2DST+1';
ST.IN                       = IN;
ST.scales                   = scales;
ST.point_spacing            = point_spacing;
ST.c                        = c;
ST.freqs                    = freqs;
ST.VerticalAveragingFWHM    = round(zfwhm,1);
% ST.AmplitudeBoosting        = allST(1).ST.AmplitudeBoosting;
% ST.GuidedFourierMode        = allST(1).ST.GuidedFourierMode;
ST.Options                  = allST(1).ST.Options;
ST.Ph                       = Ph;

% select some data properties:
props = {'C','F1','F2','F3','A','R','HA','HR','BoostFactor'};
for p = 1:length(props)
    ST.(props{p}) = nan(size(IN));
end

% now build up each level by level:
for p = 1:length(props)
    
    for i = 1:size(IN,3)-1
        
        switch props{p}
            case {'F1','F2','F3'}
                ST.(props{p})(:,:,i) = F.(props{p})(:,:,i);
                %         case 'C'
                %             ST.(props{p})(:,:,i) = F.C(:,:,i);
                %         case 'A'
                %             ST.(props{p})(:,:,i) = abs(F.C(:,:,i));
                %         case 'R'
                %             ST.(props{p})(:,:,i) = real(F.C(:,:,i));
            otherwise
                ST.(props{p})(:,:,i) = allST(i).ST.(props{p});
        end
        
    end
end

% fix the top layer:
for p = 1:length(props)
    ST.(props{p})(:,:,end) = ST.(props{p})(:,:,end-1);
end

% % % % % % Finally, I think I've got the definition of the vertical wavelength
% % % % % % flipped, so let's try:
% % % % % ST.F3 = -ST.F3;


OUT = ST;




%%%%

%
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% PLOTTY MCPLOTFACE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure; hold all; whitefig; figpos([1 0.6])
%
% %-------------------------------------------------------
% vert_gap    = 0.1;    horz_gap    = 0.075;
% lower_marg  = 0.1;    upper_marg  = 0.1;
% left_marg   = 0.05;   right_marg   = 0.05;
%
% rows = 2; cols = 4;
%
% subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
%
% %--------------------------------------------------------
%
% fs = 16;
%
% %--------------------------------------------------------
%
% switch granule
%     case 1
%         xsec = 51;
%         [~,zsec] = min(abs(selected_alt - 36));
%     case 2
%         xsec = 72;
%         [~,zsec] = min(abs(selected_alt - 36));
%     case 3
%         xsec = 47;
%         [~,zsec] = min(abs(selected_alt - 36));
% end
%
% % which LH do you want?
% LH = 1./quadadd(1./LH1,1./LH2);
% % or?
% % LH = nan(size(IN));
% % for i = 1:size(IN,3)
% %     LH(:,:,i) = 1./quadadd(allST(i).ST.F1,allST(i).ST.F2);
% % end
%
% % make some grids:
% at_vec = at_spacing .* (0:size(IN,2)-1);
% xt_vec = xt_spacing .* (0:size(IN,1)-1);
% z_vec = selected_alt(1:size(IN,3));
% [Y,X] = meshgrid(at_vec,xt_vec);
% [YY,Z] = meshgrid(at_vec,z_vec);
%
% TBP = {...
%     IN(:,:,zsec),abs(LH(:,:,zsec)),abs(LZ(:,:,zsec)),abs(1./ST3.F3(:,:,zsec)),...
%     sq(IN(xsec,:,:))',sq(abs(LH(xsec,:,:)))',sq(abs(LZ(xsec,:,:)))',sq(abs(1./ST3.F3(xsec,:,:)))'};
%
% GX = {Y,Y,Y,Y,YY,YY,YY,YY};
% GY = {X,X,X,X,Z,Z,Z,Z};
%
% nclev = 21;
%
% clims = {pm(5),[0 600],[0 75],[0 75],pm(5),[0 600],[0 75],[0 75]};
%
% clevs = {-20:0.25:20,0:50:1000,0:1:100,0:1:100,-20:0.25:20,0:50:1000,0:1:100,0:1:100};
%
% cmaps = {...
%     cbrew('RdBu',nclev),...
%     flipud(cbrew('Greens',nclev)),...
%     flipud(cbrew('Blues',nclev)),...
%     flipud(cbrew('Blues',nclev)),...
%     cbrew('RdBu',nclev),...
%     flipud(cbrew('Greens',nclev)),...
%     flipud(cbrew('Blues',nclev)),...
%     flipud(cbrew('Blues',nclev))};
%
%
% for ax = 1:8
%
%     subplot(rows,cols,ax);
%
%     axx = gca;
%
%     TBP{ax} = nph_fix2clims(TBP{ax},clims{ax});
%
%     hold on; contourf(GX{ax},GY{ax},TBP{ax},clevs{ax},'edgecolor','none'); grid on;
%
%     cbar = nph_colorbar;
%     colormap(gca,cmaps{ax})
%     clim(clims{ax})
%
%     xlim(minmax(at_vec));
%     xtick(0:500:2500)
%
%     switch ax
%         case {1,2,3,4}
%             hold on; plot(axx.XLim,xt_vec([xsec xsec]),'linest','--','linewi',1.5,'color','m');
%             ylabel('XT Distance (km)')
%         case {5,6,7,8}
%             hold on; plot(axx.XLim,z_vec([zsec zsec]),'linest','--','linewi',1.5,'color','m');
%             ylabel('Altitude (km)')
%             ylim([0 60])
%     end
%     xlabel('AT Distance (km)')
%
%     switch ax
%         case {3,4,7,8}
%             hold on; [C,h] = contour(GX{ax},GY{ax},TBP{ax},15:5:30,'edgecolor',rgb(.15));
%             hold on; clabel(C,h,'fontsize',0.6*fs,'color',rgb(.15));
%     end
%
%
%     % TITLES
%     switch ax
%         case 1
%             tit = strsplit(g1name,'_');
%             title(['AIRS ' tit{2} ' d' tit{3} ' g' tit{4}(1:3) '-' g2name(15:17)])
%         case 2
%             title('LH_{2D+1}')
%         case 3
%             title('LZ_{2D+1}')
%         case 4
%             title('LZ_{3DST}')
%     end
%
%     % letters
%     nph_text([0 1.025],['(' alphabet(ax) ')'],'fontsize',1.5*fs);
%
%     setfont(fs)
%     set(gca,'linewi',1.5,'layer','top')
%
%
% end
%
% return
%
%
%
% %% EXPORT? ================================================================
%
%
% savename = ['~/Desktop/testing_2DSTplus1_' num2str(granule)];
% disp(['Exporting to ' savename '...'])
%
% nph_saveas(gcf,savename,'png')
%
% disp('Done.')
%
%
%
%
%
%
%
% return
%
%
%
%
%
%










