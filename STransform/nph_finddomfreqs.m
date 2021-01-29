

% This function returns the N largest spectral amplitude scale combinations
% of the FFT of the input data.

% It tells you what the dominant freqs (or at least, their positions in
% fourier space) of the input data were.

% Supports 1D,2D,3D,4D inputs.

% Inputs = IN (input data), nfreqs (number of largest scales)
% Outputs: scales (a list of scale combinations corresponding to dominant freqs)

% Usage: scales = nph_finddomfreqs(IN,nfreqs)


%%%% ====================
% alright here's where we use the FFT to find the top XXXX frequencies
% present in the input data, then work out what their scales would be.
% Hold my beer...

% IN = Airs.Tp(:,1:134,16);
% nfreqs = 10;

function varargout = nph_finddomfreqs(IN,nfreqs)

% select N-D type:
sz = size(IN);
type = sum(sz ~= 1);

% take FFT:
F = fftn(IN);

ab = abs(F(:));
im = imag(F(:));

% sort by absolute spectral power
[ab,ib] = sort(ab,'descend');

% also rearrange the imaginary comps by this sorting:
im = im(ib);

% exclude the DC components with imag parts == 0
ab = ab(im ~= 0);
ib = ib(im ~= 0);

% now reshape:
% this should always work - after the zeros are taken out there should
% always be an even number of complex conjugate pairs remaining.
% note: you need the transpose ' here due to the way reshape re-lists things.
abr = reshape(ab',2,length(ab)/2);
ibr = reshape(ib',2,length(ib)/2);

% Make a coord system in fft space
v = struct;
for n = 1:type
    switch mod(sz(n),2)
        case 0 % if it's even
            N = (sz(n)/2)-1;
            v(n).vec = ifftshift([0 -N:N]);
        case 1 % if it's odd
            N = (sz(n)-1)/2;
            v(n).vec = ifftshift(-N:N);
    end
end

% what were the scales that related to these locations?
ii = struct;
switch type
    case 1
        ii(1).i = ind2sub(size(F),ibr(1,:));
    case 2
        [ii(1).i,ii(2).i] = ind2sub(size(F),ibr(1,:));
    case 3
        [ii(1).i,ii(2).i,ii(3).i] = ind2sub(size(F),ibr(1,:));
    case 4
        [ii(1).i,ii(2).i,ii(3).i,ii(4).i] = ind2sub(size(F),ibr(1,:));
end

% RESET SCALES:
scales      = cell(1,type);
scalemag    = cell(1,type);
% wavelengths = cell(1,type);
% physical_dims = nan(1,type);
goodinds = ones(1,length(ii(1).i));

% LIMIT TO MIN/MAX SCALES, NON-ZERO SCALES:
for n = 1:type
    
    scales{n} = v(n).vec(ii(n).i);
    scalemag{n} = abr(1,:);
    
    % remove zero scales:
    goodinds = all(cat(1,goodinds,scales{n} ~= 0));
    
%     % covert to wavelengths if needed:
%     physical_dims(n) = point_spacing(n) * sz(n);
%     wavelengths{n} = physical_dims(n) ./ scales{n};
%     
%     % apply MIN wavelength cutoff:
%     if minwavelengthsflag
%         goodinds = all(cat(1,goodinds,abs(wavelengths{n}) >= abs(minwavelengths(n))));
%     end
%     
%     % apply MAX wavelength cutoff:
%     if maxwavelengthsflag
%         goodinds = all(cat(1,goodinds,abs(wavelengths{n}) <= abs(maxwavelengths(n))));
%     end
    
end

% Apply this externally to the above loop so that it's the same for all
% dimensions:
for n = 1:type
    scales{n} = scales{n}(goodinds);
    scalemag{n} = scalemag{n}(goodinds);
%     wavelengths{n} = wavelengths{n}(goodinds);
end

% if the user has asked for more freqs than there are available, limit it:
for n = 1:type
    if nfreqs > length(scales{n})
        warning(['Too many frequencies requested. Limiting to ' num2str(length(scales{n})) '.'])
        nfreqs = length(scales{n});
    end
end

% Finally, select the top NFREQS from our formatted scales:
for n = 1:type
    scales{n} = scales{n}(1:nfreqs);
    scalemag{n} = scalemag{n}(1:nfreqs);
%     wavelengths{n} = wavelengths{n}(1:nfreqs);
end

% % % % % % Try this: for the collapsed spectrum we care about the order
% % % % % % (slightly) that the frequencies are called. So let's sort them in
% % % % % % order of lowest to highest scales. Not sure if it'll make a
% % % % % % difference but lets's see.
% % % % % scalemag = reshape([scales{1:type}],type,nfreqs);
% % % % % scalemag = sqrt(sum(scalemag.^2,1));
% % % % % [~,ord] = sort(scalemag,'ascend');
% % % % % 
% % % % % for n = 1:type
% % % % %     scales{n} = scales{n}(ord);
% % % % % %     wavelengths{n} = wavelengths{n}(ord);
% % % % % end

% Convert to vectors:
scales_vec      = zeros(type,nfreqs);
scalemag_vec    = zeros(type,nfreqs);
% wavelengths_vec = zeros(type,length(ord));
for n = 1:type
    scales_vec(n,:) = scales{n};
    scalemag_vec(n,:) = scalemag{n};
%     wavelengths_vec(n,:) = wavelengths{n};
end
scales = scales_vec;
scalemag = scalemag_vec(1,:);
% wavelengths = wavelengths_vec;

switch nargout
    case {1,0}
        varargout{1} = scales;
    case 2
        varargout{1} = scales;
        varargout{2} = scalemag;
end




% end





