%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ensure data is on a regular Z grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Data = func_regularise_data_z(Data,Settings)


%first, check if the data is ALREADY regular AND ascending in z. Counts if the full distribution is within 5% of the mean.
%if it's fine, we don't need to proceed
dZ_distrib = unique(diff(Data.Alt,1,2));


if min(diff(Data.Alt,1,2),[],'all') > 0 & range(dZ_distrib) < 0.1.* nanmean(dZ_distrib); return; end


%ok, we need to regularise. First, work out a scale
NewZ = nanmin(Data.Alt(:)):nanmean(abs(dZ_distrib)):nanmax(Data.Alt(:));
clear dZ_distrib


%create struct to store new data
NewData = struct();
Fields  = fieldnames(Data); Fields(ismember(Fields,'Alt')) = [];
for iF=1:1:numel(Fields); NewData.(Fields{iF}) = NaN([size(Data.Alt,1),numel(NewZ)]); end

%now interpolate EVERYTHING onto the new scale
Fields  = fieldnames(Data); Fields(ismember(Fields,'Alt')) = [];
if Settings.Verbose == 1; disp('--> Data not on a regular and common height grid, interpolating to [ min Z: mean dZ : max Z ] '); end

for iF=1:1:numel(Fields)

  %get data fields
  Fin  =    Data.(Fields{iF});
  Fout = NewData.(Fields{iF}); 

  %check we're working on a 2D profiles x heights array, and just pass it straight through if not
  %do this by chekcing if it's the same size as the 'Alt' array (arbitrary choice)
  if ~isequal(size(Fin),size(Data.Alt)); NewData.(Fields{iF}) = Fin; continue; end


  %interpolate profile-by-profile
  for iProfile=1:1:size(Data.Alt,1);
    Good = find(~isnan(Data.Alt(iProfile,:)));
    if numel(Good) < 2; continue; end
    Fout(iProfile,:) = interp1(Data.Alt(iProfile,Good),Fin(iProfile,Good),NewZ);
  end

  %store data
  NewData.(Fields{iF}) = Fout;
end

%store new altitudes
NewData.Alt = repmat(NewZ,size(Data.Alt,1),1);

%copy over, and return
Data = NewData;





return
end
