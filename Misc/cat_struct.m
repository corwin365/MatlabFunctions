function StructA = cat_struct(StructA,StructB,Dimension,VarsToIgnore)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% concatenate all variables in two structs, with specified exceptions
%
%
%StructA and StructB are two structure containing identically-formatted fields (with specified exceptions)
%dimension is the dimension to concatenate in
%VarsToIgnore is a cell array listing any variables we do not want to do this to
%
%Corwin Wright, c.wright@bath.ac.uk, 2002/JUN/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3; Dimension = 1; end
if ~exist('VarsToIgnore','var'); VarsToIgnore = {' '}; end


Fields = fieldnames(StructA);
for iField=1:1:numel(Fields);
  
  %skip specified variables
  if any(strcmp(Fields{iField},VarsToIgnore)); continue; end

  %extract variables
  F = StructA.(Fields{iField});
  G = StructB.(Fields{iField});  
  
  %if the variables are not the same size in the non-specified dimensions,
  %then nan-pad them to match
  szF = size(F); szG = size(G); 
  if numel(szF) < Dimension; aa = ones(Dimension,1); aa(1:numel(szF)) = szF; szF = aa; end
  if numel(szG) < Dimension; aa = ones(Dimension,1); aa(1:numel(szG)) = szG; szG = aa; end
  szF2 = szF; szF2(Dimension) = []; 
  szG2 = szG; szG2(Dimension) = [];
  if ~isequal(szF,szG)
    
    %work out the size of the new array
    NewSizeF = []; NewSizeG = [];
    for iDim=1:1:numel(szF);
      if iDim == Dimension; 
        NewSizeF(iDim) = szF(iDim);
        NewSizeG(iDim) = szG(iDim);        
      else
        NewSizeF(iDim) = max([szF(iDim),szG(iDim)]);
        NewSizeG(iDim) = NewSizeF(iDim);
      end
    end
    
    %make new empty arrays
    NewF = NaN(NewSizeF);
    NewG = NaN(NewSizeG);

    %and copy the data over
    %this is clunky but works. try and improve later.
    CommandF = 'NewF('; CommandG = 'NewG(';
    for iDim = 1:1:numel(NewSizeF);
      CommandF = [CommandF,'1:',num2str(szF(iDim)),','];
      CommandG = [CommandG,'1:',num2str(szG(iDim)),','];  
    end
    CommandF = CommandF(1:end-1); CommandG = CommandG(1:end-1);
    CommandF = [CommandF,') = F;']; CommandG = [CommandG,') = G;'];
    eval(CommandF); eval(CommandG);
    
    F = NewF; G = NewG;
  end
  
    
  
  %concatenate variables
  H = cat(Dimension,F,G);
  StructA.(Fields{iField}) = H;
  
  %done!
  
end

return
