function [Indices, Values] = closest(Values,Lookup)

[Values,Indices] = min(abs(Values-Lookup));

% % % % % 
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %
% % % % % %N-dimensional min(abs()) - find the index of the closest value
% % % % % %in vector 'Lookup' to every individual value in ND array 'Values'
% % % % % %
% % % % % %inputs: Values - ND array of values to find closest value in Lookup too
% % % % % %        Lookup - 1D vector of Lookup values
% % % % % %
% % % % % %outputs: idx - indices of the closest value
% % % % % %         val - closest value
% % % % % %
% % % % % %(note order of outputs is reversed from min, for backwards compat
% % % % % %with my old 'closest' function which was in 1D and just returned 
% % % % % %idx)
% % % % % %
% % % % % %Corwin Wright, c.wright@bath.ac.uk
% % % % % %2023/MAR/17
% % % % % %
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % 
% % % % % %make the ND values array and the lookup vector the same size and shape
% % % % % A = repmat(Values,[ones(ndims(Values),1)',size(Lookup)]);
% % % % % Lookup = repmat(Lookup',[1,size(Values)]);
% % % % % Lookup = permute(Lookup,[2:ndims(A),1]);
% % % % % 
% % % % % %hence, find the closest value in lookup vector Lookup to each value in ND array Values
% % % % % [Values,Indices] = min(abs(A-Lookup),[],ndims(A));
% % % % % 
% % % % % 

end

