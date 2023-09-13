function OutStruct = spawn_uniform_struct(Fields,Size)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spawn a struct containing many identically-shaped NaN fields
%Corwin Wright, c.wright@bath.ac.uk, 2023/09/06 
%
%inputs:
%  Fields    -cell array containing list of fields, e.g. {'A','B'}
%  Size      - size of each array
%
%outputs:
%  OutStruct - the output struct
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OutStruct = struct();
for iField=1:1:numel(Fields)

  OutStruct.(Fields{iField}) = NaN(Size);

end


end
