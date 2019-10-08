function A = np2mat(X, tp)
% Convert NumPy ndarray to a Matlab array
if nargin < 2
    tp = 'd';
end
sz = int32(py.array.array('i', X.shape));
if strcmp(tp, 'd')
    A = reshape(double(py.array.array(tp,X.flatten('F'))), sz);
elseif strcmp(tp, 'i')
    A = reshape(int32(py.array.array(tp,X.flatten('F'))), sz);
else
    error('Unknown data type')
end