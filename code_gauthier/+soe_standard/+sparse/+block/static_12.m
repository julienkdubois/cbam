function [y, T, residual, g1] = static_12(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(1, 1);
  T(11)=log(y(48));
  residual(1)=(T(11))-(T(11)*params(51)+x(13));
if nargout > 3
    g1_v = NaN(1, 1);
g1_v(1)=1/y(48)-params(51)*1/y(48);
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 1, 1);
end
end
