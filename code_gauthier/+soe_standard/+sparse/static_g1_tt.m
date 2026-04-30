function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = soe_standard.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 34
    T = [T; NaN(34 - size(T, 1), 1)];
end
T(25) = getPowerDeriv((y(2))*T(17)*T(18),1-params(20),1);
T(26) = getPowerDeriv(y(8),1-params(33),1);
T(27) = y(1)*params(24)*getPowerDeriv(y(9),(-params(16)),1);
T(28) = getPowerDeriv(y(9)/y(26),(-params(16)),1);
T(29) = getPowerDeriv((y(14))*T(21)*T(22),1-params(20),1);
T(30) = getPowerDeriv(y(20),1-params(33),1);
T(31) = getPowerDeriv(y(26)*y(21),(-params(16)),1);
T(32) = T(28)*(-y(9))/(y(26)*y(26));
T(33) = params(31)*getPowerDeriv(y(31),params(32),1);
T(34) = params(31)*getPowerDeriv(y(32),params(32),1);
end
