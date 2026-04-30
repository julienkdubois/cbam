function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = soe_standard.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 48
    T = [T; NaN(48 - size(T, 1), 1)];
end
T(33) = getPowerDeriv(y(62)-params(7)*y(1),(-params(1)),1);
T(34) = getPowerDeriv((y(123)-y(62)*params(7))/(y(62)-params(7)*y(1)),(-params(1)),1);
T(35) = (-((-y(62))/(y(1)*y(1))/(y(62)/y(1))));
T(36) = (-(1/y(1)/(y(62)/y(1))));
T(37) = getPowerDeriv(T(24),1-params(20),1);
T(38) = getPowerDeriv(y(69),1-params(33),1);
T(39) = y(62)*params(24)*getPowerDeriv(y(70),(-params(16)),1);
T(40) = getPowerDeriv(y(70)/y(87),(-params(16)),1);
T(41) = getPowerDeriv(y(74)-params(8)*y(13),(-params(2)),1);
T(42) = getPowerDeriv((y(135)-y(74)*params(8))/(y(74)-params(8)*y(13)),(-params(2)),1);
T(43) = getPowerDeriv(T(29),1-params(20),1);
T(44) = getPowerDeriv(y(81),1-params(33),1);
T(45) = getPowerDeriv(y(87)*y(82),(-params(16)),1);
T(46) = T(40)*(-y(70))/(y(87)*y(87));
T(47) = params(31)*getPowerDeriv(y(92),params(32),1);
T(48) = params(31)*getPowerDeriv(y(93),params(32),1);
end
