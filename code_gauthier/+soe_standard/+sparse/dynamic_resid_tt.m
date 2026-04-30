function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 32
    T = [T; NaN(32 - size(T, 1), 1)];
end
T(1) = params(5)*params(12)*((y(123)-y(62)*params(7))/(y(62)-params(7)*y(1)))^(-params(1));
T(2) = T(1)*y(126);
T(3) = params(5)*params(13)*((y(135)-y(74)*params(8))/(y(74)-params(8)*y(13)))^(-params(2));
T(4) = T(3)*y(138);
T(5) = T(4)*(y(138)-params(36))*y(142);
T(6) = y(68)^params(6);
T(7) = y(80)^params(6);
T(8) = y(69)^(1-params(33));
T(9) = y(81)^(1-params(33));
T(10) = params(31)*y(92)^params(32);
T(11) = y(69)^(-params(33));
T(12) = params(31)*y(93)^params(32);
T(13) = y(81)^(-params(33));
T(14) = params(24)*y(70)^(-params(16));
T(15) = y(62)*T(14);
T(16) = (y(70)/y(87))^(-params(16));
T(17) = params(12)*0.5*(y(65)-params(36))^2;
T(18) = params(25)*y(82)^(-params(16));
T(19) = (y(87)*y(82))^(-params(16));
T(20) = params(13)*0.5*(y(77)-params(36))^2;
T(21) = y(2)^params(20);
T(22) = (steady_state(2))*(y(64)/(steady_state(3)))^params(21);
T(23) = (y(69)/(steady_state(8)))^params(22);
T(24) = T(22)*T(23);
T(25) = T(24)^(1-params(20));
T(26) = y(14)^params(20);
T(27) = (steady_state(14))*(y(76)/(steady_state(15)))^params(21);
T(28) = (y(81)/(steady_state(20)))^params(22);
T(29) = T(27)*T(28);
T(30) = T(29)^(1-params(20));
T(31) = (1-params(19))*(1-params(25))*y(107)/params(19);
T(32) = y(87)/y(26);
end
