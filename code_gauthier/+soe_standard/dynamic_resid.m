function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = soe_standard.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(61, 1);
    residual(1) = (y(36)) - ((y(26)-params(7)*y(1))^(-params(1)));
    residual(2) = (y(48)) - ((y(38)-params(8)*y(7))^(-params(2)));
    residual(3) = (y(36)) - (params(5)*y(91)*y(27)/y(88));
    residual(4) = (y(48)) - (params(5)*y(96)*y(39)/y(93));
    residual(5) = (params(10)*y(32)^params(3)) - (y(36)*y(31));
    residual(6) = (params(11)*y(44)^params(4)) - (y(48)*y(43));
residual(7) = 1-params(14)+params(14)*y(63)*y(30)-params(12)*y(29)*(y(29)-params(36))+T(2)*(y(89)-params(36))*y(90)/y(33);
residual(8) = 1-params(15)+params(15)*y(69)*y(42)-params(13)*y(41)*(y(41)-params(36))+T(5)/y(33);
    residual(9) = (y(60)) - (y(31)*y(32)/(y(33)*params(6)));
    residual(10) = (y(61)) - (y(43)*y(44)/(params(6)*y(45)));
    residual(11) = (y(33)) - (params(34)*y(62)*T(6));
    residual(12) = (y(45)) - (params(34)*y(68)*T(7));
    residual(13) = (1) - (params(24)*y(34)^(1-params(16))+(1-params(24))*y(51)^(1-params(16)));
    residual(14) = (1) - (params(25)*y(46)^(1-params(16))+(1-params(25))*(1/y(51))^(1-params(16)));
    residual(15) = (y(34)/y(4)) - (y(29)/y(28));
    residual(16) = (y(46)/y(10)) - (y(41)/y(40));
    residual(17) = (y(52)) - (params(19)*params(29)*(1-y(56))*T(8));
    residual(18) = (y(53)) - ((1-params(19))*params(30)*(1-y(57))*T(9));
    residual(19) = (y(60)) - (y(30)-T(10)-(1-y(56))*params(29)*(1-params(33))*y(54)*T(11));
    residual(20) = (y(61)) - (y(42)-T(12)-(1-y(57))*params(30)*(1-params(33))*y(55)*T(13));
    residual(21) = (T(8)*params(29)*y(54)) - (params(31)*params(32)*y(56)^(params(32)-1));
    residual(22) = (T(9)*params(30)*y(55)) - (params(31)*params(32)*y(57)^(params(32)-1));
    residual(23) = (y(33)) - (T(15)+(1-params(19))*y(38)*(1-params(25))*y(71)*T(16)/params(19)+y(58)+y(33)*T(10)+y(33)*T(17)+0.5*params(9)*(y(35)-(steady_state(10)))^2);
    residual(24) = (y(45)) - (y(38)*T(18)+params(19)*y(26)*(1-params(24))*y(65)*T(19)/(1-params(19))+y(59)+y(45)*T(12)+y(45)*T(20)-0.5*params(9)*(y(47)-(steady_state(22)))^2);
    residual(25) = (y(27)) - (T(21)*T(25)*y(64));
    residual(26) = (y(39)) - (T(26)*T(30)*y(70));
    residual(27) = (y(58)) - ((steady_state(8))*params(38)*y(67));
    residual(28) = (y(59)) - ((steady_state(20))*params(39)*y(73));
    residual(29) = (y(54)) - (params(26)*y(66));
    residual(30) = (y(55)) - (params(27)*y(72));
    residual(31) = (y(35)) - (y(8)/y(28)*y(5)*y(50)+y(34)*(T(15)+y(38)*T(16)*T(31))-y(26));
residual(32) = params(19)*y(35)+(1-params(19))*y(47);
    residual(33) = (y(97)) - (y(27)*(1+params(9)*(y(35)-(steady_state(10))))/y(39)/y(86));
    residual(34) = (T(32)) - (y(41)*y(50)/y(29));
    residual(35) = (y(37)) - ((1-params(19))*y(38)*T(16)*(1-params(25))*y(65));
    residual(36) = (y(49)) - (params(19)*y(26)*T(19)*(1-params(24))*y(71));
    residual(37) = (y(74)) - (log(y(33)/y(3)));
    residual(38) = (y(75)) - (log(y(45)/y(9)));
    residual(39) = (y(76)) - (log(y(26)/y(1)));
    residual(40) = (y(77)) - (log(y(26)/y(1)));
    residual(41) = (y(78)) - (y(29)-(steady_state(4)));
    residual(42) = (y(79)) - (y(41)-(steady_state(16)));
    residual(43) = (y(80)) - (y(27)-(steady_state(2)));
    residual(44) = (y(81)) - (y(39)-(steady_state(14)));
    residual(45) = (y(82)) - (log(y(50)));
    residual(46) = (y(83)) - (log(T(32)));
    residual(47) = (y(84)) - (log(y(37)/y(6)));
    residual(48) = (y(85)) - (log(y(49)/y(11)));
    residual(49) = (log(y(62))) - (params(40)*log(y(13))+x(it_, 1));
    residual(50) = (log(y(63))) - (params(42)*log(y(14))+x(it_, 2));
    residual(51) = (log(y(64))) - (params(41)*log(y(15))+x(it_, 3));
    residual(52) = (log(y(65))) - (params(43)*log(y(16))+x(it_, 4));
    residual(53) = (log(y(67))) - (params(43)*log(y(18))+x(it_, 12));
    residual(54) = (log(y(66))) - (params(44)*log(y(17))+x(it_, 10));
    residual(55) = (log(y(68))) - (params(46)*log(y(19))+x(it_, 5));
    residual(56) = (log(y(69))) - (params(48)*log(y(20))+x(it_, 6));
    residual(57) = (log(y(70))) - (params(47)*log(y(21))+x(it_, 7));
    residual(58) = (log(y(71))) - (params(49)*log(y(22))+x(it_, 8));
    residual(59) = (log(y(73))) - (params(51)*log(y(24))+x(it_, 13));
    residual(60) = (log(y(72))) - (params(50)*log(y(23))+x(it_, 11));
    residual(61) = (log(y(86))) - (params(23)*log(y(25))+x(it_, 9));

end
