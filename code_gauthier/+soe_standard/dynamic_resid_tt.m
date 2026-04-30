function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 32);

T(1) = params(5)*params(12)*((y(87)-y(26)*params(7))/(y(26)-params(7)*y(1)))^(-params(1));
T(2) = T(1)*y(89);
T(3) = params(5)*params(13)*((y(92)-y(38)*params(8))/(y(38)-params(8)*y(7)))^(-params(2));
T(4) = T(3)*y(94);
T(5) = T(4)*(y(94)-params(36))*y(95);
T(6) = y(32)^params(6);
T(7) = y(44)^params(6);
T(8) = y(33)^(1-params(33));
T(9) = y(45)^(1-params(33));
T(10) = params(31)*y(56)^params(32);
T(11) = y(33)^(-params(33));
T(12) = params(31)*y(57)^params(32);
T(13) = y(45)^(-params(33));
T(14) = params(24)*y(34)^(-params(16));
T(15) = y(26)*T(14);
T(16) = (y(34)/y(51))^(-params(16));
T(17) = params(12)*0.5*(y(29)-params(36))^2;
T(18) = params(25)*y(46)^(-params(16));
T(19) = (y(51)*y(46))^(-params(16));
T(20) = params(13)*0.5*(y(41)-params(36))^2;
T(21) = y(2)^params(20);
T(22) = (steady_state(2))*(y(28)/(steady_state(3)))^params(21);
T(23) = (y(33)/(steady_state(8)))^params(22);
T(24) = T(22)*T(23);
T(25) = T(24)^(1-params(20));
T(26) = y(8)^params(20);
T(27) = (steady_state(14))*(y(40)/(steady_state(15)))^params(21);
T(28) = (y(45)/(steady_state(20)))^params(22);
T(29) = T(27)*T(28);
T(30) = T(29)^(1-params(20));
T(31) = (1-params(19))*(1-params(25))*y(71)/params(19);
T(32) = y(51)/y(12);

end
