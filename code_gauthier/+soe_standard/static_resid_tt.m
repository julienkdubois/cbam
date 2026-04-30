function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 24);

T(1) = y(7)^params(6);
T(2) = y(19)^params(6);
T(3) = y(8)^(1-params(33));
T(4) = y(20)^(1-params(33));
T(5) = params(31)*y(31)^params(32);
T(6) = y(8)^(-params(33));
T(7) = params(31)*y(32)^params(32);
T(8) = y(20)^(-params(33));
T(9) = params(24)*y(9)^(-params(16));
T(10) = y(1)*T(9);
T(11) = (y(9)/y(26))^(-params(16));
T(12) = params(12)*0.5*(y(4)-params(36))^2;
T(13) = params(25)*y(21)^(-params(16));
T(14) = (y(26)*y(21))^(-params(16));
T(15) = params(13)*0.5*(y(16)-params(36))^2;
T(16) = y(2)^params(20);
T(17) = (y(3)/(y(3)))^params(21);
T(18) = (y(8)/(y(8)))^params(22);
T(19) = ((y(2))*T(17)*T(18))^(1-params(20));
T(20) = y(14)^params(20);
T(21) = (y(15)/(y(15)))^params(21);
T(22) = (y(20)/(y(20)))^params(22);
T(23) = ((y(14))*T(21)*T(22))^(1-params(20));
T(24) = (1-params(19))*(1-params(25))*y(46)/params(19);

end
