function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
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

assert(length(T) >= 34);

T = soe_standard.static_resid_tt(T, y, x, params);

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
