function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 48);

T = soe_standard.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(33) = getPowerDeriv(y(26)-params(7)*y(1),(-params(1)),1);
T(34) = getPowerDeriv((y(87)-y(26)*params(7))/(y(26)-params(7)*y(1)),(-params(1)),1);
T(35) = (-((-y(26))/(y(1)*y(1))/(y(26)/y(1))));
T(36) = (-(1/y(1)/(y(26)/y(1))));
T(37) = getPowerDeriv(T(24),1-params(20),1);
T(38) = getPowerDeriv(y(33),1-params(33),1);
T(39) = y(26)*params(24)*getPowerDeriv(y(34),(-params(16)),1);
T(40) = getPowerDeriv(y(34)/y(51),(-params(16)),1);
T(41) = getPowerDeriv(y(38)-params(8)*y(7),(-params(2)),1);
T(42) = getPowerDeriv((y(92)-y(38)*params(8))/(y(38)-params(8)*y(7)),(-params(2)),1);
T(43) = getPowerDeriv(T(29),1-params(20),1);
T(44) = getPowerDeriv(y(45),1-params(33),1);
T(45) = getPowerDeriv(y(51)*y(46),(-params(16)),1);
T(46) = T(40)*(-y(34))/(y(51)*y(51));
T(47) = params(31)*getPowerDeriv(y(56),params(32),1);
T(48) = params(31)*getPowerDeriv(y(57),params(32),1);

end
