function [p_int_H, p_int_F, PG_H, PG_F, mc_H, mc_F] = ss_no_cbam(eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, zeta, alpha_h, alpha_F, p_H, p_F)
% SS_NO_CBAM  Steady-state prices for the no-CBAM model.
%
% Solves the 2x2 fixed point p_int = PG(p_int_H, p_int_F) in log-space
% to avoid conditioning issues with eta > 1 (p^{1-eta} = p^{-9} for eta=10).
%
% Equations (e=1 at SS):
%   R1: p_int_H = PG_H(p_int_H, p_int_F)
%   R2: p_int_F = PG_F(p_int_H, p_int_F)

e_ss = 1;

% Initial guesses in log-space (log prices close to 0)
log_p0 = [log(p_H); log(p_F)];

opts = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxIterations', 5000, 'MaxFunctionEvaluations', 100000);

[log_p_sol, fval, exitflag] = fsolve(@(lp) residuals_log(lp, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, e_ss), log_p0, opts);

if exitflag <= 0
    warning('ss_no_cbam: fsolve did not converge (exitflag=%d). Residuals: [%.2e, %.2e]', exitflag, fval(1), fval(2));
end

p_int_H = exp(log_p_sol(1));
p_int_F = exp(log_p_sol(2));

p_tilde_F = (1+tau_i) * p_int_F / e_ss;
p_tilde_H = (1+tau_i) * p_int_H / e_ss;

PG_H = ( (1-gamma_y_H)*l_H*p_int_H^(1-eta) + gamma_y_H*l_F*p_tilde_F^(1-eta) )^(1/(1-eta));
PG_F = ( gamma_y_F*l_H*p_tilde_H^(1-eta)   + (1-gamma_y_F)*l_F*p_int_F^(1-eta) )^(1/(1-eta));

mc_H = (eta-1)/eta;
mc_F = (eta-1)/eta;

fprintf('ss_no_cbam: p_int_H=%.6f  p_int_F=%.6f  PG_H=%.6f  PG_F=%.6f\n', p_int_H, p_int_F, PG_H, PG_F);

end


function res = residuals_log(log_p, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, e_ss)

p_int_H = exp(log_p(1));
p_int_F = exp(log_p(2));

p_tilde_F = (1+tau_i) * p_int_F / e_ss;
p_tilde_H = (1+tau_i) * p_int_H / e_ss;

PG_H = ( (1-gamma_y_H)*l_H*p_int_H^(1-eta) + gamma_y_H*l_F*p_tilde_F^(1-eta) )^(1/(1-eta));
PG_F = ( gamma_y_F*l_H*p_tilde_H^(1-eta)   + (1-gamma_y_F)*l_F*p_int_F^(1-eta) )^(1/(1-eta));

% Residuals in log-space: log(p_int) - log(PG) = 0
R1 = log_p(1) - log(PG_H);
R2 = log_p(2) - log(PG_F);

res = [R1; R2];

end
