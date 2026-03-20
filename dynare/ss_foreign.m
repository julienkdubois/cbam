function [mu_F, p_int_F, mc_F] = ss_foreign(tau_H, tau_F, sig_F, ...
    theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, ...
    l_H, l_F, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F)
% SS_FOREIGN  Solve the foreign intermediate firm steady-state fixed point.
%
% Solves simultaneously for (mu_F, p_int_F) using fsolve.
%
% The two residual equations are:
%
%   R1: FOC abatement
%       theta1_F*theta2_F*mu_F^{theta2_F-1}
%           - tau_F*sig_F
%           - eta*chi_F*((p_int_F - mc_F)/p_tilde_F)*CBAM*sig_F  = 0
%
%   R2: NKPC at SS (Rotemberg = 0):
%       p_int_F - [eta_eff/(eta_eff-1)] * mc_F * PG_F  = 0
%
% All auxiliary quantities (PG_H, PG_F, G_H, G_F, chi_F, eta_eff, mc_F)
% are computed inside the residual function given (mu_F, p_int_F).
%
% Note: p_int_H is pinned by the normalisation p_int_H = PG_H(p_int_H).
% That fixed point is also iterated inside the residual so that PG_F is
% evaluated consistently.

CBAM = tau_H - tau_F;   % = 0 in symmetric baseline, > 0 with CBAM

% ---- initial guesses ----
mu_F0    = (tau_F * sig_F / (theta1_F * theta2_F))^(1/(theta2_F - 1));
p_int_F0 = p_F;

x0 = [mu_F0; p_int_F0];

% ---- fsolve options ----
opts = optimoptions('fsolve', ...
    'Display',       'off',      ...
    'TolFun',        1e-12,      ...
    'TolX',          1e-12,      ...
    'MaxIterations', 2000,       ...
    'MaxFunctionEvaluations', 50000);

[x_sol, fval, exitflag] = fsolve(@(x) residuals(x, ...
    tau_H, tau_F, sig_F, theta1_F, theta2_F, eta, ...
    gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_F, alpha, ...
    zeta, alpha_h, alpha_F, p_H, p_F, CBAM), x0, opts);

if exitflag <= 0
    warning('ss_foreign: fsolve did not converge (exitflag=%d). Residuals: [%.2e, %.2e]', ...
        exitflag, fval(1), fval(2));
end

mu_F    = x_sol(1);
p_int_F = x_sol(2);

% Recover mc_F at the solution
[~, mc_F] = residuals(x_sol, tau_H, tau_F, sig_F, theta1_F, theta2_F, eta, ...
    gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_F, alpha, ...
    zeta, alpha_h, alpha_F, p_H, p_F, CBAM);

fprintf('ss_foreign: mu_F=%.6f  p_int_F=%.6f  mc_F=%.6f\n', mu_F, p_int_F, mc_F);

end % main function


% ================================================================
function [res, mc_F_out] = residuals(x, ...
    tau_H, tau_F, sig_F, theta1_F, theta2_F, eta, ...
    gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_F, alpha, ...
    zeta, alpha_h, alpha_F, p_H, p_F, CBAM)
% Returns [R1; R2] as described above.
% Also returns mc_F as second output (used by caller to recover mc_F).

mu_F    = x(1);
p_int_F = x(2);

% Guard bounds
mu_F    = max(min(mu_F,    1-1e-8), 1e-8);
p_int_F = max(p_int_F, 1e-8);

% e = 1 at SS
e_ss = 1;

% --- Effective export price ---
p_tilde_F = (1+tau_i)*p_int_F + CBAM*sig_F*(1 - mu_F);

% --- PG_H: closed-form from normalisation p_int_H = PG_H ---
% p_int_H^{1-eta}*(1-(1-g_y)*l_H) = g_y*l_F*(p_tilde_F)^{1-eta}
% => p_int_H = (g_y*l_F / (1-(1-g_y)*l_H))^{1/(1-eta)} * p_tilde_F
p_int_H = ( gamma_y_H*l_F / (1 - (1-gamma_y_H)*l_H) )^(1/(1-eta)) * p_tilde_F;
PG_H    = p_int_H;

% --- PG_F ---
PG_F = ( gamma_y_F*l_H*((1+tau_i)*p_int_H/e_ss)^(1-eta) ...
       + (1-gamma_y_F)*l_F*p_int_F^(1-eta) )^(1/(1-eta));

% --- G_H and G_F from FOC of final good firm ---
G_H = (zeta * alpha_h * p_H * l_H^(1-zeta) / PG_H)^(1/(1-zeta));
G_F = (zeta * alpha_F * p_F * l_F^(1-zeta) / PG_F)^(1/(1-zeta));

% --- Per-firm demands (consistent with model.mod resource constraints) ---
y_dom_F = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F / l_F;
y_exp_F = gamma_y_H * l_H * (e_ss * p_tilde_F / PG_H)^(-eta) * G_H / l_F;
y_int_F = y_dom_F + y_exp_F;

if y_int_F <= 0
    res = [1e6; 1e6];
    mc_F_out = NaN;
    return
end

chi_F = y_exp_F / y_int_F;

% --- Effective elasticity ---
eta_eff_F = eta * (1 - chi_F + chi_F * p_int_F*(1+tau_i) / p_tilde_F);
eta_eff_F = max(eta_eff_F, 1 + 1e-8);

% --- mc_F from NKPC at SS: p_int_F = [eta_eff/(eta_eff-1)] * mc_F * PG_F ---
mc_F = p_int_F * (eta_eff_F - 1) / eta_eff_F / PG_F;

% --- Labour and wage implied by mc_F and production function ---
n_F   = (y_int_F / (Gamma_F * l_F^(1-alpha)))^(1/alpha);
w_F   = (mc_F - theta1_F*mu_F^theta2_F - tau_F*sig_F*(1-mu_F)) * alpha * y_int_F / n_F;

% Verify mc consistency (used for R2; these should match by construction)
mc_F_check = (1/alpha)*w_F*(n_F/y_int_F) + theta1_F*mu_F^theta2_F + tau_F*sig_F*(1-mu_F);

% --- Residuals ---
% R1: FOC abatement
markup_ratio = (p_int_F - mc_F) / p_tilde_F;
R1 = theta1_F*theta2_F*mu_F^(theta2_F-1) ...
   - tau_F*sig_F ...
   - eta * chi_F * markup_ratio * CBAM * sig_F;

% R2: NKPC markup pricing (expressed as relative deviation for scaling)
R2 = p_int_F - (eta_eff_F / (eta_eff_F - 1)) * mc_F * PG_F;

res      = [R1; R2];
mc_F_out = mc_F;

end % residuals