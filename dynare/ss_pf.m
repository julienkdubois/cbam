function [mu_F, p_int_H, p_int_F, PG_H, PG_F, mc_H, mc_F, n_H, n_F, c_H, c_F, w_H, w_F, G_H, G_F] = ss_pf(cbam_val, tau_H, tau_F, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H_fixed, psi_F_fixed, r_H_in, r_F_in)
% ss_pf  Steady-state solver for the two-country NK+CBAM model.
% r_H_in and r_F_in are optional and default to pi_star/beta.

    % Default rates when not provided
    if nargin < 36 || isempty(r_H_in)
        r_H_in = pi_star / beta;
    end
    if nargin < 37 || isempty(r_F_in)
        r_F_in = pi_star / beta;
    end

    % Unknowns: [mu_F, p_int_F, p_int_H, n_H, n_F, c_H, c_F]
    x0 = [0.1; 0.7; 0.7; 0.8; 0.8; 0.6; 0.8];
    opts = optimoptions('fsolve', ...
        'Display', 'off', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'MaxFunctionEvaluations', 100000, ...
        'MaxIterations', 5000);

    [sol, ~, exitflag] = fsolve(@(x) global_system(x, cbam_val, tau_H, tau_F, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H_fixed, psi_F_fixed, r_H_in, r_F_in), x0, opts);

    if exitflag <= 0
        warning('ss_pf: fsolve did not converge (exitflag=%d)', exitflag);
    end

    [res_check, ~, ~, ~, ~, ~, ~, ~, ~] = global_system(sol, cbam_val, tau_H, tau_F, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H_fixed, psi_F_fixed, r_H_in, r_F_in);
    max_res = max(abs(res_check));
    if ~isfinite(max_res) || max_res > 1e-8
        error('ss_pf:bad_ss', 'ss_pf residual too large (max |res| = %.3e) for cbam=%.6f', max_res, cbam_val);
    end

    % Extract solution
    mu_F    = sol(1);
    p_int_F = sol(2);
    p_int_H = sol(3);
    n_H     = sol(4);
    n_F     = sol(5);
    c_H     = sol(6);
    c_F     = sol(7);

    [~, PG_H, PG_F, mc_H, mc_F, w_H, w_F, G_H, G_F] = global_system(sol, cbam_val, tau_H, tau_F, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H_fixed, psi_F_fixed, r_H_in, r_F_in);
end


function [res, PG_H, PG_F, mc_H, mc_F, wH, wF, G_H, G_F] = global_system(x, cbam_val, tau_H, tau_F, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H, psi_F, r_H, r_F)

    mu_F    = x(1);
    p_int_F = x(2);
    p_int_H = x(3);
    nH      = x(4);
    nF      = x(5);
    cH      = x(6);
    cF      = x(7);

    % Export effective price and aggregates
    p_tilde_F = (1+tau_i)*p_int_F + cbam_val*sig_F*(1 - mu_F);
    PG_H = ( (1-gamma_y_H)*l_H*p_int_H^(1-eta) + gamma_y_H*l_F*(p_tilde_F)^(1-eta) )^(1/(1-eta));
    PG_F = ( gamma_y_F*l_H*((1+tau_i)*p_int_H)^(1-eta) + (1-gamma_y_F)*l_F*p_int_F^(1-eta) )^(1/(1-eta));

    % Consumption
    cH_h = (1-gamma_c_H) * cH;
    cH_f = gamma_c_H     * cH;
    cF_f = (1-gamma_c_F) * cF;
    cF_h = gamma_c_F     * cF;

    % Final-good sector
    G_H = (zeta * alpha_h * p_H * l_H^(1-zeta) / PG_H)^(1/(1-zeta));
    G_F = (zeta * alpha_F * p_F * l_F^(1-zeta) / PG_F)^(1/(1-zeta));
    y_H = alpha_h * G_H^zeta * l_H^(1-zeta);
    y_F = alpha_F * G_F^zeta * l_F^(1-zeta);

    % Intermediate production
    y_int_H = Gamma_H * l_H^(1-alpha) * nH^alpha;
    y_int_F = Gamma_F * l_F^(1-alpha) * nF^alpha;

    % Labour supply 
    lb_H = cH^(-sigmaC_H);
    lb_F = cF^(-sigmaC_F);
    wH = (psi_H * nH^sigmaH_H) / lb_H;
    wF = (psi_F * nF^sigmaH_F) / lb_F;

    mu_H  = (tau_H * sig_H / (theta1_H * theta2_H))^(1/(theta2_H - 1));
    mc_H  = (1/alpha)*wH*(nH/y_int_H) + theta1_H*mu_H^theta2_H + tau_H*sig_H*(1-mu_H);
    mc_F  = (1/alpha)*wF*(nF/y_int_F) + theta1_F*mu_F^theta2_F + tau_F*sig_F*(1-mu_F);

    % Export share and effective foreign elasticity
    y_dom_F_per = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F / l_F;
    y_exp_F_per = gamma_y_H * l_H * (p_tilde_F / PG_H)^(-eta) * G_H / l_F;
    chi_F       = y_exp_F_per / (y_dom_F_per + y_exp_F_per);
    eta_eff_F   = eta * (1 - chi_F + chi_F * p_int_F*(1+tau_i) / p_tilde_F);

    % 7-equation system
    % 1. Foreign abatement FOC
    res(1,1) = theta1_F*theta2_F*mu_F^(theta2_F-1) - tau_F*sig_F - eta*chi_F*((p_int_F-mc_F)/p_tilde_F)*cbam_val*sig_F;
    % 2. NKPC Foreign
    res(2,1) = (1-eta_eff_F) * (p_int_F / PG_F) + eta_eff_F * mc_F;
    % 3. NKPC Home
    res(3,1) = (1-eta) * (p_int_H / PG_H) + eta * mc_H;
    % 4. Home intermediate market
    res(4,1) = y_int_H - ( (1-gamma_y_H)*(p_int_H/PG_H)^(-eta)*G_H + gamma_y_F*((1+tau_i)*p_int_H/PG_F)^(-eta)*G_F * l_F/l_H );
    % 5. Foreign intermediate market
    res(5,1) = y_int_F - ( (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F + gamma_y_H * (p_tilde_F / PG_H)^(-eta) * G_H * l_H/l_F );
    % 6. Home final-goods market
    res(6,1) = y_H - ( l_H * cH_h + l_F * cF_h );
    % 7. Foreign final-goods market
    res(7,1) = y_F - ( l_F * cF_f + l_H * cH_f );


end