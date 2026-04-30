function [mu_F, p_int_H, p_int_F, PG_H, PG_F, mc_H, mc_F, n_H, n_F, c_H, c_F, w_H, w_F, G_H, G_F] = ss_pf(cbam_val, tau_H, tau_F, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H_fixed, psi_F_fixed, r_H_in, r_F_in)

    % Valeurs par défaut pour les taux
    if nargin < 35 || isempty(r_H_in), r_H_in = pi_star / beta; end
    if nargin < 36 || isempty(r_F_in), r_F_in = pi_star / beta; end

    % --- Initialisation Robuste ---
    % On estime la consommation par tête initiale
    y_guess_H = alpha_h * l_H; 
    y_guess_F = alpha_F * l_F;
    x0 = [0.1; 1.0; 1.0; 0.33; 0.33; y_guess_H/l_H; y_guess_F/l_F];

    % --- Stratégie Multi-départ avec Transformation ---
    % On utilise pack_state pour transformer le problème en espace non contraint
    candidate_x0 = [x0, [0.05; 0.5; 0.5; 0.2; 0.2; 0.5; 0.5], x0 * 1.2, x0 * 0.8];
    
    algs = {'levenberg-marquardt', 'trust-region-dogleg'};
    best_sol = [];
    best_max_res = Inf;

    for a = 1:numel(algs)
        opts_try = optimoptions('fsolve', 'Display', 'off', 'Algorithm', algs{a}, ...
            'FunctionTolerance', 1e-14, 'StepTolerance', 1e-14);
        
        for s = 1:size(candidate_x0,2)
            z0 = pack_state(candidate_x0(:,s));
            try
                % Résolution dans l'espace transformé (z) pour garantir x > 0
                [z_sol, ~, exitflag] = fsolve(@(z) global_system_transformed(z, cbam_val, tau_H, tau_F, sig_F, sig_H, ...
                    theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, ...
                    Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, ...
                    sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H_fixed, ...
                    psi_F_fixed, r_H_in, r_F_in), z0, opts_try);
                
                sol_try = unpack_state(z_sol);
                [res_vec] = global_system(sol_try, cbam_val, tau_H, tau_F, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H_fixed, psi_F_fixed, r_H_in, r_F_in);
                
                max_res_try = max(abs(res_vec));
                if isfinite(max_res_try) && max_res_try < best_max_res
                    best_max_res = max_res_try;
                    best_sol = sol_try;
                end
                if best_max_res <= 1e-10, break; end
            catch
            end
        end
        if best_max_res <= 1e-10, break; end
    end

    if isempty(best_sol) || best_max_res > 1e-5
        error('ss_pf: Failed to converge (Res = %e). Check sig_H/sig_F values.', best_max_res);
    end

    % --- Extraction des variables finales ---
    sol = best_sol;
    mu_F = sol(1); p_int_F = sol(2); p_int_H = sol(3); n_H = sol(4); n_F = sol(5); c_H = sol(6); c_F = sol(7);
    [~, PG_H, PG_F, mc_H, mc_F, w_H, w_F, G_H, G_F] = global_system(sol, cbam_val, tau_H, tau_F, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H_fixed, psi_F_fixed, r_H_in, r_F_in);
end

% --- Fonctions Auxiliaires de Sécurité ---[cite: 5]
function z = pack_state(x)
    z = zeros(size(x));
    mu = min(max(x(1), 1e-8), 1 - 1e-8);
    z(1) = log(mu / (1 - mu)); % Transformation logit pour mu[cite: 5]
    z(2:end) = log(max(x(2:end), 1e-8)); % Transformation log pour variables > 0[cite: 5]
end

function x = unpack_state(z)
    x = zeros(size(z));
    x(1) = 1 / (1 + exp(-z(1)));
    x(2:end) = exp(z(2:end));
end

function res = global_system_transformed(z, varargin)
    x = unpack_state(z);
    res = global_system(x, varargin{:});
    % Protection contre les valeurs complexes ou infinies[cite: 5]
    if ~all(isreal(res)) || any(~isfinite(res(:))), res = 1e8 * ones(size(res)); end
end

function [res, PG_H, PG_F, mc_H, mc_F, wH, wF, G_H, G_F] = global_system(x, cbam_val, tau_H, tau_F, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H, p_F, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H, psi_F, r_H, r_F)
    mu_F = x(1); p_int_F = x(2); p_int_H = x(3); nH = x(4); nF = x(5); cH = x(6); cF = x(7);

    % Sécurisation des arguments de puissance pour éviter les nombres complexes[cite: 5]
    p_tilde_F = (1+tau_i)*p_int_F + cbam_val*sig_F*(1 - mu_F);
    PG_H = (max((1-gamma_y_H)*l_H*p_int_H^(1-eta) + gamma_y_H*l_F*p_tilde_F^(1-eta), 1e-12))^(1/(1-eta));
    PG_F = (max(gamma_y_F*l_H*((1+tau_i)*p_int_H)^(1-eta) + (1-gamma_y_F)*l_F*p_int_F^(1-eta), 1e-12))^(1/(1-eta));

    G_H = (max(zeta * alpha_h * p_H * l_H^(1-zeta) / PG_H, 1e-12))^(1/(1-zeta));
    G_F = (max(zeta * alpha_F * p_F * l_F^(1-zeta) / PG_F, 1e-12))^(1/(1-zeta));
    y_H = alpha_h * G_H^zeta * l_H^(1-zeta);
    y_F = alpha_F * G_F^zeta * l_F^(1-zeta);

    y_int_H = Gamma_H * l_H^(1-alpha) * nH^alpha;
    y_int_F = Gamma_F * l_F^(1-alpha) * nF^alpha;

    wH = (psi_H * nH^sigmaH_H) / (cH^(-sigmaC_H));
    wF = (psi_F * nF^sigmaH_F) / (cF^(-sigmaC_F));

    theta1_H_eff = theta1_H;
    theta1_F_eff = theta1_F;
    mu_H = (max(tau_H * sig_H / (theta1_H_eff * theta2_H), 1e-12))^(1/(theta2_H - 1));
    
    mc_H = (1/alpha)*wH*(nH/y_int_H) + theta1_H_eff*mu_H^theta2_H + tau_H*sig_H*(1-mu_H);
    mc_F = (1/alpha)*wF*(nF/y_int_F) + theta1_F_eff*mu_F^theta2_F + tau_F*sig_F*(1-mu_F);

    y_exp_F_per = gamma_y_H * l_H * (p_tilde_F / PG_H)^(-eta) * G_H / l_F;
    y_dom_F_per = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F / l_F;
    chi_F = y_exp_F_per / (y_dom_F_per + y_exp_F_per);
    eta_eff_F = eta * (1 - chi_F + chi_F * p_int_F*(1+tau_i) / p_tilde_F);

    % Système d'équations (Résidus)[cite: 5]
    res(1,1) = theta1_F_eff*theta2_F*mu_F^(theta2_F-1) - tau_F*sig_F - eta*chi_F*((p_int_F-mc_F)/p_tilde_F)*cbam_val*sig_F;
    res(2,1) = (1-eta_eff_F) * (p_int_F / PG_F) + eta_eff_F * mc_F;
    res(3,1) = (1-eta) * (p_int_H / PG_H) + eta * mc_H;
    res(4,1) = y_int_H - ( (1-gamma_y_H)*(p_int_H/PG_H)^(-eta)*G_H + gamma_y_F*((1+tau_i)*p_int_H/PG_F)^(-eta)*G_F * l_F/l_H );
    res(5,1) = y_int_F - ( (1-gamma_y_F)*(p_int_F/PG_F)^(-eta)*G_F + gamma_y_H*(p_tilde_F/PG_H)^(-eta)*G_H * l_H/l_F );
    res(6,1) = y_H - ( l_H * (1-gamma_c_H)*cH + l_F * gamma_c_F*cF );
    res(7,1) = y_F - ( l_F * (1-gamma_c_F)*cF + l_H * gamma_c_H*cH );
end