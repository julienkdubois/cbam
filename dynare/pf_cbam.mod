% ================================================================
% VARIABLES
% ================================================================

@#ifndef POLICY_HOME
    @#define POLICY_HOME = 0
@#endif

@#ifndef POLICY_FOREIGN
    @#define POLICY_FOREIGN = 0
@#endif

% POLICY_HOME / POLICY_FOREIGN:
% 0 = Taylor rule
% 1 = Ramsey instrument

var
    %--- HOME HOUSEHOLDS ---
    c_H lb_H r_H w_H n_H cpi_H b_F c_H_h c_H_f

    %--- FOREIGN HOUSEHOLDS ---
    c_F lb_F r_F w_F n_F cpi_F b_H c_F_f c_F_h

    %--- HOME PRICES ---
    p_H PG_H p_int_H pi_H

    %--- FOREIGN PRICES ---
    p_F PG_F p_int_F pi_F

    %--- HOME FINAL GOOD SECTOR ---
    y_H G_H

    %--- FOREIGN FINAL GOOD SECTOR ---
    y_F G_F

    %--- HOME INTERMEDIATE SECTOR ---
    y_int_H mc_H mu_H

    %--- FOREIGN INTERMEDIATE SECTOR ---
    y_int_F mc_F mu_F

    %--- EMISSIONS ---
    e_H e_F

    %--- INTERNATIONAL ---
    de rer NFA_H NFA_F

    %--- TRADE ---
    ex_H ex_F

    %--- POLICIES ---
    tau_H tau_F T_H T_F

    %--- SHOCKS / AUX ---
    e_z_H e_z_F e_p_H e_p_F e_r_H e_r_F e_x_H e_x_F e_t_H e_t_F e_e
    ;

varexo 
    cbam % Transition variable for CBAM policy
    eta_z_H eta_z_F eta_p_H eta_p_F eta_r_H eta_r_F eta_x_H eta_x_F eta_t_H eta_t_F eta_e
    ;

% ================================================================
% PARAMETERS
% ================================================================

parameters
    l_H l_F sigmaC_H sigmaC_F sigmaH_H sigmaH_F beta
    psi_H psi_F % Fixed structural parameters
    phi gamma_c_H gamma_c_F psi_B alpha_h alpha_F zeta
    eta gamma_y_H gamma_y_F alpha Gamma_H Gamma_F
    sig_H sig_F theta1_H theta1_F theta2_H theta2_F
    tau_H_ss tau_F_ss kappa_H kappa_F tau_i
    rho phi_pi phi_y pi_star p_H_ss p_F_ss
    rho_z_H rho_z_F rho_p_H rho_p_F rho_r_H rho_r_F
    rho_x_H rho_x_F rho_t_H rho_t_F rho_e
    ;

% --- Calibration ---
l_H       = 0.4;
l_F       = 0.6;
sigmaC_H  = 1.5;
sigmaC_F  = 1.5;
sigmaH_H  = 2.0;
sigmaH_F  = 2.0;
beta      = 0.994;
phi       = 2.0;
gamma_c_H = 0.1;
gamma_c_F = 0.1;
psi_B     = 0.007;
zeta      = 0.33;
alpha_h   = 1.0;
alpha_F   = 1.0;
eta       = 2.0;
gamma_y_H = 0.3;
gamma_y_F = 0.3;
alpha     = 0.7;
Gamma_H   = 1.0;
Gamma_F   = 1.0;
sig_H     = 0.2;
sig_F     = 0.2;
theta1_H  = 0.05;
theta1_F  = 0.05;
theta2_H  = 2.6;
theta2_F  = 2.6;
tau_H_ss  = 0.1;
tau_F_ss  = 0.02;
kappa_H   = 100;
kappa_F   = 100;
tau_i     = 0.05;
rho       = 0.8;
phi_pi    = 1.5;
phi_y     = 0.05;
pi_star   = 1;
p_H_ss    = 1;
p_F_ss    = 1;
rho_z_H   = 0.95;
rho_z_F   = 0.95;
rho_p_H   = 0.95;
rho_p_F   = 0.95;
rho_r_H   = 0.4;
rho_r_F   = 0.4;
rho_x_H   = 0.4;
rho_x_F   = 0.4;
rho_t_H   = 0.4;
rho_t_F   = 0.4;
rho_e     = 0.1;

% --- PRE-CALIBRATION OF PSI ---
[~, ~, ~, ~, ~, ~, ~, nH_init, nF_init, cH_init, cF_init, wH_init, wF_init] = ss_pf(0, tau_H_ss, tau_F_ss, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H_ss, p_F_ss, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, 1, 1);
psi_H = (cH_init^(-sigmaC_H) * wH_init) / (nH_init^sigmaH_H);
psi_F = (cF_init^(-sigmaC_F) * wF_init) / (nF_init^sigmaH_F);

% ================================================================
% MODEL
% ================================================================

model;
% --- Local definitions ---
#e_t = rer * cpi_H / cpi_F;
#p_tilde_F = (1+tau_i)*p_int_F + cbam*sig_F*(1-mu_F);
#y_dom_F   = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F / l_F;
#y_exp_F   = gamma_y_H * l_H * (e_t * p_tilde_F / PG_H)^(-eta) * G_H / l_F;
#chi_F     = y_exp_F / (y_dom_F + y_exp_F);
#eta_eff_F = eta * (1 - chi_F + chi_F * p_int_F*(1+tau_i) / p_tilde_F);

% --- HOUSEHOLDS ---
[name='MU Home']    lb_H = c_H^(-sigmaC_H);
[name='MU Foreign'] lb_F = c_F^(-sigmaC_F);
[name='Euler H']    lb_H = beta * lb_H(+1) * r_H / cpi_H(+1);
[name='Euler F']    lb_F = beta * lb_F(+1) * r_F / cpi_F(+1);
[name='Labor H']    psi_H * n_H^sigmaH_H = lb_H * w_H;
[name='Labor F']    psi_F * n_F^sigmaH_F = lb_F * w_F;

% --- CONSUMPTION DEMANDS ---
[name='C_H_h'] c_H_h = (1-gamma_c_H) * (cpi_H / p_H)^phi * c_H;
[name='C_H_f'] c_H_f = gamma_c_H * (cpi_H / (e_t * p_F))^phi * c_H;
[name='C_F_f'] c_F_f = (1-gamma_c_F) * (cpi_F / p_F)^phi * c_F;
[name='C_F_h'] c_F_h = gamma_c_F * (cpi_F / (p_H / e_t))^phi * c_F;

% --- INDICES & DYNAMICS ---
[name='CPI H'] cpi_H^(1-phi) = (1-gamma_c_H) * p_H^(1-phi) + gamma_c_H * (e_t * p_F)^(1-phi);
[name='CPI F'] cpi_F^(1-phi) = (1-gamma_c_F) * p_F^(1-phi) + gamma_c_F * (p_H / e_t)^(1-phi);
[name='Price Dyn H'] p_H / p_H(-1) = pi_H / cpi_H;
[name='Price Dyn F'] p_F / p_F(-1) = pi_F / cpi_F;
[name='RER Dyn']     rer / rer(-1) = de * cpi_F / cpi_H;

% --- FINAL GOOD SECTOR ---
[name='Y Final H'] y_H = alpha_h * G_H^zeta * l_H^(1-zeta);
[name='Y Final F'] y_F = alpha_F * G_F^zeta * l_F^(1-zeta);
[name='FOC G_H']   PG_H = zeta * p_H * y_H / G_H;
[name='FOC G_F']   PG_F = zeta * p_F * y_F / G_F;
[name='Price PG H'] PG_H^(1-eta) = (1-gamma_y_H) * l_H * p_int_H^(1-eta) + gamma_y_H * l_F * (e_t * p_tilde_F)^(1-eta);
[name='Price PG F'] PG_F^(1-eta) = gamma_y_F * l_H * ((1+tau_i) * p_int_H / e_t)^(1-eta) + (1-gamma_y_F) * l_F * p_int_F^(1-eta);

% --- INTERMEDIATE SECTOR ---
[name='Prod Int H'] y_int_H = e_z_H * Gamma_H * l_H^(1-alpha) * n_H^alpha;
[name='Prod Int F'] y_int_F = e_z_F * Gamma_F * l_F^(1-alpha) * n_F^alpha;
[name='Abatement H'] theta1_H * theta2_H * mu_H^(theta2_H-1) = tau_H * sig_H;
[name='Abatement F'] theta1_F * theta2_F * mu_F^(theta2_F-1) = tau_F * sig_F + eta * chi_F * ((p_int_F - mc_F) / p_tilde_F) * cbam * sig_F;
[name='MC H'] mc_H = (1/alpha) * w_H * (n_H / y_int_H) + theta1_H * mu_H^theta2_H + tau_H * sig_H * (1 - mu_H);
[name='MC F'] mc_F = (1/alpha) * w_F * (n_F / y_int_F) + theta1_F * mu_F^theta2_F + tau_F * sig_F * (1 - mu_F);
[name='NKPC H'] kappa_H * (pi_H - pi_star) * pi_H = (1-eta) * (p_int_H / PG_H) + eta * e_p_H * mc_H + beta * kappa_H * (pi_H(+1) - pi_star) * pi_H(+1) * y_int_H(+1) / y_int_H;
[name='NKPC F'] kappa_F * (pi_F - pi_star) * pi_F = (1-eta_eff_F) * (p_int_F / PG_F) + eta_eff_F * e_p_F * mc_F + beta * kappa_F * (pi_F(+1) - pi_star) * pi_F(+1) * y_int_F(+1) / y_int_F;

% --- POLICIES ---
[name='Emissions H'] e_H = l_H * sig_H * (1 - mu_H) * y_int_H;
[name='Emissions F'] e_F = l_F * sig_F * (1 - mu_F) * y_int_F;
@#if POLICY_HOME == 0
[name='Taylor H'] r_H = r_H(-1)^rho * (STEADY_STATE(r_H) * (cpi_H / pi_star)^phi_pi * (y_H / STEADY_STATE(y_H))^phi_y)^(1-rho) * e_r_H;
@#endif
@#if POLICY_FOREIGN == 0
[name='Taylor F'] r_F = r_F(-1)^rho * (STEADY_STATE(r_F) * (cpi_F / pi_star)^phi_pi * (y_F / STEADY_STATE(y_F))^phi_y)^(1-rho) * e_r_F;
@#endif
[name='Gov H']    T_H = tau_H * e_H + e_t * cbam * sig_F * (1-mu_F) * gamma_y_H * l_F * (e_t * p_tilde_F / PG_H)^(-eta) * G_H;
[name='Gov F']    T_F = tau_F * e_F;
[name='Tax H']    tau_H = tau_H_ss * e_t_H;
[name='Tax F']    tau_F = tau_F_ss * e_t_F;

% --- MARKET CLEARING ---
[name='Resources Final H'] y_H = l_H * c_H_h + l_F * c_F_h + (kappa_H/2) * (pi_H - pi_star)^2 * y_H;
[name='Resources Final F'] y_F = l_F * c_F_f + l_H * c_H_f + (kappa_F/2) * (pi_F - pi_star)^2 * y_F;
[name='Resources Int H']   y_int_H = (1-gamma_y_H) * (p_int_H / PG_H)^(-eta) * G_H + gamma_y_F * ((1+tau_i) * p_int_H / e_t / PG_F)^(-eta) * G_F * l_F / l_H;
[name='Resources Int F']   y_int_F = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F + gamma_y_H * (e_t * p_tilde_F / PG_H)^(-eta) * G_H * l_H / l_F;

% --- INTERNATIONAL ---
[name='Exports Home'] ex_H = gamma_y_F * l_F * ((1+tau_i)*p_int_H/e_t/PG_F)^(-eta) * G_F / l_H;
[name='Exports Foreign'] ex_F = gamma_y_H * l_H * (e_t*p_tilde_F/PG_H)^(-eta) * G_H / l_F;

% --- INTERNATIONAL FINANCE ---
[name='NFA Accumulation'] NFA_H = r_F(-1) / cpi_H * (rer/rer(-1) * cpi_H/cpi_F) * NFA_H(-1) + p_H * l_F * c_F_h - e_t * p_F * l_H * c_H_f + p_int_H * ex_H * l_H/l_F - e_t * (1+tau_i) * p_int_F * ex_F * l_F/l_H;

[name='NFA Clearing'] l_H * NFA_H + l_F * NFA_F = 0;

[name='UIP Condition'] (rer(+1)*cpi_H(+1)/cpi_F(+1)) / (rer*cpi_H/cpi_F) = (1 + psi_B * (NFA_H - STEADY_STATE(NFA_H))) * r_H / r_F / e_e;

[name='NFA Identity'] NFA_H = e_t * b_F;

[name='Bond Clearing'] l_H * b_F + l_F * b_H = 0;

% --- Shocks ---
log(e_z_H) = rho_z_H * log(e_z_H(-1)) + eta_z_H;
log(e_z_F) = rho_z_F * log(e_z_F(-1)) + eta_z_F;
log(e_p_H) = rho_p_H * log(e_p_H(-1)) + eta_p_H;
log(e_p_F) = rho_p_F * log(e_p_F(-1)) + eta_p_F;
log(e_r_H) = rho_r_H * log(e_r_H(-1)) + eta_r_H;
log(e_r_F) = rho_r_F * log(e_r_F(-1)) + eta_r_F;
log(e_x_H) = rho_x_H * log(e_x_H(-1)) + eta_x_H;
log(e_x_F) = rho_x_F * log(e_x_F(-1)) + eta_x_F;
log(e_t_H) = rho_t_H * log(e_t_H(-1)) + eta_t_H;
log(e_t_F) = rho_t_F * log(e_t_F(-1)) + eta_t_F;
log(e_e)   = rho_e   * log(e_e(-1))   + eta_e;

end;

% ================================================================
% STEADY STATE
% ================================================================

steady_state_model;
    % Interest rate
    r_H_ss = pi_star / beta;
    r_F_ss = pi_star / beta;

    % Unified numeric steady-state solve
    [mu_F, p_int_H, p_int_F, PG_H, PG_F, mc_H, mc_F, n_H, n_F, c_H, c_F, w_H, w_F, G_H, G_F] = ss_pf(cbam, tau_H_ss, tau_F_ss, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H_ss, p_F_ss, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H, psi_F, r_H_ss, r_F_ss);

    % Nominal block
    cpi_H = pi_star; cpi_F = pi_star; pi_H = pi_star; pi_F = pi_star;
    rer = 1; de = 1;
    % Rates from Euler-consistent SS
    r_H = r_H_ss; r_F = r_F_ss;
    tau_H = tau_H_ss; tau_F = tau_F_ss; p_H = p_H_ss; p_F = p_F_ss;

    % Abatement and production
    mu_H = (tau_H_ss * sig_H / (theta1_H * theta2_H))^(1/(theta2_H - 1));
    y_int_H = Gamma_H * l_H^(1-alpha) * n_H^alpha;
    y_int_F = Gamma_F * l_F^(1-alpha) * n_F^alpha;
    
    y_H = alpha_h * G_H^zeta * l_H^(1-zeta);
    y_F = alpha_F * G_F^zeta * l_F^(1-zeta);
    lb_H = c_H^(-sigmaC_H); lb_F = c_F^(-sigmaC_F);

    % Consumption components
    c_H_h = (1-gamma_c_H) * (cpi_H/p_H)^phi * c_H;
    c_H_f = gamma_c_H     * (cpi_H/p_F)^phi * c_H;
    c_F_f = (1-gamma_c_F) * (cpi_F/p_F)^phi * c_F;
    c_F_h = gamma_c_F     * (cpi_F/p_H)^phi * c_F;

    % Trade and NFA
    p_tilde_F_ss = (1+tau_i)*p_int_F + cbam*sig_F*(1 - mu_F);
    ex_H = gamma_y_F * l_F * ((1+tau_i)*p_int_H/PG_F)^(-eta) * G_F / l_H;
    ex_F = gamma_y_H * l_H * (p_tilde_F_ss/PG_H)^(-eta) * G_H / l_F;

    NFA_H = (p_H*l_F*c_F_h - p_F*l_H*c_H_f + p_int_H*ex_H*l_H/l_F - (1+tau_i)*p_int_F*ex_F*l_F/l_H) / (1 - r_F/cpi_H);
    NFA_F = -(l_H/l_F)*NFA_H; b_F = NFA_H; b_H = NFA_F;

    % Emissions, fiscal variables, and shocks
    e_H = l_H * sig_H * (1-mu_H) * y_int_H;
    e_F = l_F * sig_F * (1-mu_F) * y_int_F;
    T_H = tau_H * e_H + cbam * sig_F * (1-mu_F) * gamma_y_H * l_F * (p_tilde_F_ss/PG_H)^(-eta) * G_H;
    T_F = tau_F * e_F;
    e_z_H=1; e_z_F=1; e_p_H=1; e_p_F=1; e_r_H=1; e_r_F=1; e_x_H=1; e_x_F=1; e_t_H=1; e_t_F=1; e_e=1;
end;

% ================================================================
% TRANSITION SETUP
% ================================================================

@#if POLICY_HOME == 0 && POLICY_FOREIGN == 0

initval;
cbam = 0;
end;
steady;

endval;
cbam = tau_H_ss - tau_F_ss;
end;
steady;

shocks;
var cbam; periods 13:40; values (tau_H_ss - tau_F_ss);
end;

perfect_foresight_setup(periods=40);
perfect_foresight_solver(maxit=400, tolf=1e-8);

@#else

initval;
    cbam  = 0;
    r_H   = pi_star / beta;
    r_F   = pi_star / beta;
    tau_H = tau_H_ss;
    tau_F = tau_F_ss;
    p_H   = p_H_ss;
    p_F   = p_F_ss;
    pi_H  = pi_star;
    pi_F  = pi_star;
    cpi_H = pi_star;
    cpi_F = pi_star;
    rer   = 1;
    de    = 1;
    NFA_H = 0;
    NFA_F = 0;
    b_F   = 0;
    b_H   = 0;
    e_z_H = 1; e_z_F = 1; e_p_H = 1; e_p_F = 1;
    e_r_H = 1; e_r_F = 1; e_x_H = 1; e_x_F = 1;
    e_t_H = 1; e_t_F = 1; e_e   = 1;
end;

@#if POLICY_HOME == 1 && POLICY_FOREIGN == 0
% CORRECTION: planner_discount = beta (was 0.99, inconsistent with beta=0.994)
planner_objective(log(c_H) - psi_H*(n_H^(1+sigmaH_H))/(1+sigmaH_H));
ramsey_model(planner_discount=beta, instruments=(r_H));
@#elseif POLICY_HOME == 0 && POLICY_FOREIGN == 1
planner_objective(log(c_F) - psi_F*(n_F^(1+sigmaH_F))/(1+sigmaH_F));
ramsey_model(planner_discount=beta, instruments=(r_F));
@#elseif POLICY_HOME == 1 && POLICY_FOREIGN == 1
planner_objective((log(c_H) - psi_H*(n_H^(1+sigmaH_H))/(1+sigmaH_H)) + (log(c_F) - psi_F*(n_F^(1+sigmaH_F))/(1+sigmaH_F)));
ramsey_model(planner_discount=beta, instruments=(r_H, r_F));
@#endif

@#if POLICY_HOME == 1 && POLICY_FOREIGN == 0
% Continuation path: start from a milder external block, then return to baseline.
gamma_y_H_target = gamma_y_H;
gamma_y_F_target = gamma_y_F;
psi_B_target     = psi_B;

check;

shocks;
var cbam; periods 13:40; values (tau_H_ss - tau_F_ss);
end;

perfect_foresight_setup(periods=40);
perfect_foresight_solver(maxit=600, tolf=1e-8);

@#endif

/*
% ================================================================
% PLOTS - ONE FIGURE PER VARIABLE
% ================================================================
rplot tau_H tau_F;
rplot mu_H mu_F;
rplot e_H e_F;
rplot y_H y_F;
rplot c_H c_F;
rplot n_H n_F;
rplot w_H w_F;
rplot p_int_H p_int_F;
rplot PG_H PG_F;
rplot rer;
rplot ex_H ex_F;
rplot NFA_H NFA_F;
rplot mc_H mc_F;
rplot r_H r_F;
rplot pi_H pi_F;
rplot T_H T_F;
*/
