% ================================================================
% cbam_estimation.mod — Corrected version
% ================================================================

% ================================================================
% ENDOGENOUS VARIABLES
% ================================================================
var
    %--- HOME HOUSEHOLDS ---
    c_H lb_H_E lb_H_U D_H r_H w_H n_H cpi_H b_F c_H_h c_H_f

    %--- FOREIGN HOUSEHOLDS ---
    c_F lb_F_E lb_F_U D_F r_F w_F n_F cpi_F b_H c_F_f c_F_h

    %--- HOME PRICES ---
    p_H PG_H p_int_H pi_H

    %--- FOREIGN PRICES ---
    p_F PG_F p_int_F pi_F

    %--- HOME FINAL GOOD SECTOR ---
    y_H G_H

    %--- FOREIGN FINAL GOOD SECTOR ---
    y_F G_F

    %--- HOME INTERMEDIATE SECTOR ---
    y_int_H mc_H mu_H Pi_H

    %--- FOREIGN INTERMEDIATE SECTOR ---
    y_int_F mc_F mu_F Pi_F

    %--- EMISSIONS ---
    e_H e_F

    %--- INTERNATIONAL ---
    de rer NFA_H NFA_F

    %--- TRADE ---
    ex_H ex_F tb_H

    %--- POLICIES ---
    tau_H tau_F T_H T_F

    %--- SHOCK PROCESSES ---
    e_z_H e_z_F e_p_H e_p_F e_r_H e_r_F e_x_H e_x_F e_t_H e_t_F e_e e_b_H e_b_F

    %--- OBSERVABLES ---
    obs_dy_h obs_dy_f
    obs_dc_h obs_dc_f
    obs_pi_h obs_pi_f
    obs_r_h  obs_r_f
    obs_de
    obs_tb_h
    ;

% ================================================================
% EXOGENOUS VARIABLES
% ================================================================
varexo
    cbam
    eta_z_H eta_z_F eta_p_H eta_p_F eta_r_H eta_r_F
    eta_x_H eta_x_F eta_t_H eta_t_F eta_e eta_b_H eta_b_F
    ;

% ================================================================
% PARAMETERS
% ================================================================
parameters
    l_H l_F sigmaC_H sigmaC_F sigmaH_H sigmaH_F beta
    psi_H psi_F omega Dc vartheta
    phi gamma_c_H gamma_c_F psi_B alpha_h alpha_F zeta
    eta gamma_y_H gamma_y_F alpha Gamma_H Gamma_F
    sig_H sig_F theta1_H theta1_F theta2_H theta2_F
    tau_H_ss tau_F_ss kappa_H kappa_F tau_i
    rho_H rho_F phi_pi_H phi_pi_F phi_y_H phi_y_F pi_star p_H_ss p_F_ss
    rho_z_H rho_z_F rho_p_H rho_p_F rho_r_H rho_r_F
    rho_x_H rho_x_F rho_t_H rho_t_F rho_e rho_b_H rho_b_F
    trend_g_H trend_g_F trend_c_H trend_c_F
    pi_bar_H pi_bar_F r_bar_H r_bar_F tb_bar_H de_bar
    ;

% --- Constants and trends (same as pf_cbam) ---
trend_g_H = 0.3429; trend_g_F = 0.3690;
trend_c_H = -0.1845; trend_c_F = -1.6273;
pi_bar_H  = 0.5814; pi_bar_F  = 1.9121;
r_bar_H   = 0.9270;
r_bar_F   = 0.7822;
de_bar    = 0.2541; tb_bar_H  = 0.9228;

% --- Structural calibration (same as pf_cbam) ---
N_total = 1430;
l_H = 450 / N_total;
l_F = 980 / N_total;

sigmaC_H = 1.5;
sigmaC_F = 1.5;
sigmaH_H = 2.0;
sigmaH_F = 2.0;
beta     = 0.994;

psi_B    = 1e-4;

omega    = 0.02;
Dc       = 0.97;
vartheta = 0.05;

%---Calibrated block (quarterly targets — same as pf_cbam) ---
alpha_h   = 10.525429;
alpha_F   = 16.467935;
gamma_y_H = 0.300120;
gamma_y_F = 0.042785;
gamma_c_H = 0.195369;
gamma_c_F = 0.037566;
sig_H     = 2.669698;
sig_F     = 8.375147;

phi   = 3.0;
zeta  = 0.33;
eta   = 3.0;
alpha = 0.7;
Gamma_H = 1.0;
Gamma_F = 1.0;

theta1_H = 0.5;
theta1_F = 0.1;
theta2_H = 2.6;
theta2_F = 2.6;

% Carbon taxes
tau_H_ss = 0.007;
tau_F_ss = 1e-6;

kappa_H = 100;
kappa_F = 100;
tau_i   = 0.05;

% Taylor rule — country-specific parameters
rho_H     = 0.8;
rho_F     = 0.8;
phi_pi_H  = 1.5;
phi_pi_F  = 1.5;
phi_y_H   = 0.05;
phi_y_F   = 0.05;
pi_star = 1.00;
p_H_ss  = 1;
p_F_ss  = 1;

% --- Shock persistence (same as pf_cbam) ---
rho_z_H = 0.95; rho_z_F = 0.95;
rho_p_H = 0.70; rho_p_F = 0.70;
rho_r_H = 0.50; rho_r_F = 0.50;
rho_t_H = 0.70; rho_t_F = 0.70;
rho_e   = 0.80;
rho_x_H = 0.80;
rho_x_F = 0.80;
rho_b_H = 0.80; rho_b_F = 0.80;

% Initialize psi (recalibrated via ss_pf below)
psi_H = 1;
psi_F = 1;

% --- PRE-CALIBRATION OF PSI (same as pf_cbam) ---
[~, ~, ~, ~, ~, ~, ~, nH_init, nF_init, cH_init, cF_init, wH_init, wF_init] = ...
    ss_pf(0, tau_H_ss, tau_F_ss, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, ...
          eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, ...
          alpha_h, alpha_F, p_H_ss, p_F_ss, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, ...
          phi, gamma_c_H, gamma_c_F, pi_star, beta, omega, Dc, vartheta, 1, 1);

psi_H = ( ((cH_init * (1 - omega*Dc)/(1-omega))^(-sigmaC_H)) * wH_init ) / (nH_init^sigmaH_H);
psi_F = ( ((cF_init * (1 - omega*Dc)/(1-omega))^(-sigmaC_F)) * wF_init ) / (nF_init^sigmaH_F);

% ================================================================
% MODEL
% ================================================================
model;
% --- Auxiliary local variables (structurally updated with e_x_F) ---
#e_t       = rer * cpi_H / cpi_F;
#p_tilde_F = (1+tau_i)*p_int_F + cbam*sig_F*(1-mu_F);
#y_dom_F   = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F / l_F;
#y_exp_F   = e_x_F * gamma_y_H * l_H * (e_t * p_tilde_F / PG_H)^(-eta) * G_H / l_F;

% FIX: e_x_F added
#chi_F     = y_exp_F / (y_dom_F + y_exp_F);
#eta_eff_F = eta * (1 - chi_F + chi_F * p_int_F*(1+tau_i) / p_tilde_F);

% ----------------------------------------------------------------
% HOUSEHOLDS (with employed/unemployed heterogeneity)
% ----------------------------------------------------------------
[name='Transfer H']      D_H = Dc * c_H;
[name='Transfer F']      D_F = Dc * c_F;

[name='MU Employed H']   lb_H_E = e_b_H * ((c_H - omega*D_H)/(1-omega))^(-sigmaC_H);
[name='MU Unemployed H'] lb_H_U = e_b_H * D_H^(-sigmaC_H);

[name='MU Employed F']   lb_F_E = e_b_F * ((c_F - omega*D_F)/(1-omega))^(-sigmaC_F);
[name='MU Unemployed F'] lb_F_U = e_b_F * D_F^(-sigmaC_F);

[name='Euler H'] lb_H_E = beta * r_H / cpi_H(+1) * ( (1-omega)*lb_H_E(+1) + omega*lb_H_U(+1) );
[name='Euler F'] lb_F_E = beta * r_F / cpi_F(+1) * ( (1-omega)*lb_F_E(+1) + omega*lb_F_U(+1) );

[name='Labor H'] psi_H * n_H^sigmaH_H = lb_H_E * w_H;
[name='Labor F'] psi_F * n_F^sigmaH_F = lb_F_E * w_F;

% ----------------------------------------------------------------
% CONSUMPTION DEMAND
% ----------------------------------------------------------------
[name='C_H_h'] c_H_h = (1-gamma_c_H) * (cpi_H / p_H)^phi * c_H;
[name='C_H_f'] c_H_f = gamma_c_H * (cpi_H / (e_t * p_F))^phi * c_H;

[name='C_F_f'] c_F_f = (1-gamma_c_F) * (cpi_F / p_F)^phi * c_F;
[name='C_F_h'] c_F_h = gamma_c_F * (cpi_F / (p_H / e_t))^phi * c_F;

% ----------------------------------------------------------------
% PRICE INDICES AND NOMINAL DYNAMICS
% ----------------------------------------------------------------
[name='CPI H']       cpi_H^(1-phi) = (1-gamma_c_H) * p_H^(1-phi) + gamma_c_H * (e_t * p_F)^(1-phi);
[name='CPI F']       cpi_F^(1-phi) = (1-gamma_c_F) * p_F^(1-phi) + gamma_c_F * (p_H / e_t)^(1-phi);

[name='Price Dyn H'] p_H / p_H(-1) = pi_H / cpi_H;
[name='Price Dyn F'] p_F / p_F(-1) = pi_F / cpi_F;

[name='RER Dyn']     rer / rer(-1) = de * cpi_F / cpi_H;

% ----------------------------------------------------------------
% FINAL SECTOR
% ----------------------------------------------------------------
[name='Y Final H'] y_H = alpha_h * G_H^zeta * l_H^(1-zeta);
[name='Y Final F'] y_F = alpha_F * G_F^zeta * l_F^(1-zeta);

[name='FOC G_H']   PG_H = zeta * p_H * y_H / G_H;
[name='FOC G_F']   PG_F = zeta * p_F * y_F / G_F;

[name='Price PG H'] PG_H^(1-eta) = (1-gamma_y_H) * p_int_H^(1-eta) + gamma_y_H * (e_t * p_tilde_F)^(1-eta);
[name='Price PG F'] PG_F^(1-eta) = gamma_y_F * ((1+tau_i) * p_int_H / e_t)^(1-eta) + (1-gamma_y_F) * p_int_F^(1-eta);

% ----------------------------------------------------------------
% INTERMEDIATE SECTOR
% ----------------------------------------------------------------
[name='Prod Int H'] y_int_H = e_z_H * Gamma_H * l_H^(1-alpha) * n_H^alpha;
[name='Prod Int F'] y_int_F = e_z_F * Gamma_F * l_F^(1-alpha) * n_F^alpha;

[name='Abatement H'] theta1_H * theta2_H * mu_H^(theta2_H-1) = tau_H * sig_H;
[name='Abatement F'] theta1_F * theta2_F * mu_F^(theta2_F-1) = tau_F * sig_F + eta * chi_F * ((p_int_F - mc_F) / p_tilde_F) * cbam * sig_F;

[name='MC H'] mc_H = (1/alpha) * w_H * (n_H / y_int_H) + theta1_H * mu_H^theta2_H + tau_H * sig_H * (1 - mu_H);
[name='MC F'] mc_F = (1/alpha) * w_F * (n_F / y_int_F) + theta1_F * mu_F^theta2_F + tau_F * sig_F * (1 - mu_F);

[name='Profit Int H'] Pi_H = p_int_H * y_int_H - w_H * n_H * l_H - theta1_H * mu_H^theta2_H * y_int_H - tau_H * sig_H * (1-mu_H) * y_int_H;
[name='Profit Int F'] Pi_F = p_int_F * y_int_F - w_F * n_F * l_F - theta1_F * mu_F^theta2_F * y_int_F - tau_F * sig_F * (1-mu_F) * y_int_F;

% NKPC
[name='NKPC H'] kappa_H * (pi_H - pi_star) * pi_H = (1-eta) * (p_int_H / PG_H) + eta * e_p_H * mc_H + beta * (1-vartheta) * kappa_H * (pi_H(+1) - pi_star) * pi_H(+1) * y_int_H(+1) / y_int_H;
[name='NKPC F'] kappa_F * (pi_F - pi_star) * pi_F = (1-eta_eff_F) * (p_int_F / PG_F) + eta_eff_F * e_p_F * mc_F + beta * (1-vartheta) * kappa_F * (pi_F(+1) - pi_star) * pi_F(+1) * y_int_F(+1) / y_int_F;

% ----------------------------------------------------------------
% POLICIES
% ----------------------------------------------------------------
[name='Emissions H'] e_H = sig_H * (1 - mu_H) * y_int_H;
[name='Emissions F'] e_F = sig_F * (1 - mu_F) * y_int_F;

[name='Taylor H'] r_H = r_H(-1)^rho_H * (STEADY_STATE(r_H) * (cpi_H / pi_star)^phi_pi_H * (y_H / STEADY_STATE(y_H))^phi_y_H)^(1-rho_H) * e_r_H;
[name='Taylor F'] r_F = r_F(-1)^rho_F * (STEADY_STATE(r_F) * (cpi_F / pi_star)^phi_pi_F * (y_F / STEADY_STATE(y_F))^phi_y_F)^(1-rho_F) * e_r_F;

% FIX: Added e_x_F to CBAM revenue to strictly match actual imports
[name='Gov H'] T_H = tau_H * e_H + e_t * cbam * sig_F * (1-mu_F) * e_x_F * gamma_y_H * l_H * (e_t * p_tilde_F / PG_H)^(-eta) * G_H;
[name='Gov F'] T_F = tau_F * e_F;

[name='Tax H'] tau_H = tau_H_ss * e_t_H;
[name='Tax F'] tau_F = tau_F_ss * e_t_F;

% ----------------------------------------------------------------
% MARKET CLEARING
% ----------------------------------------------------------------
[name='Resources Final H'] y_H = l_H * c_H_h + l_F * c_F_h + (kappa_H/2) * (pi_H - pi_star)^2 * y_H + theta1_H * mu_H^theta2_H * y_int_H + vartheta * Pi_H;
[name='Resources Final F'] y_F = l_F * c_F_f + l_H * c_H_f + (kappa_F/2) * (pi_F - pi_star)^2 * y_F + theta1_F * mu_F^theta2_F * y_int_F + vartheta * Pi_F;

% FIX: Added structural trade shocks e_x_H and e_x_F to intermediate constraints
[name='Resources Int H']   y_int_H = l_H * (1-gamma_y_H) * (p_int_H / PG_H)^(-eta) * G_H + e_x_H * gamma_y_F * ((1+tau_i) * p_int_H / e_t / PG_F)^(-eta) * G_F * l_F;
[name='Resources Int F']   y_int_F = l_F * (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F + e_x_F * gamma_y_H * (e_t * p_tilde_F / PG_H)^(-eta) * G_H * l_H;

% ----------------------------------------------------------------
% INTERNATIONAL TRADE & FINANCE
% ----------------------------------------------------------------
[name='Exports Home']    ex_H = e_x_H * gamma_y_F * l_F * ((1+tau_i)*p_int_H/e_t/PG_F)^(-eta) * G_F / l_H;
[name='Exports Foreign'] ex_F = e_x_F * gamma_y_H * l_H * (e_t*p_tilde_F/PG_H)^(-eta) * G_H / l_F;

[name='Trade Balance H'] tb_H = (p_int_H * ex_H * l_H - e_t * (1+tau_i) * p_int_F * ex_F * l_F + p_H * l_F * c_F_h - e_t * p_F * l_H * c_H_f) / y_H;

[name='NFA Accumulation']
    NFA_H = (r_F(-1) / cpi_H * (rer/rer(-1) * cpi_H/cpi_F)) * NFA_H(-1)
          + p_H * l_F * c_F_h - e_t * p_F * l_H * c_H_f
          + p_int_H * ex_H * l_H - e_t * (1+tau_i) * p_int_F * ex_F * l_F;

[name='NFA Clearing']  NFA_H + NFA_F = 0;
[name='UIP Condition'] (rer(+1)*cpi_H(+1)/cpi_F(+1)) / (rer*cpi_H/cpi_F) = (1 + psi_B * (NFA_H - STEADY_STATE(NFA_H))) * r_H / r_F / e_e;

[name='NFA Identity']  NFA_H = e_t * b_F * l_H;
[name='Bond Clearing'] l_H * b_F + l_F * b_H = 0;

% ----------------------------------------------------------------
% SHOCKS
% ----------------------------------------------------------------
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
log(e_b_H) = rho_b_H * log(e_b_H(-1)) + eta_b_H;
log(e_b_F) = rho_b_F * log(e_b_F(-1)) + eta_b_F;

% ----------------------------------------------------------------
% MEASUREMENT EQUATIONS
% ----------------------------------------------------------------
obs_dy_h = trend_g_H + 100 * (y_H / y_H(-1) - 1);
obs_dy_f = trend_g_F + 100 * (y_F / y_F(-1) - 1);
obs_dc_h = trend_c_H + 100 * (c_H / c_H(-1) - 1);
obs_dc_f = trend_c_F + 100 * (c_F / c_F(-1) - 1);
obs_pi_h = pi_bar_H  + 100 * (pi_H - 1);
obs_pi_f = pi_bar_F  + 100 * (pi_F - 1);
obs_r_h  = r_bar_H   + 100 * (r_H - 1);
obs_r_f  = r_bar_F   + 100 * (r_F - 1);

% FIX: Sign inverted (-) so Eurostat appreciation data matches Model depreciation (rer up)
obs_de   = de_bar    - 100 * (rer / rer(-1) - 1);
obs_tb_h = tb_bar_H + 100 * (tb_H - STEADY_STATE(tb_H));
end;

% ================================================================
% STEADY STATE
% ================================================================
steady_state_model;
    r_H_ss = (pi_star / beta) / (1 - omega + omega * ((Dc * (1 - omega)) / (1 - omega * Dc))^(-sigmaC_H));
    r_F_ss = (pi_star / beta) / (1 - omega + omega * ((Dc * (1 - omega)) / (1 - omega * Dc))^(-sigmaC_F));

    [mu_F, p_int_H, p_int_F, PG_H, PG_F, mc_H, mc_F, n_H, n_F, c_H, c_F, w_H, w_F, G_H, G_F, NFA_H] = 
        ss_pf(cbam, tau_H_ss, tau_F_ss, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, 
              eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, 
              alpha_h, alpha_F, p_H_ss, p_F_ss, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, 
              phi, gamma_c_H, gamma_c_F, pi_star, beta, omega, 
              Dc, vartheta, psi_H, psi_F, r_H_ss, r_F_ss);

    cpi_H = pi_star;
    cpi_F = pi_star;
    pi_H  = pi_star; pi_F  = pi_star;
    rer   = 1;
    de    = 1;
    r_H   = r_H_ss;
    r_F   = r_F_ss;
    tau_H = tau_H_ss;
    tau_F = tau_F_ss;
    p_H   = p_H_ss;
    p_F   = p_F_ss;
    
    mu_H    = (tau_H_ss * sig_H / (theta1_H * theta2_H))^(1/(theta2_H - 1));
    y_int_H = Gamma_H * l_H^(1-alpha) * n_H^alpha;
    y_int_F = Gamma_F * l_F^(1-alpha) * n_F^alpha;
    y_H     = alpha_h * G_H^zeta * l_H^(1-zeta);
    y_F     = alpha_F * G_F^zeta * l_F^(1-zeta);

    D_H    = Dc * c_H;
    D_F    = Dc * c_F;
    lb_H_E = ((c_H - omega*D_H)/(1-omega))^(-sigmaC_H);
    lb_H_U = D_H^(-sigmaC_H);
    lb_F_E = ((c_F - omega*D_F)/(1-omega))^(-sigmaC_F);
    lb_F_U = D_F^(-sigmaC_F);

    Pi_H = p_int_H * y_int_H - w_H * n_H * l_H - theta1_H * mu_H^theta2_H * y_int_H - tau_H * sig_H * (1-mu_H) * y_int_H;
    Pi_F = p_int_F * y_int_F - w_F * n_F * l_F - theta1_F * mu_F^theta2_F * y_int_F - tau_F * sig_F * (1-mu_F) * y_int_F;
    
    c_H_h = (1-gamma_c_H) * (cpi_H/p_H)^phi * c_H;
    c_H_f = gamma_c_H     * (cpi_H/p_F)^phi  * c_H;
    c_F_f = (1-gamma_c_F) * (cpi_F/p_F)^phi  * c_F;
    c_F_h = gamma_c_F     * (cpi_F/p_H)^phi  * c_F;
    
    p_tilde_F_ss = (1+tau_i)*p_int_F + cbam*sig_F*(1 - mu_F);
    ex_H = gamma_y_F * l_F * ((1+tau_i)*p_int_H/PG_F)^(-eta) * G_F / l_H;
    ex_F = gamma_y_H * l_H * (p_tilde_F_ss/PG_H)^(-eta) * G_H / l_F;

    NFA_F = -NFA_H;
    b_F   = NFA_H / l_H;
    b_H   = -b_F * l_H / l_F;
    
    e_H = sig_H * (1-mu_H) * y_int_H;
    e_F = sig_F * (1-mu_F) * y_int_F;
    T_H = tau_H * e_H + cbam * sig_F * (1-mu_F) * gamma_y_H * l_H * (p_tilde_F_ss/PG_H)^(-eta) * G_H;
    T_F = tau_F * e_F;

    tb_H = (p_int_H * ex_H * l_H - (1+tau_i) * p_int_F * ex_F * l_F + p_H * l_F * c_F_h - p_F * l_H * c_H_f) / y_H;
    
    e_z_H=1; e_z_F=1; e_p_H=1; e_p_F=1;
    e_r_H=1; e_r_F=1; e_x_H=1; e_x_F=1;
    e_t_H=1; e_t_F=1; e_e=1; e_b_H=1; e_b_F=1;

    obs_dy_h = trend_g_H;
    obs_dy_f = trend_g_F;
    obs_dc_h = trend_c_H;
    obs_dc_f = trend_c_F;
    obs_pi_h = pi_bar_H + 100*(pi_H - 1);
    obs_pi_f = pi_bar_F + 100*(pi_F - 1);
    obs_r_h  = r_bar_H  + 100*(r_H - 1);
    obs_r_f  = r_bar_F  + 100*(r_F - 1);
    obs_de   = de_bar;
    obs_tb_h = tb_bar_H;   
end;

steady;
check;

% ================================================================
% ESTIMATION
% ================================================================
shocks;
    var eta_z_H = 0.01^2;
    var eta_z_F = 0.01^2;
    var eta_p_H = 0.005^2; 
    var eta_p_F = 0.005^2;
    var eta_r_H = 0.002^2; 
    var eta_r_F = 0.002^2;
    var eta_t_H = 0.005^2; 
    var eta_t_F = 0.005^2;
    var eta_e   = 0.005^2;
    var eta_x_H = 0.01^2;
    var eta_x_F = 0.01^2;
    var eta_b_H = 0.01^2;
    var eta_b_F = 0.01^2;
end;

varobs obs_dy_h obs_dy_f obs_pi_h obs_pi_f obs_r_h obs_r_f obs_tb_h obs_de;

estimated_params;
    kappa_H, 100, inv_gamma_pdf, 100, 50;
    kappa_F, 100, inv_gamma_pdf, 100, 50;
    
    phi, 3.0, inv_gamma_pdf, 3.0, 1.5;
    eta, 3.0, inv_gamma_pdf, 3.0, 1.5;

    rho_H, 0.8, beta_pdf, 0.8, 0.05;
    phi_pi_H, 1.5, inv_gamma_pdf, 1.5, 0.75;
    phi_y_H, 0.05, inv_gamma_pdf, 0.05, 0.025;
    
    stderr eta_z_H, 0.010, inv_gamma_pdf, 0.010, 0.05;
    stderr eta_z_F, 0.010, inv_gamma_pdf, 0.010, 0.10;

    // stderr eta_p_H, 0.005, inv_gamma_pdf, 0.005, 0.003;
    // stderr eta_p_F, 0.005, inv_gamma_pdf, 0.005, 0.003;
    
    stderr eta_r_H, 0.002, inv_gamma_pdf, 0.002, 0.001;
    stderr eta_r_F, 0.002, inv_gamma_pdf, 0.002, 0.001;
    
    stderr eta_e,   0.005, inv_gamma_pdf, 0.005, 0.0025;
    
    stderr eta_x_H, 0.010, inv_gamma_pdf, 0.010, 0.01;
    // stderr eta_x_F, 0.010, inv_gamma_pdf, 0.010, 0.01;

    stderr eta_b_H, 0.010, inv_gamma_pdf, 0.010, 0.005;
    stderr eta_b_F, 0.010, inv_gamma_pdf, 0.010, 0.005;
end;

estimation(
    datafile='../data/clean/dynare_data_eu_row.csv',
    first_obs=1,
    mode_check,
    mode_compute=4,
    mh_replic=50000,
    mh_nblocks=2,
    mh_drop=0.4
) obs_dy_h obs_dy_f obs_pi_h obs_pi_f obs_r_h obs_r_f obs_de obs_tb_h;