
% ================================================================
% VARIABLES
% ================================================================

@#ifndef POLICY_HOME
    @#define POLICY_HOME = 0
@#endif

@#ifndef POLICY_FOREIGN
    @#define POLICY_FOREIGN = 0
@#endif

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
    e_z_H e_z_F e_p_H e_p_F e_r_H e_r_F e_t_H e_t_F e_e e_x_H e_x_F
    e_b_H e_b_F
    e_ex_H

    %--- OBSERVABLES ---
    obs_dy_h obs_dy_f obs_dc_h obs_dc_f obs_pi_h obs_pi_f obs_r_h obs_r_f obs_de obs_tb_f
    ;

varexo
    cbam 
    eta_z_H eta_z_F eta_p_H eta_p_F eta_r_H eta_r_F eta_t_H eta_t_F eta_e eta_x_H eta_x_F
    eta_b_H eta_b_F  
    eta_ex_H         
    ;

% ================================================================
% PARAMETERS
% ================================================================

parameters
    l_H l_F sigmaC_H sigmaC_F sigmaH_H sigmaH_F beta
    psi_H psi_F
    phi gamma_c_H gamma_c_F psi_B alpha_h alpha_F zeta
    eta gamma_y_H gamma_y_F alpha Gamma_H Gamma_F
    sig_H sig_F theta1_H theta1_F theta2_H theta2_F
    tau_H_ss tau_F_ss kappa_H kappa_F tau_i
    rho_H rho_F phi_pi_H phi_pi_F phi_y_H phi_y_F pi_star p_H_ss p_F_ss
    rho_z_H rho_z_F rho_p_H rho_p_F rho_r_H rho_r_F
    rho_t_H rho_t_F rho_e rho_x_H rho_x_F
    rho_b_H rho_b_F   
    rho_ex_H          
    trend_g_H trend_g_F trend_c_H trend_c_F  
    pi_bar_H pi_bar_F r_bar_H r_bar_F rer_bar tb_bar_F 
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
gamma_c_H = 0.17;
gamma_c_F = 0.1;
psi_B     = 0.007;
zeta      = 0.33;
alpha_h   = 1.0;
alpha_F   = 1.0;
eta       = 2.0;      
gamma_y_H = 0.3;
gamma_y_F = 0.38;
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
rho_H     = 0.8;
rho_F     = 0.8;
phi_pi_H  = 1.5;
phi_pi_F  = 1.5;
phi_y_H   = 0.05;
phi_y_F   = 0.05;
pi_star   = 1.00;    
p_H_ss    = 1;
p_F_ss    = 1;

% Structural shock persistence
rho_z_H   = 0.95;
rho_z_F   = 0.95;
rho_p_H   = 0.60;
rho_p_F   = 0.60;
rho_r_H   = 0.4;
rho_r_F   = 0.4;
rho_t_H   = 0.4;
rho_t_F   = 0.4;
rho_e     = 0.1;
rho_x_H   = 0.80;
rho_x_F   = 0.80;
rho_b_H   = 0.80;     
rho_b_F   = 0.80;
rho_ex_H  = 0.50;   

% --- Measurement-equation constants ---
pi_bar_H  = 0.448;    % mean obs_pi_h
pi_bar_F  = 0.520;    % mean obs_pi_f
r_bar_H   = 0.168;    % mean obs_r_h
r_bar_F   = 0.410;    % mean obs_r_f
rer_bar   = 1.00;     % mean obs_de
tb_bar_F  = -0.759;  

% --- Growth trends ---
trend_g_H = 0.349;   
trend_g_F = 0.540;    
trend_c_H = 0.265;    
trend_c_F = 0.482;    

% --- PRE-CALIBRATION OF PSI ---
[~, ~, ~, ~, ~, ~, ~, nH_init, nF_init, cH_init, cF_init, wH_init, wF_init] = ss_pf(0, tau_H_ss, tau_F_ss, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F, eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, p_H_ss, p_F_ss, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, phi, gamma_c_H, gamma_c_F, pi_star, beta, kappa_H, kappa_F, 1, 1, pi_star / beta, pi_star / beta);
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
[name='MU Home']    lb_H = e_b_H * c_H^(-sigmaC_H);
[name='MU Foreign'] lb_F = e_b_F * c_F^(-sigmaC_F);
[name='Euler H']    lb_H = beta * lb_H(+1) * r_H / cpi_H(+1);
[name='Euler F']    lb_F = beta * lb_F(+1) * r_F / cpi_F(+1);
[name='Labor H']    psi_H * n_H^sigmaH_H = lb_H * w_H;
[name='Labor F']    psi_F * n_F^sigmaH_F = lb_F * w_F;

% --- CONSUMPTION DEMANDS ---
[name='C_H_h'] c_H_h = (1-gamma_c_H) * (cpi_H / p_H)^phi * c_H;
[name='C_H_f'] c_H_f = e_x_H * gamma_c_H * (cpi_H / (e_t * p_F))^phi * c_H;
[name='C_F_f'] c_F_f = (1-gamma_c_F) * (cpi_F / p_F)^phi * c_F;
[name='C_F_h'] c_F_h = e_x_F * gamma_c_F * (cpi_F / (p_H / e_t))^phi * c_F;

% --- INDICES & DYNAMICS ---
[name='CPI H'] 1 = (1-gamma_c_H) * p_H^(1-phi) + gamma_c_H * (rer)^(1-phi);
[name='CPI F'] 1 = (1-gamma_c_F) * p_F^(1-phi) + gamma_c_F * (1/rer)^(1-phi);
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
[name='NKPC H'] kappa_H * (pi_H - 1) * pi_H = (1-eta) * (p_int_H / PG_H) + eta * e_p_H * mc_H + beta * kappa_H * (pi_H(+1) - 1) * pi_H(+1) * y_int_H(+1) / y_int_H;
[name='NKPC F'] kappa_F * (pi_F - 1) * pi_F = (1-eta_eff_F) * (p_int_F / PG_F) + eta_eff_F * e_p_F * mc_F + beta * kappa_F * (pi_F(+1) - 1) * pi_F(+1) * y_int_F(+1) / y_int_F;

% --- POLICIES ---
[name='Emissions H'] e_H = l_H * sig_H * (1 - mu_H) * y_int_H;
[name='Emissions F'] e_F = l_F * sig_F * (1 - mu_F) * y_int_F;
[name='Taylor H'] r_H = r_H(-1)^rho_H * (STEADY_STATE(r_H) * (cpi_H / steady_state(cpi_H))^phi_pi_H * (y_H / STEADY_STATE(y_H))^phi_y_H)^(1-rho_H) * e_r_H;
[name='Taylor F'] r_F = r_F(-1)^rho_F * (STEADY_STATE(r_F) * (cpi_F / steady_state(cpi_F))^phi_pi_F * (y_F / STEADY_STATE(y_F))^phi_y_F)^(1-rho_F) * e_r_F;
[name='Gov H']    T_H = tau_H * e_H + e_t * cbam * sig_F * (1-mu_F) * gamma_y_H * l_F * (e_t * p_tilde_F / PG_H)^(-eta) * G_H;
[name='Gov F']    T_F = tau_F * e_F;
[name='Tax H']    tau_H = tau_H_ss * e_t_H;
[name='Tax F']    tau_F = tau_F_ss * e_t_F;

% --- MARKET CLEARING ---
[name='Resources Final H'] y_H = l_H * c_H_h + l_F * c_F_h + (kappa_H/2) * (pi_H - 1)^2 * y_H;
[name='Resources Final F'] y_F = l_F * c_F_f + l_H * c_H_f + (kappa_F/2) * (pi_F - 1)^2 * y_F;
[name='Resources Int H']   y_int_H = (1-gamma_y_H) * (p_int_H / PG_H)^(-eta) * G_H + gamma_y_F * ((1+tau_i) * p_int_H / e_t / PG_F)^(-eta) * G_F * l_F / l_H;
[name='Resources Int F']   y_int_F = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F + gamma_y_H * (e_t * p_tilde_F / PG_H)^(-eta) * G_H * l_H / l_F;

% --- INTERNATIONAL ---
[name='Exports Home']    ex_H = e_ex_H * gamma_y_F * l_F * ((1+tau_i)*p_int_H/e_t/PG_F)^(-eta) * G_F / l_H;
[name='Exports Foreign'] ex_F = gamma_y_H * l_H * (e_t*p_tilde_F/PG_H)^(-eta) * G_H / l_F;

% --- INTERNATIONAL FINANCE ---
[name='NFA Accumulation'] NFA_H = r_F(-1) / cpi_H * (rer/rer(-1) * cpi_H/cpi_F) * NFA_H(-1) + p_H * l_F * c_F_h - e_t * p_F * l_H * c_H_f + p_int_H * ex_H * l_H/l_F - e_t * (1+tau_i) * p_int_F * ex_F * l_F/l_H;
[name='NFA Clearing']     l_H * NFA_H + l_F * NFA_F = 0;
[name='UIP Condition']    (rer(+1)*cpi_H(+1)/cpi_F(+1)) / (rer*cpi_H/cpi_F) = (1 + psi_B * (NFA_H - STEADY_STATE(NFA_H))) * r_H / r_F / e_e;
[name='NFA Identity']     NFA_H = e_t * b_F;
[name='Bond Clearing']    l_H * b_F + l_F * b_H = 0;

% --- Shocks ---
log(e_z_H)  = rho_z_H  * log(e_z_H(-1))  + eta_z_H;
log(e_z_F)  = rho_z_F  * log(e_z_F(-1))  + eta_z_F;
log(e_p_H)  = rho_p_H  * log(e_p_H(-1))  + eta_p_H;
log(e_p_F)  = rho_p_F  * log(e_p_F(-1))  + eta_p_F;
log(e_r_H)  = rho_r_H  * log(e_r_H(-1))  + eta_r_H;
log(e_r_F)  = rho_r_F  * log(e_r_F(-1))  + eta_r_F;
log(e_t_H)  = rho_t_H  * log(e_t_H(-1))  + eta_t_H;
log(e_t_F)  = rho_t_F  * log(e_t_F(-1))  + eta_t_F;
log(e_e)    = rho_e    * log(e_e(-1))    + eta_e;
log(e_x_H)  = rho_x_H  * log(e_x_H(-1))  + eta_x_H;
log(e_x_F)  = rho_x_F  * log(e_x_F(-1))  + eta_x_F;
log(e_b_H)  = rho_b_H  * log(e_b_H(-1))  + eta_b_H;
log(e_b_F)  = rho_b_F  * log(e_b_F(-1))  + eta_b_F;
log(e_ex_H) = rho_ex_H * log(e_ex_H(-1)) + eta_ex_H;

% ================================================================
% MEASUREMENT EQUATIONS
% ================================================================

    % --- GDP growth ---
    obs_dy_h = trend_g_H + 100 * (y_H / y_H(-1) - 1);
    obs_dy_f = trend_g_F + 100 * (y_F / y_F(-1) - 1);

    % --- Consumption growth --- [FIX-4] separate trends
    obs_dc_h = trend_c_H + 100 * (c_H / c_H(-1) - 1);
    obs_dc_f = trend_c_F + 100 * (c_F / c_F(-1) - 1);

    % --- Inflation ---
    obs_pi_h = pi_bar_H + 100*(pi_H - 1);
    obs_pi_f = pi_bar_F + 100*(pi_F - 1);

    % --- Interest rates ---
    obs_r_h  = r_bar_H  + 100*(r_H - 1);
    obs_r_f  = r_bar_F  + 100*(r_F - 1);

    % --- Bilateral real exchange rate
    obs_de = rer_bar + 100*(rer - 1);

    % --- Trade Balance Foreign (% GDP) ---
    obs_tb_f = tb_bar_F + 100 * ( (p_int_F*ex_F*l_F/l_H - (1/rer)*p_int_H*ex_H*l_H/l_F + p_F*l_H*c_H_f/l_F - (1/rer)*p_H*c_F_h) / y_F );

end;

% ================================================================
% STEADY STATE
% ================================================================

steady_state_model;
    r_H = pi_star / beta;
    r_F = pi_star / beta;

    rer = 1; de = 1;
    cpi_H = r_H * beta;
    cpi_F = r_F * beta;
    pi_H  = cpi_H;
    pi_F  = cpi_F;
    p_H   = ((1 - gamma_c_H * (rer)^(1-phi))/(1-gamma_c_H))^(1/(1-phi));
    p_F   = ((1 - gamma_c_F * (1/rer)^(1-phi))/(1-gamma_c_F))^(1/(1-phi));

    [mu_F, p_int_H, p_int_F, PG_H, PG_F, mc_H, mc_F, n_H, n_F, c_H, c_F, w_H, w_F, G_H, G_F] =
        ss_pf(0, tau_H_ss, tau_F_ss, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F,
              eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta,
              alpha_h, alpha_F, p_H_ss, p_F_ss, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F,
              phi, gamma_c_H, gamma_c_F, pi_star, beta, kappa_H, kappa_F, psi_H, psi_F, r_H, r_F);

    tau_H = tau_H_ss; tau_F = tau_F_ss;
    mu_H  = (tau_H_ss * sig_H / (theta1_H * theta2_H))^(1/(theta2_H - 1));
    y_int_H = Gamma_H * l_H^(1-alpha) * n_H^alpha;
    y_int_F = Gamma_F * l_F^(1-alpha) * n_F^alpha;
    y_H = alpha_h * G_H^zeta * l_H^(1-zeta);
    y_F = alpha_F * G_F^zeta * l_F^(1-zeta);
    lb_H = c_H^(-sigmaC_H); lb_F = c_F^(-sigmaC_F);  % e_b=1 at steady state

    e_t_ss = rer * cpi_H / cpi_F;
    c_H_h  = (1-gamma_c_H) * (cpi_H / p_H)^phi * c_H;
    c_H_f  = gamma_c_H     * (cpi_H / (e_t_ss * p_F))^phi * c_H;
    c_F_f  = (1-gamma_c_F) * (cpi_F / p_F)^phi * c_F;
    c_F_h  = gamma_c_F     * (cpi_F / (p_H / e_t_ss))^phi * c_F;

    p_tilde_F_ss = (1+tau_i)*p_int_F;   % cbam=0
    % [FIX-3] e_ex_H = 1 at steady state
    ex_H = gamma_y_F * l_F * ((1+tau_i)*p_int_H/e_t_ss/PG_F)^(-eta) * G_F / l_H;
    ex_F = gamma_y_H * l_H * (e_t_ss*p_tilde_F_ss/PG_H)^(-eta) * G_H / l_F;

    NFA_H = (p_H*l_F*c_F_h - e_t_ss*p_F*l_H*c_H_f + p_int_H*ex_H*l_H/l_F - e_t_ss*(1+tau_i)*p_int_F*ex_F*l_F/l_H) / (1 - r_F/cpi_H);
    NFA_F = -(l_H/l_F)*NFA_H; b_F = NFA_H; b_H = NFA_F;

    e_H = l_H * sig_H * (1-mu_H) * y_int_H;
    e_F = l_F * sig_F * (1-mu_F) * y_int_F;
    T_H = tau_H * e_H;
    T_F = tau_F * e_F;

    % Shocks = 1 at steady state
    e_z_H=1; e_z_F=1; e_p_H=1; e_p_F=1;
    e_r_H=1; e_r_F=1; e_t_H=1; e_t_F=1; e_e=1;
    e_x_H=1; e_x_F=1;
    e_b_H=1; e_b_F=1; 
    e_ex_H=1; 

    % --- Measurement equations at steady state ---
    obs_dy_h = trend_g_H;
    obs_dy_f = trend_g_F;
    obs_dc_h = trend_c_H;
    obs_dc_f = trend_c_F; 

    obs_pi_h = pi_bar_H + 100*(cpi_H - 1);  
    obs_pi_f = pi_bar_F + 100*(cpi_F - 1);

    obs_r_h  = r_bar_H  + 100*(r_H - 1);
    obs_r_f  = r_bar_F  + 100*(r_F - 1);

    obs_de = rer_bar + 100*(de - 1);  

    obs_tb_f = tb_bar_F + 100 * ( (p_int_F*ex_F*l_F/l_H - (1/rer)*p_int_H*ex_H*l_H/l_F + p_F*l_H*c_H_f/l_F - (1/rer)*p_H*c_F_h) / y_F );
end;

% ================================================================
% SHOCKS
% ================================================================

shocks;
    var cbam = 0;

    var eta_z_H  = 0.01^2;
    var eta_z_F  = 0.01^2;
    var eta_p_H  = 0.005^2;
    var eta_p_F  = 0.005^2;
    var eta_r_H  = 0.002^2;
    var eta_r_F  = 0.002^2;
    var eta_t_H  = 0.005^2;
    var eta_t_F  = 0.005^2;
    var eta_e    = 0.005^2;
    var eta_x_H  = 0.01^2;
    var eta_x_F  = 0.01^2;
    var eta_b_H  = 0.01^2;   
    var eta_b_F  = 0.01^2;   
    var eta_ex_H = 0.01^2;   
end;

steady;
resid;
check;

% ================================================================
% BAYESIAN ESTIMATION
% ================================================================


varobs obs_dy_h obs_dy_f obs_dc_h obs_dc_f obs_pi_h obs_pi_f obs_r_h obs_r_f obs_tb_f obs_de;

estimated_params;
    % ----------------------------------------------------------------
    % STRUCTURAL PARAMETERS
    % ----------------------------------------------------------------
    % Nominal rigidities
    kappa_H, 100, inv_gamma_pdf, 100, 25;
    kappa_F, 100, inv_gamma_pdf, 100, 25;

    % Armington elasticities 
    eta, 2.0, gamma_pdf, 2.0, 0.75;
    phi, 2.0, gamma_pdf, 2.0, 0.5;

    % Monetary policy
    // rho_H, 0.8, beta_pdf, 0.8, 0.1;          % Taylor-rule persistence (Home)
    // rho_F, 0.8, beta_pdf, 0.8, 0.1;          % Taylor-rule persistence (Foreign)
    phi_pi_H, 1.8, gamma_pdf, 1.8, 0.2;      % inflation response (Home)
    phi_pi_F, 1.8, gamma_pdf, 1.8, 0.2;      % inflation response (Foreign)
    phi_y_H, 0.05, gamma_pdf, 0.05, 0.05;    % output response (Home)
    phi_y_F, 0.05, gamma_pdf, 0.05, 0.05;    % output response (Foreign)

    % ----------------------------------------------------------------
    % SHOCK VARIANCES
    % ----------------------------------------------------------------
    % TFP
    stderr eta_z_H, 0.010, inv_gamma_pdf, 0.010, 0.005;
    stderr eta_z_F, 0.010, inv_gamma_pdf, 0.010, 0.010;

    % Consumption preference
    stderr eta_b_H, 0.010, inv_gamma_pdf, 0.010, 0.005;
    stderr eta_b_F, 0.010, inv_gamma_pdf, 0.010, 0.005;

    % Markup prix
    stderr eta_p_H, 0.005, inv_gamma_pdf, 0.005, 0.0025;
    stderr eta_p_F, 0.005, inv_gamma_pdf, 0.005, 0.0025;

    % Monetary policy
    stderr eta_r_H, 0.002, inv_gamma_pdf, 0.002, 0.001;
    stderr eta_r_F, 0.002, inv_gamma_pdf, 0.002, 0.001;

    % UIP 
    stderr eta_e, 0.005, inv_gamma_pdf, 0.005, 0.0025;

    % Demand export
    stderr eta_ex_H, 0.010, inv_gamma_pdf, 0.010, 0.005;

end;

estimated_params_init(use_calibration);
end;

identification;

estimation(
    datafile   = '../data/clean/dynare_data_bilateral_2002Q1_2019Q4.csv',
    first_obs  = 1,
    nobs       = 68,
    mode_compute = 4,
    mh_jscale = 0.4,
    mode_check,
    mh_drop = 0.4,
    mh_nblocks = 2,
    mh_replic  = 10000
) obs_dy_h obs_dy_f obs_dc_h obs_dc_f obs_pi_h obs_pi_f obs_r_h obs_r_f obs_de;
