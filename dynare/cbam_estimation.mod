% ================================================================
% VARIABLES ENDOGÈNES
% ================================================================
var
    %--- HOUSEHOLDS HOME (EU) ---
    c_H    % household consumption (home)
    lb_H   % Lagrange multiplier (home budget constraint)
    r_H    % nominal interest rate (home)
    w_H    % real wage (home)
    n_H    % hours worked / labour supply (home)
    cpi_H  % consumer price index (home)
    b_F    % foreign bond position (home liabilities)
    c_H_h  % consumption of home-produced goods by home households
    c_H_f  % consumption of foreign-produced goods by home households

    %--- HOUSEHOLDS FOREIGN (ROW) ---
    c_F    % household consumption (foreign/ROW)
    lb_F   % Lagrange multiplier (foreign budget constraint)
    r_F    % nominal interest rate (foreign)
    w_F    % real wage (foreign)
    n_F    % hours worked / labour supply (foreign)
    cpi_F  % consumer price index (foreign)
    b_H    % home bond position (foreign liabilities)
    c_F_f  % consumption of foreign-produced goods by foreign households
    c_F_h  % consumption of home-produced goods by foreign households

    %--- PRICES ---
    p_H     % price level (home)
    PG_H    % price index / aggregate price (home)
    p_int_H % internal price / producer price (home)
    pi_H    % inflation factor (home)

    p_F     % price level (foreign)
    PG_F    % price index / aggregate price (foreign)
    p_int_F % internal price / producer price (foreign)
    pi_F    % inflation factor (foreign)

    %--- PRODUCTION ---
    y_H   % output / GDP (home)
    G_H   % intermediate productivity (home)
    y_F   % output / GDP (foreign)
    G_F   % intermediate productivity (foreign)

    y_int_H % intermediate / tradable sector output (home)
    mc_H    % marginal cost (home)
    mu_H    % abatement (home)

    y_int_F % intermediate / tradable sector output (foreign)
    mc_F    % marginal cost (foreign)
    mu_F    % abatement (foreign)

    %--- ENVIRONMENT & INTERNATIONAL ---
    e_H  % emissions / environmental externality (home)
    e_F  % emissions / environmental externality (foreign)

    de   % change in real exchange rate / de-trending term
    rer  % real exchange rate
    NFA_H % net foreign assets (home)
    NFA_F % net foreign assets (foreign)

    ex_H % exports of home goods
    ex_F % exports of foreign goods

    %--- POLICIES ---
    tau_H % carbon / tariff parameter (home)
    tau_F % carbon / tariff parameter (foreign)
    T_H   % tax revenue (home)
    T_F   % tax revenue (foreign)

    %--- SHOCK PROCESSES ---
    e_z_H % productivity shock (home)
    e_z_F % productivity shock (foreign)
    e_p_H % price shock (home)
    e_p_F % price shock (foreign)
    e_r_H % monetary / interest shock (home)
    e_r_F % monetary / interest shock (foreign)
    e_t_H % technology / tax shock (home)
    e_t_F % technology / tax shock (foreign)
    e_e   % exchange rate shock / FX shock
    e_x_H % demand / external shock affecting home
    e_x_F % demand / external shock affecting foreign

    e_b_H % household-specific shock (home, e.g., preference)
    e_b_F % household-specific shock (foreign)
    e_ex_H % export shock (home)

    %--- OBSERVABLES ---
    obs_dy_h  % observed output growth (home)
    obs_dy_f  % observed output growth (foreign)
    obs_dc_h  % observed consumption growth (home)
    obs_dc_f  % observed consumption growth (foreign)
    obs_pi_h  % observed inflation (home)
    obs_pi_f  % observed inflation (foreign)
    obs_r_h   % observed interest rate (home)
    obs_r_f   % observed interest rate (foreign)
    obs_de    % observed real exchange rate change
    obs_tb_h  % observed trade balance / net exports (home)
    ;

% ================================================================
% VARIABLES EXOGÈNES
% ================================================================
varexo
    cbam      % carbon border adjustment policy parameter (policy shock / level)
    %--- exogenous shocks (eta_*) are shock innovations for AR(1) processes ---
    eta_z_H   % productivity shock innovation (home)
    eta_z_F   % productivity shock innovation (foreign)
    eta_p_H   % price shock innovation (home)
    eta_p_F   % price shock innovation (foreign)
    eta_r_H   % monetary / interest shock innovation (home)
    eta_r_F   % monetary / interest shock innovation (foreign)
    eta_t_H   % technology / tax shock innovation (home)
    eta_t_F   % technology / tax shock innovation (foreign)
    eta_e     % exchange rate shock innovation
    eta_x_H   % external demand shock innovation (home)
    eta_x_F   % external demand shock innovation (foreign)
    eta_b_H   % household preference / idiosyncratic shock (home)
    eta_b_F   % household preference / idiosyncratic shock (foreign)
    eta_ex_H  % export shock innovation (home)
    ;

% ================================================================
% PARAMÈTRES ET CALIBRATION
% ================================================================
parameters
    l_H l_F sigmaC_H sigmaC_F sigmaH_H sigmaH_F beta
    psi_H psi_F phi gamma_c_H gamma_c_F psi_B alpha_h alpha_F zeta
    eta gamma_y_H gamma_y_F alpha Gamma_H Gamma_F
    sig_H sig_F theta1_H theta1_F theta2_H theta2_F
    tau_H_ss tau_F_ss kappa_H kappa_F tau_i
    rho_H rho_F phi_pi_H phi_pi_F phi_y_H phi_y_F pi_star p_H_ss p_F_ss
    rho_z_H rho_z_F rho_p_H rho_p_F rho_r_H rho_r_F
    rho_t_H rho_t_F rho_e rho_x_H rho_x_F
    rho_b_H rho_b_F rho_ex_H     
    trend_g_H trend_g_F trend_c_H trend_c_F  
    pi_bar_H pi_bar_F r_bar_H r_bar_F tb_bar_H de_bar
    ;

% --- Constantes et Trends ---
trend_g_H = 0.3492; trend_g_F = 2.468;
trend_c_H = 0.276386; trend_c_F = 0.2499;
pi_bar_H  = 0.5828; pi_bar_F  = 0.7943;
r_bar_H   = 0.963618; r_bar_F   = 1.4454;
de_bar    = 0.3596; tb_bar_H  = 1.0883;

% --- Structural calibration ---
l_H=0.15;    % population / weight for home 
l_F=0.85;    % population / weight for foreign
sigmaC_H=1.5; % consumption CRRA (risk aversion) - home
sigmaC_F=1.5; % consumption CRRA (risk aversion) - foreign
sigmaH_H=2.0; % labour supply curvature / disutility parameter - home
sigmaH_F=2.0; % labour supply curvature / disutility parameter - foreign
beta=0.994;   % discount factor
psi_B=0.007;  % cost of holding foreign assets

phi=3.0;        % Armington elasticity of subsitution in consumption (between home and foreign goods)
gamma_c_H=0.17; % consumption share of foreign goods  - home
gamma_c_F=0.03;  % consumption share of foreign goods - foreign
gamma_y_H=0.3;  % openness / imported input share in production - home
gamma_y_F=0.03; % openness / imported input share in production - foreign


%Calibrated parameters to match steady-state targets for output and emissions
alpha_h = 327.421909;
alpha_F = 231.112095;
sig_H = 4.0;
sig_F = 1.0;

zeta=0.33;     % production elasticity on the G/technology component
eta=3.0;       % Armington elasticity in production (between home and foreign intermediate goods)
alpha=0.7;     % production labor share (Cobb-Douglas exponent on labor)
Gamma_H=1.0;   % productivity scale - home
Gamma_F=1.0;   % productivity scale - foreign


theta1_H=0.1; % abatement cost function parameter (scale) - home
theta1_F=0.1; % abatement cost function parameter (scale) - foreign

theta2_H=2.6;  % abatement cost function curvature - home
theta2_F=2.6;  % abatement cost function curvature - foreign
tau_H_ss=0.03;  % steady-state carbon tax rate - home
tau_F_ss=1e-6; % steady-state carbon tax rate - foreign

kappa_H=100;   % Phillips curve / price stickiness slope - home
kappa_F=100;   % Phillips curve / price stickiness slope - foreign

tau_i=0.05;    % iceberg cost rate

rho_H=0.8;     % monetary policy AR(1) persistence - home
rho_F=0.8;     % monetary policy AR(1) persistence - foreign
phi_pi_H=1.5;  % Taylor rule inflation response - home
phi_pi_F=1.5;  % Taylor rule inflation response - foreign
phi_y_H=0.05;  % Taylor rule output response - home
phi_y_F=0.05;  % Taylor rule output response - foreign
pi_star=1.00;  % steady-state inflation target
p_H_ss=1;      % steady-state price level (normalization) - home
p_F_ss=1;      % steady-state price level (normalization) - foreign


% --- Persistance des Chocs ---
rho_z_H  = 0.95;
rho_z_F  = 0.95; 
rho_p_H  = 0.70; 
rho_p_F  = 0.70; 
rho_r_H  = 0.50;
rho_r_F  = 0.50;
rho_t_H  = 0.70;
rho_t_F  = 0.70; 
rho_e    = 0.80;
rho_x_H  = 0.80;
rho_x_F  = 0.80; 
rho_b_H  = 0.80;
rho_b_F  = 0.80; 
rho_ex_H = 0.70;

% Initialisation psi
psi_H = 1; 
psi_F = 1;

% ================================================================
% MODÈLE
% ================================================================
model;
#e_t = rer * cpi_H / cpi_F; 
#p_tilde_F = (1+tau_i)*p_int_F + cbam*sig_F*(1-mu_F); 
#y_dom_F   = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F / l_F; 
#y_exp_F   = gamma_y_H * l_H * (e_t * p_tilde_F / PG_H)^(-eta) * G_H / l_F; 
#chi_F     = y_exp_F / (y_dom_F + y_exp_F); 
#eta_eff_F = eta * (1 - chi_F + chi_F * p_int_F*(1+tau_i) / p_tilde_F); 
#theta1_H_eff = theta1_H * (y_H / l_H);
#theta1_F_eff = theta1_F * (y_F / l_F);

% --- MÉNAGES ---
lb_H = e_b_H * c_H^(-sigmaC_H);
lb_F = e_b_F * c_F^(-sigmaC_F);
lb_H = beta * lb_H(+1) * r_H / cpi_H(+1);
lb_F = beta * lb_F(+1) * r_F / cpi_F(+1);
psi_H * n_H^sigmaH_H = lb_H * w_H;
psi_F * n_F^sigmaH_F = lb_F * w_F;

% --- DEMANDE ---
c_H_h = (1-gamma_c_H) * (cpi_H / p_H)^phi * c_H;
c_H_f = e_x_H * gamma_c_H * (cpi_H / (e_t * p_F))^phi * c_H;
c_F_f = (1-gamma_c_F) * (cpi_F / p_F)^phi * c_F;
c_F_h = e_x_F * gamma_c_F * (cpi_F / (p_H / e_t))^phi * c_F;

% --- PRIX & RER ---
1 = (1-gamma_c_H) * p_H^(1-phi) + gamma_c_H * (rer)^(1-phi);
1 = (1-gamma_c_F) * p_F^(1-phi) + gamma_c_F * (1/rer)^(1-phi);
p_H / p_H(-1) = pi_H / cpi_H;
p_F / p_F(-1) = pi_F / cpi_F;
rer / rer(-1) = de * cpi_F / cpi_H;

% --- PRODUCTION ---
y_H = alpha_h * G_H^zeta * l_H^(1-zeta);
y_F = alpha_F * G_F^zeta * l_F^(1-zeta);
PG_H = zeta * p_H * y_H / G_H;
PG_F = zeta * p_F * y_F / G_F;
PG_H^(1-eta) = (1-gamma_y_H) * l_H * p_int_H^(1-eta) + gamma_y_H * l_F * (e_t * p_tilde_F)^(1-eta);
PG_F^(1-eta) = gamma_y_F * l_H * ((1+tau_i) * p_int_H / e_t)^(1-eta) + (1-gamma_y_F) * l_F * p_int_F^(1-eta);
y_int_H = e_z_H * Gamma_H * l_H^(1-alpha) * n_H^alpha;
y_int_F = e_z_F * Gamma_F * l_F^(1-alpha) * n_F^alpha;

% --- ABATTEMENT & NKPC ---
theta1_H_eff * theta2_H * mu_H^(theta2_H-1) = tau_H * sig_H;
theta1_F_eff * theta2_F * mu_F^(theta2_F-1) = tau_F * sig_F + eta * chi_F * ((p_int_F - mc_F) / p_tilde_F) * cbam * sig_F;
mc_H = (1/alpha) * w_H * (n_H / y_int_H) + theta1_H_eff * mu_H^theta2_H + tau_H * sig_H * (1 - mu_H);
mc_F = (1/alpha) * w_F * (n_F / y_int_F) + theta1_F_eff * mu_F^theta2_F + tau_F * sig_F * (1 - mu_F);
kappa_H * (pi_H - 1) * pi_H = (1-eta) * (p_int_H / PG_H) + eta * mc_H + beta * kappa_H * (pi_H(+1) - 1) * pi_H(+1) * y_int_H(+1) / y_int_H + e_p_H;
kappa_F * (pi_F - 1) * pi_F = (1-eta_eff_F) * (p_int_F / PG_F) + eta_eff_F * mc_F + beta * kappa_F * (pi_F(+1) - 1) * pi_F(+1) * y_int_F(+1) / y_int_F + e_p_F;

% --- POLITIQUES ---
e_H = l_H * sig_H * (1 - mu_H) * y_int_H;
e_F = l_F * sig_F * (1 - mu_F) * y_int_F;
r_H = r_H(-1)^rho_H * (STEADY_STATE(r_H) * (cpi_H / steady_state(cpi_H))^phi_pi_H * (y_H / STEADY_STATE(y_H))^phi_y_H)^(1-rho_H) * e_r_H;
r_F = r_F(-1)^rho_F * (STEADY_STATE(r_F) * (cpi_F / steady_state(cpi_F))^phi_pi_F * (y_F / STEADY_STATE(y_F))^phi_y_F)^(1-rho_F) * e_r_F;
T_H = tau_H * e_H + e_t * cbam * sig_F * (1-mu_F) * gamma_y_H * l_F * (e_t * p_tilde_F / PG_H)^(-eta) * G_H;
T_F = tau_F * e_F;
tau_H = tau_H_ss * e_t_H;
tau_F = tau_F_ss * e_t_F;

% --- MARCHÉS & NFA ---
y_H = l_H * c_H_h + l_F * c_F_h + (kappa_H/2) * (pi_H - 1)^2 * y_H;
y_F = l_F * c_F_f + l_H * c_H_f + (kappa_F/2) * (pi_F - 1)^2 * y_F;
y_int_H = (1-gamma_y_H) * (p_int_H / PG_H)^(-eta) * G_H + gamma_y_F * ((1+tau_i) * p_int_H / e_t / PG_F)^(-eta) * G_F * l_F / l_H;
y_int_F = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F + gamma_y_H * (e_t * p_tilde_F / PG_H)^(-eta) * G_H * l_H / l_F;
ex_H = e_ex_H * gamma_y_F * l_F * ((1+tau_i)*p_int_H/e_t/PG_F)^(-eta) * G_F / l_H;
ex_F = gamma_y_H * l_H * (e_t*p_tilde_F/PG_H)^(-eta) * G_H / l_F;
NFA_H = r_F(-1) / cpi_H * (rer/rer(-1) * cpi_H/cpi_F) * NFA_H(-1) + p_H * l_F * c_F_h - e_t * p_F * l_H * c_H_f + p_int_H * ex_H * l_H/l_F - e_t * (1+tau_i) * p_int_F * ex_F * l_F/l_H;
l_H * NFA_H + l_F * NFA_F = 0;
(rer(+1)*cpi_H(+1)/cpi_F(+1)) / (rer*cpi_H/cpi_F) = (1 + psi_B * (NFA_H - STEADY_STATE(NFA_H))) * r_H / r_F / e_e;

% --- LIAISONS OBLIGATOIRES POUR b_F et b_H ---
NFA_H = e_t * b_F;
l_H * b_F + l_F * b_H = 0;

% --- CHOCS ---
log(e_z_H)=rho_z_H*log(e_z_H(-1))+eta_z_H; log(e_z_F)=rho_z_F*log(e_z_F(-1))+eta_z_F;
e_p_H = rho_p_H * e_p_H(-1) + eta_p_H; 
e_p_F = rho_p_F * e_p_F(-1) + eta_p_F;
log(e_r_H)=rho_r_H*log(e_r_H(-1))+eta_r_H; log(e_r_F)=rho_r_F*log(e_r_F(-1))+eta_r_F;
log(e_t_H)=rho_t_H*log(e_t_H(-1))+eta_t_H; log(e_t_F)=rho_t_F*log(e_t_F(-1))+eta_t_F;
log(e_e)=rho_e*log(e_e(-1))+eta_e;
log(e_x_H)=rho_x_H*log(e_x_H(-1))+eta_x_H; log(e_x_F)=rho_x_F*log(e_x_F(-1))+eta_x_F;
log(e_b_H)=rho_b_H*log(e_b_H(-1))+eta_b_H; log(e_b_F)=rho_b_F*log(e_b_F(-1))+eta_b_F;
log(e_ex_H)=rho_ex_H*log(e_ex_H(-1))+eta_ex_H;

% --- MESURE ---
obs_dy_h=trend_g_H+100*(y_H/y_H(-1)-1); obs_dy_f=trend_g_F+100*(y_F/y_F(-1)-1);
obs_dc_h=trend_c_H+100*(c_H/c_H(-1)-1); obs_dc_f=trend_c_F+100*(c_F/c_F(-1)-1);
obs_pi_h=pi_bar_H+100*(pi_H-1); obs_pi_f=pi_bar_F+100*(pi_F-1);
obs_r_h=r_bar_H+100*(r_H-1); obs_r_f=r_bar_F+100*(r_F-1);
obs_de = de_bar + 100 * (rer / rer(-1) - 1);

#tb_h_ratio = (p_int_H*ex_H*l_H/l_F - p_int_F*ex_F*l_F/l_H + p_H*l_F*c_F_h/l_H - p_F*c_H_f) / y_H;
obs_tb_h = tb_bar_H + 100 * (tb_h_ratio - STEADY_STATE(tb_h_ratio));

end;

% ================================================================
% STEADY STATE
% ================================================================

steady_state_model;
    r_H = pi_star / beta;
    r_F = pi_star / beta;
    rer = 1;
    de = 1;
    cpi_H = r_H * beta;
    cpi_F = r_F * beta;
    pi_H  = cpi_H;
    pi_F  = cpi_F;
    p_H   = ((1 - gamma_c_H * (rer)^(1-phi))/(1-gamma_c_H))^(1/(1-phi));
    p_F   = ((1 - gamma_c_F * (1/rer)^(1-phi))/(1-gamma_c_F))^(1/(1-phi));

        % Récupération des variables via ss_pf 
        [mu_F, p_int_H, p_int_F, PG_H, PG_F, mc_H, mc_F, n_H, n_F, c_H, c_F, w_H, w_F, G_H, G_F] =
                ss_pf(0, tau_H_ss, tau_F_ss, sig_F, sig_H, theta1_H, theta2_H, theta1_F, theta2_F,
                            eta, gamma_y_H, gamma_y_F, tau_i, l_H, l_F, Gamma_H, Gamma_F, alpha, zeta,
                            alpha_h, alpha_F, p_H_ss, p_F_ss, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F,
                            phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H, psi_F, r_H, r_F);

    tau_H = tau_H_ss; tau_F = tau_F_ss;
    y_int_H = Gamma_H * l_H^(1-alpha) * n_H^alpha;
    y_int_F = Gamma_F * l_F^(1-alpha) * n_F^alpha;
    y_H = alpha_h * G_H^zeta * l_H^(1-zeta);
    y_F = alpha_F * G_F^zeta * l_F^(1-zeta);
    mu_H  = (tau_H_ss * sig_H / ((theta1_H * (y_H / l_H)) * theta2_H))^(1/(theta2_H - 1));
    lb_H = c_H^(-sigmaC_H); lb_F = c_F^(-sigmaC_F);

    e_t_ss = rer * cpi_H / cpi_F;
    c_H_h  = (1-gamma_c_H) * (cpi_H / p_H)^phi * c_H;
    c_H_f  = gamma_c_H     * (cpi_H / (e_t_ss * p_F))^phi * c_H;
    c_F_f  = (1-gamma_c_F) * (cpi_F / p_F)^phi * c_F;
    c_F_h  = gamma_c_F     * (cpi_F / (p_H / e_t_ss))^phi * c_F;

    % Définition de p_tilde pour ex_F (Fixe le résidu sur ex_F)
    p_tilde_F_ss = (1+tau_i)*p_int_F;
    ex_H = gamma_y_F * l_F * ((1+tau_i)*p_int_H/e_t_ss/PG_F)^(-eta) * G_F / l_H;
    ex_F = gamma_y_H * l_H * (e_t_ss*p_tilde_F_ss/PG_H)^(-eta) * G_H / l_F;
    
    % Équilibre international
    NFA_H = (p_H*l_F*c_F_h - e_t_ss*p_F*l_H*c_H_f + p_int_H*ex_H*l_H/l_F - e_t_ss*(1+tau_i)*p_int_F*ex_F*l_F/l_H) / (1 - r_F/cpi_H);
    NFA_F = -(l_H/l_F)*NFA_H; 
    b_F = NFA_H;
    b_H = NFA_F;

    e_H = l_H * sig_H * (1-mu_H) * y_int_H; e_F = l_F * sig_F * (1-mu_F) * y_int_F;
    T_H = tau_H * e_H; T_F = tau_F * e_F;
    
    % Chocs unitaires
    e_z_H=1; e_z_F=1; e_p_H=0; e_p_F=0; e_r_H=1; e_r_F=1; e_t_H=1; e_t_F=1; e_e=1;
    e_x_H=1; e_x_F=1; e_b_H=1; e_b_F=1; e_ex_H=1; 

    % Équations de mesure
    obs_dy_h = trend_g_H;
    obs_dy_f = trend_g_F;
    obs_dc_h = trend_c_H;
    obs_dc_f = trend_c_F; 
    obs_pi_h = pi_bar_H + 100*(pi_H - 1);
    obs_pi_f = pi_bar_F + 100*(pi_F - 1);
    obs_r_h  = r_bar_H  + 100*(r_H - 1);
    obs_r_f  = r_bar_F  + 100*(r_F - 1);
    obs_de = de_bar;

    % Calcul du ratio pour obs_tb_h (Identique à la formule du bloc model)
    tb_h_ratio_ss = (p_int_H*ex_H*l_H/l_F - (1/rer)*p_int_F*ex_F*l_F/l_H + p_H*l_F*c_F_h/l_H - (1/rer)*p_F*c_H_f) / y_H;
    obs_tb_h = tb_bar_H + 100 * (tb_h_ratio_ss - tb_h_ratio_ss);
end;

steady;
check;


% ================================================================
% ESTIMATION
% ================================================================
shocks;
    var eta_z_H = 0.01^2; var eta_z_F = 0.01^2;
    var eta_p_H = 0.005^2; var eta_p_F = 0.005^2;
    var eta_r_H = 0.002^2; var eta_r_F = 0.002^2;
    var eta_t_H = 0.005^2; var eta_t_F = 0.005^2;
    var eta_e = 0.005^2; var eta_x_H = 0.01^2; var eta_x_F = 0.01^2;
    var eta_b_H = 0.01^2; var eta_b_F = 0.01^2; var eta_ex_H = 0.01^2;
end;

varobs obs_dy_h obs_dy_f obs_pi_h obs_pi_f obs_r_h obs_r_f obs_de obs_tb_h;

estimated_params;
    % ----------------------------------------------------------------
    % PARAMÈTRES STRUCTURELS 
    % ----------------------------------------------------------------
    % Valeur initiale, Prior, Distribution, Moyenne Prior, Écart-type Prior
    kappa_H, 100, inv_gamma_pdf, 100, 25;
    kappa_F, 100, inv_gamma_pdf, 100, 25;
    
    // eta, 2.0, inv_gamma_pdf, 2.0, 0.5;
    // phi, 3.0, inv_gamma_pdf, 3.0, 0.25;

    // phi_pi_H, 1.5, inv_gamma_pdf, 1.5, 0.5;
    // phi_pi_F, 1.5, inv_gamma_pdf, 1.5,
    phi_y_H, 0.05, inv_gamma_pdf, 0.05, 0.02;
    // phi_y_F, 0.05, inv_gamma_pdf, 0.05, 0.02;
    // rho_H, 0.8, inv_gamma_pdf, 0.8, 0.1;
    // rho_F, 0.8, inv_gamma_pdf, 0.8, 0.1;

    %----------------------------------------------------------------
    % Persistance des chocs (contraints entre 0 et 1)
    %----------------------------------------------------------------
    % Chocs de productivité (Généralement proches de 0.9)
    // rho_z_H, 0.85, 0, 1, beta_pdf, 0.80, 0.10;
    // rho_z_F, 0.85, 0, 1, beta_pdf, 0.80, 0.10;
    
    % Chocs de prix / markup
    // rho_p_H, 0.70, 0, 1, beta_pdf, 0.70, 0.15;
    // rho_p_F, 0.70, 0, 1, beta_pdf, 0.70, 0.15;
    
    % Chocs de politique monétaire (Généralement moins persistants)
    rho_r_H, 0.50, 0, 1, beta_pdf, 0.50, 0.20;
    rho_r_F, 0.50, 0, 1, beta_pdf, 0.50, 0.20;
    
    % Chocs de préférence / demande
    // rho_b_H, 0.70, 0, 1, beta_pdf, 0.70, 0.15;
    // rho_b_F, 0.70, 0, 1, beta_pdf, 0.70, 0.15;
    
    % Chocs de change et externe
    rho_e,    0.80, 0, 1, beta_pdf, 0.80, 0.10;
    rho_ex_H, 0.70, 0, 1, beta_pdf, 0.70, 0.15;

    % ----------------------------------------------------------------
    % ÉCARTS-TYPES DES CHOCS
    % ----------------------------------------------------------------
    
    % Chocs de productivité : obs_dy_h et obs_dy_f
    stderr eta_z_H, 0.010, inv_gamma_pdf, 0.010, 0.05;
    stderr eta_z_F, 0.010, inv_gamma_pdf, 0.010, 0.10;
    
    % Chocs de demande : obs_dc_h
    stderr eta_b_H, 0.010, inv_gamma_pdf, 0.010, 0.005;
    stderr eta_b_F, 0.010, inv_gamma_pdf, 0.010, 0.005;
    
    % Chocs de coûts : obs_pi_h et obs_pi_f
    // stderr eta_p_H, 0.005, inv_gamma_pdf, 0.005, 0.01;
    // stderr eta_p_F, 0.005, inv_gamma_pdf, 0.005, 0.05;
    
    % Chocs monétaires : obs_r_h et obs_r_f
    stderr eta_r_H, 0.002, inv_gamma_pdf, 0.002, 0.001;
    stderr eta_r_F, 0.002, inv_gamma_pdf, 0.002, 0.001;
    
    % Choc de change (Taux de change réel) - Nécessaire pour obs_de
    stderr eta_e,   0.005, inv_gamma_pdf, 0.005, 0.0025;
    
    % Choc d'exportation - Nécessaire pour obs_tb_h
    stderr eta_ex_H, 0.010, inv_gamma_pdf, 0.010, 0.01;
end;

estimation(
    datafile='../data/clean/dynare_data_eu_row.csv',

    first_obs=1,  
    mode_check,

    mode_compute=4, 
    mh_replic=10000, 
    mh_nblocks=2, 
    mh_drop=0.4
) obs_dy_h obs_dy_f  obs_pi_h obs_pi_f obs_r_h obs_r_f obs_de obs_tb_h obs_dc_h obs_dc_f ;
