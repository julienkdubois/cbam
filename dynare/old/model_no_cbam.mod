close all;

% ================================================================
% VARIABLES
% ================================================================

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

    %--- SHOCKS ---
    e_z_H e_z_F
    e_p_H e_p_F
    e_r_H e_r_F
    e_x_H e_x_F
    e_t_H e_t_F
    e_e
    ;

varexo
    eta_z_H eta_z_F
    eta_p_H eta_p_F
    eta_r_H eta_r_F
    eta_x_H eta_x_F
    eta_t_H eta_t_F
    eta_e
    ;

parameters
    l_H l_F
    sigmaC_H sigmaC_F
    sigmaH_H sigmaH_F
    beta
    psi_H psi_F
    phi
    gamma_c_H gamma_c_F
    psi_B
    alpha_h alpha_F
    zeta
    eta
    gamma_y_H gamma_y_F
    alpha
    Gamma_H Gamma_F
    sig_H sig_F
    theta1_H theta1_F
    theta2_H theta2_F
    tau_H_ss tau_F_ss
    kappa_H kappa_F
    tau_i
    rho phi_pi phi_y
    pi_star
    p_H_ss p_F_ss
    rho_z_H rho_z_F
    rho_p_H rho_p_F
    rho_r_H rho_r_F
    rho_x_H rho_x_F
    rho_t_H rho_t_F
    rho_e
    ;

l_H       = 0.4;
l_F       = 0.6;
sigmaC_H  = 1.5;
sigmaC_F  = 1.5;
sigmaH_H  = 2.0;
sigmaH_F  = 2.0;
beta      = 0.994;
psi_H     = 1.0;
psi_F     = 1.0;
phi       = 2.0;
gamma_c_H = 0.1;
gamma_c_F = 0.1;
psi_B     = 0.007;
zeta      = 0.7;
alpha_h   = 1.0;
alpha_F   = 1.0;
eta       = 10.0;
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
tau_H_ss  = 0.05;
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

model;

% ----------------------------------------------------------------
% Local variables
% ----------------------------------------------------------------
#e_t = rer * cpi_H / cpi_F;

% No CBAM: effective export price = iceberg cost only
#p_tilde_F = (1+tau_i) * p_int_F;

% ================================================================
% HOUSEHOLDS  [6 equations]
% ================================================================

[name='Marginal utility, Home']
lb_H = c_H^(-sigmaC_H);

[name='Marginal utility, Foreign']
lb_F = c_F^(-sigmaC_F);

[name='Euler equation, Home']
lb_H = beta * lb_H(+1) * r_H / cpi_H(+1);

[name='Euler equation, Foreign']
lb_F = beta * lb_F(+1) * r_F / cpi_F(+1);

[name='Labor supply, Home']
psi_H * n_H^sigmaH_H = lb_H * w_H;

[name='Labor supply, Foreign']
psi_F * n_F^sigmaH_F = lb_F * w_F;

% ================================================================
% CONSUMPTION DEMANDS & PRICE INDICES  [9 equations]
% ================================================================

[name='Demand home goods, Home HH']
c_H_h = (1-gamma_c_H) * (cpi_H / p_H)^phi * c_H;

[name='Demand foreign goods, Home HH']
c_H_f = gamma_c_H * (cpi_H / (e_t * p_F))^phi * c_H;

[name='Demand foreign goods, Foreign HH']
c_F_f = (1-gamma_c_F) * (cpi_F / p_F)^phi * c_F;

[name='Demand home goods, Foreign HH']
c_F_h = gamma_c_F * (cpi_F / (p_H / e_t))^phi * c_F;

[name='CPI index, Home']
cpi_H^(1-phi) = (1-gamma_c_H) * p_H^(1-phi) + gamma_c_H * (e_t * p_F)^(1-phi);

[name='CPI index, Foreign']
cpi_F^(1-phi) = (1-gamma_c_F) * p_F^(1-phi) + gamma_c_F * (p_H / e_t)^(1-phi);

[name='Relative price dynamics, Home']
p_H / p_H(-1) = pi_H / cpi_H;

[name='Relative price dynamics, Foreign']
p_F / p_F(-1) = pi_F / cpi_F;

[name='Real exchange rate']
rer / rer(-1) = de * cpi_F / cpi_H;

% ================================================================
% FINAL GOOD SECTOR  [6 equations]
% ================================================================

[name='Final good production, Home']
y_H = alpha_h * G_H^zeta * l_H^(1-zeta);

[name='Final good production, Foreign']
y_F = alpha_F * G_F^zeta * l_F^(1-zeta);

[name='FOC intermediate bundle, Home']
PG_H = zeta * p_H * y_H / G_H;

[name='FOC intermediate bundle, Foreign']
PG_F = zeta * p_F * y_F / G_F;

% No CBAM: p_tilde_F = (1+tau_i)*p_int_F, symmetric with home
[name='Intermediate bundle price index, Home']
PG_H^(1-eta) = (1-gamma_y_H) * l_H * p_int_H^(1-eta) + gamma_y_H * l_F * (e_t * p_tilde_F)^(1-eta);

[name='Intermediate bundle price index, Foreign']
PG_F^(1-eta) = gamma_y_F * l_H * ((1+tau_i) * p_int_H / e_t)^(1-eta) + (1-gamma_y_F) * l_F * p_int_F^(1-eta);

% ================================================================
% INTERMEDIATE GOOD SECTOR  [8 equations]
% ================================================================

[name='Intermediate production, Home']
y_int_H = e_z_H * Gamma_H * l_H^(1-alpha) * n_H^alpha;

[name='Intermediate production, Foreign']
y_int_F = e_z_F * Gamma_F * l_F^(1-alpha) * n_F^alpha;

% No CBAM: abatement FOC symmetric for both countries
[name='Optimal abatement, Home']
theta1_H * theta2_H * mu_H^(theta2_H-1) = tau_H * sig_H;

[name='Optimal abatement, Foreign']
theta1_F * theta2_F * mu_F^(theta2_F-1) = tau_F * sig_F;

[name='Marginal cost, Home']
mc_H = (1/alpha) * w_H * (n_H / y_int_H) + theta1_H * mu_H^theta2_H + tau_H * sig_H * (1 - mu_H);

[name='Marginal cost, Foreign']
mc_F = (1/alpha) * w_F * (n_F / y_int_F) + theta1_F * mu_F^theta2_F + tau_F * sig_F * (1 - mu_F);

% No CBAM: eta_eff_F = eta (no wedge on foreign demand elasticity)
[name='NKPC, Home']
kappa_H * (pi_H - pi_star) * pi_H = (1-eta) * (p_int_H / PG_H) + eta * e_p_H * mc_H + beta * kappa_H * (pi_H(+1) - pi_star) * pi_H(+1) * y_int_H(+1) / y_int_H;

[name='NKPC, Foreign']
kappa_F * (pi_F - pi_star) * pi_F = (1-eta) * (p_int_F / PG_F) + eta * e_p_F * mc_F + beta * kappa_F * (pi_F(+1) - pi_star) * pi_F(+1) * y_int_F(+1) / y_int_F;

% ================================================================
% EMISSIONS  [2 equations]
% ================================================================

[name='Emissions, Home']
e_H = l_H * sig_H * (1 - mu_H) * y_int_H;

[name='Emissions, Foreign']
e_F = l_F * sig_F * (1 - mu_F) * y_int_F;

% ================================================================
% PUBLIC SECTOR & MONETARY POLICY  [6 equations]
% ================================================================

[name='Taylor rule, Home']
r_H = r_H(-1)^rho * (STEADY_STATE(r_H) * (cpi_H / pi_star)^phi_pi * (y_H / STEADY_STATE(y_H))^phi_y)^(1-rho) * e_r_H;

[name='Taylor rule, Foreign']
r_F = r_F(-1)^rho * (STEADY_STATE(r_F) * (cpi_F / pi_star)^phi_pi * (y_F / STEADY_STATE(y_F))^phi_y)^(1-rho) * e_r_F;

% No CBAM: government collects only domestic carbon tax
[name='Government budget, Home']
T_H = tau_H * e_H;

[name='Government budget, Foreign']
T_F = tau_F * e_F;

[name='Carbon tax, Home']
tau_H = tau_H_ss * e_t_H;

[name='Carbon tax, Foreign']
tau_F = tau_F_ss * e_t_F;

% ================================================================
% GOODS MARKET CLEARING  [4 equations]
% ================================================================

[name='Resource constraint, Final good Home']
y_H = l_H * c_H_h + e_x_F * l_F * c_F_h + (kappa_H/2) * (pi_H - pi_star)^2 * y_H;

[name='Resource constraint, Final good Foreign']
y_F = l_F * c_F_f + e_x_H * l_H * c_H_f + (kappa_F/2) * (pi_F - pi_star)^2 * y_F;

[name='Resource constraint, Intermediate good Home']
y_int_H = (1-gamma_y_H) * (p_int_H / PG_H)^(-eta) * G_H + gamma_y_F * ((1+tau_i) * p_int_H / e_t / PG_F)^(-eta) * G_F * l_F / l_H;

[name='Resource constraint, Intermediate good Foreign']
y_int_F = (1-gamma_y_F) * (p_int_F / PG_F)^(-eta) * G_F + gamma_y_H * (e_t * p_tilde_F / PG_H)^(-eta) * G_H * l_H / l_F;

% ================================================================
% INTERNATIONAL FINANCIAL MARKETS  [5 equations]
% ================================================================

[name='NFA accumulation, Home']
NFA_H = r_F(-1) / cpi_H * de * NFA_H(-1) + p_H * l_F * c_F_h - e_t * p_F * l_H * c_H_f + p_int_H * gamma_y_F * l_F * ((1+tau_i)*p_int_H/e_t/PG_F)^(-eta) * G_F - e_t * (1+tau_i) * p_int_F * gamma_y_H * l_H * (e_t*p_tilde_F/PG_H)^(-eta) * G_H;

[name='International bond market clearing']
l_H * NFA_H + l_F * NFA_F = 0;

[name='UIP, Home']
de(+1) = (1 + psi_B * (NFA_H - STEADY_STATE(NFA_H))) * r_H / r_F / e_e;

[name='NFA identity, Home']
NFA_H = e_t * b_F;

[name='Bond market clearing']
l_H * b_F + l_F * b_H = 0;

% ================================================================
% TRADE  [2 equations]
% ================================================================

[name='Exports, Home']
ex_H = gamma_y_F * l_F * ((1+tau_i)*p_int_H/e_t/PG_F)^(-eta) * G_F / l_H + e_x_F * l_F * c_F_h / l_H;

[name='Exports, Foreign']
ex_F = gamma_y_H * l_H * (e_t*p_tilde_F/PG_H)^(-eta) * G_H / l_F + e_x_H * l_H * c_H_f / l_F;

% ================================================================
% STOCHASTIC PROCESSES  [11 equations]
% ================================================================

[name='TFP shock, Home']
log(e_z_H) = rho_z_H * log(e_z_H(-1)) + eta_z_H;

[name='TFP shock, Foreign']
log(e_z_F) = rho_z_F * log(e_z_F(-1)) + eta_z_F;

[name='Cost-push shock, Home']
log(e_p_H) = rho_p_H * log(e_p_H(-1)) + eta_p_H;

[name='Cost-push shock, Foreign']
log(e_p_F) = rho_p_F * log(e_p_F(-1)) + eta_p_F;

[name='Monetary policy shock, Home']
log(e_r_H) = rho_r_H * log(e_r_H(-1)) + eta_r_H;

[name='Monetary policy shock, Foreign']
log(e_r_F) = rho_r_F * log(e_r_F(-1)) + eta_r_F;

[name='Import demand shock, Home']
log(e_x_H) = rho_x_H * log(e_x_H(-1)) + eta_x_H;

[name='Import demand shock, Foreign']
log(e_x_F) = rho_x_F * log(e_x_F(-1)) + eta_x_F;

[name='Carbon tax shock, Home']
log(e_t_H) = rho_t_H * log(e_t_H(-1)) + eta_t_H;

[name='Carbon tax shock, Foreign']
log(e_t_F) = rho_t_F * log(e_t_F(-1)) + eta_t_F;

[name='UIP shock']
log(e_e) = rho_e * log(e_e(-1)) + eta_e;

end;

% ================================================================
% STEADY STATE
% ================================================================

steady_state_model;

% --- 1. Nominal anchors ---
cpi_H = pi_star;
cpi_F = pi_star;
pi_H  = pi_star;
pi_F  = pi_star;
de    = 1;
rer   = 1;
r_H   = pi_star / beta;
r_F   = pi_star / beta;
tau_H = tau_H_ss;
tau_F = tau_F_ss;

% --- 2. Shocks at SS = 1 ---
e_z_H = 1; e_z_F = 1;
e_p_H = 1; e_p_F = 1;
e_r_H = 1; e_r_F = 1;
e_x_H = 1; e_x_F = 1;
e_t_H = 1; e_t_F = 1;
e_e   = 1;

% --- 3. Numeraire ---
p_H = p_H_ss;
p_F = p_F_ss;
e_ss = 1;

% --- 4. Abatement: analytical for both countries (no CBAM wedge) ---
mu_H = (tau_H_ss * sig_H / (theta1_H * theta2_H))^(1/(theta2_H - 1));
mu_F = (tau_F_ss * sig_F / (theta1_F * theta2_F))^(1/(theta2_F - 1));

% --- 5. Intermediate prices ---
p_int_H = 1;
p_int_F = 1;
p_tilde_F_ss = (1+tau_i) * p_int_F;
p_tilde_H_ss = (1+tau_i) * p_int_H;

% --- 6. PG aggregators ---
PG_H = ( (1-gamma_y_H)*l_H*p_int_H^(1-eta) + gamma_y_H*l_F*p_tilde_F_ss^(1-eta) )^(1/(1-eta));
PG_F = ( gamma_y_F*l_H*p_tilde_H_ss^(1-eta) + (1-gamma_y_F)*l_F*p_int_F^(1-eta) )^(1/(1-eta));

% --- 7. Marginal costs ---
mc_H = (eta-1)/eta * (p_int_H / PG_H);
mc_F = (eta-1)/eta * (p_int_F / PG_F);

% --- 8. Final good sector ---
G_H = (zeta * alpha_h * p_H_ss * l_H^(1-zeta) / PG_H)^(1/(1-zeta));
G_F = (zeta * alpha_F * p_F_ss * l_F^(1-zeta) / PG_F)^(1/(1-zeta));
y_H = alpha_h * G_H^zeta * l_H^(1-zeta);
y_F = alpha_F * G_F^zeta * l_F^(1-zeta);

% --- 9. Intermediate output per firm ---
y_int_H = (1-gamma_y_H)*(p_int_H/PG_H)^(-eta)*G_H + gamma_y_F*((1+tau_i)*p_int_H/e_ss/PG_F)^(-eta)*G_F * l_F/l_H;
y_int_F = (1-gamma_y_F)*(p_int_F/PG_F)^(-eta)*G_F + gamma_y_H*(e_ss*p_tilde_F_ss/PG_H)^(-eta)*G_H * l_H/l_F;

% --- 10. Labour ---
n_H = (y_int_H / (Gamma_H * l_H^(1-alpha)))^(1/alpha);
n_F = (y_int_F / (Gamma_F * l_F^(1-alpha)))^(1/alpha);

% --- 11. Wages ---
w_H = (mc_H - theta1_H*mu_H^theta2_H - tau_H_ss*sig_H*(1-mu_H)) * alpha * y_int_H / n_H;
w_F = (mc_F - theta1_F*mu_F^theta2_F - tau_F_ss*sig_F*(1-mu_F)) * alpha * y_int_F / n_F;

% --- 12. Consumption from goods market clearing ---
A11 = l_H*(1-gamma_c_H)*(cpi_H/p_H_ss)^phi;
A12 = l_F*gamma_c_F    *(cpi_F/p_H_ss)^phi;
A21 = l_H*gamma_c_H    *(cpi_H/p_F_ss)^phi;
A22 = l_F*(1-gamma_c_F)*(cpi_F/p_F_ss)^phi;
c_F = (y_F - A21/A11*y_H) / (A22 - A21*A12/A11);
c_H = (y_H - A12*c_F) / A11;

% --- 13. Consumption demands ---
c_H_h = (1-gamma_c_H) * (cpi_H/p_H_ss)^phi * c_H;
c_H_f = gamma_c_H     * (cpi_H/p_F_ss)^phi * c_H;
c_F_f = (1-gamma_c_F) * (cpi_F/p_F_ss)^phi * c_F;
c_F_h = gamma_c_F     * (cpi_F/p_H_ss)^phi * c_F;

% --- 14. Marginal utilities + psi calibrated residually ---
lb_H  = c_H^(-sigmaC_H);
lb_F  = c_F^(-sigmaC_F);
psi_H = lb_H * w_H / n_H^sigmaH_H;
psi_F = lb_F * w_F / n_F^sigmaH_F;

% --- 15. Emissions and government ---
e_H = l_H * sig_H * (1-mu_H) * y_int_H;
e_F = l_F * sig_F * (1-mu_F) * y_int_F;
T_H = tau_H_ss * e_H;
T_F = tau_F_ss * e_F;

% --- 16. NFA ---
TB_final_H = p_H_ss*l_F*c_F_h - e_ss*p_F_ss*l_H*c_H_f;
TB_inter_H = p_int_H * gamma_y_F*l_F*((1+tau_i)*p_int_H/e_ss/PG_F)^(-eta)*G_F - e_ss*(1+tau_i)*p_int_F * gamma_y_H*l_H*(e_ss*p_tilde_F_ss/PG_H)^(-eta)*G_H;
NFA_H = (TB_final_H + TB_inter_H) / (1 - r_F/cpi_H);
NFA_F = -(l_H/l_F) * NFA_H;
b_F   = NFA_H / e_ss;
b_H   = NFA_F * e_ss;

% --- 17. Exports ---
ex_H = gamma_y_F*l_F*((1+tau_i)*p_int_H/e_ss/PG_F)^(-eta)*G_F/l_H + l_F*c_F_h/l_H;
ex_F = gamma_y_H*l_H*(e_ss*p_tilde_F_ss/PG_H)^(-eta)*G_H/l_F + l_H*c_H_f/l_F;

end;

resid;
check;

% ================================================================
% SHOCKS & SIMULATION
% ================================================================

shocks;
var eta_z_H;  stderr 0.01;
var eta_z_F;  stderr 0.01;
var eta_p_H;  stderr 0.01;
var eta_r_H;  stderr 0.01;
var eta_x_H;  stderr 0.01;
var eta_t_H;  stderr 0.01;
var eta_t_F;  stderr 0.01;
var eta_e;    stderr 0.01;
end;

stoch_simul(order=1, irf=20);
