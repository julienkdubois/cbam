close all;

%----------------------------------------------------------------
% 1. Variables
%----------------------------------------------------------------
var
    %--- HOME HOUSEHOLDS ---
    c_H         % aggregate consumption 
    lb_H        % marginal utility (lambda)
    r_H         % interest rate on home bonds
    w_H         % real wage
    n_H         % hours worked
    cpi_H       % CPI inflation (P_t/P_{t-1})
    b_F         % home household stock of foreign bonds
    P_H         % price level from home (CPI)

    %--- FOREIGN HOUSEHOLDS ---
    c_F
    lb_F
    r_F
    w_F
    n_F
    cpi_F
    b_star_H    % foreign household stock of home bonds
    P_F

    %--- HOME FINAL GOOD SECTOR ---
    y_H         % final home output
    G_H         % CES intermediate aggregator
    PG_H        % price index of intermediates
    p_H         % relative price of home final good (p_t(h)/P_t)

    %--- FOREIGN FINAL GOOD SECTOR ---
    y_F
    G_F
    PG_F
    p_F

    %--- HOME INTERMEDIATE SECTOR ---
    mc_H        % marginal cost of home intermediate firms
    pi_H        % PPI home inflation p_t(h)/p_t-1(h)
    mu_H        % abatement rate for home intermediate firms

    %--- FOREIGN INTERMEDIATE SECTOR ---
    mc_F
    pi_F
    mu_F

    %--- EMISSIONS ---
    e_H         % home emissions
    e_F

    %--- INTERNATIONAL ---
    de          % nominal exchange rate growth e_t/e_{t-1}
    rer         % real exchange rate e_t*P*_t/P_t
    NFA_H       % home NFA = de * b_H     
    NFA_F       % foreign NFA = b_F / de

    %--- TRADE ---
    ex_H      % home exports (final good + intermediates)
    ex_F

    %--- POLICIES ---
    tau_H       % home carbon tax                         
    tau_F
    T_H         % lump-sum redistribution
    T_F
    ;


%----------------------------------------------------------------
% 2. Parameters
%----------------------------------------------------------------
parameters
    % Country size
        l_H l_F   % country size

    % Households
        % Preferences 
            sigmaC_H sigmaC_F   % risk aversion
            sigmaH_H sigmaH_F   % inverse Frisch elasticity
            beta                % discount factor
            hc_H hc_F           % habit persistence
            psi_H psi_F         % labor disutility weight

        % CES aggreator
            phi                 % elasticity home/foreign final consumption CES
            gamma_c_H           % import share in home consumption
            gamma_c_F           % import share in foreign consumption
        
        % other
            psi_B               % AC parameter for bond holdings


    % Final good production
        alpha_h             % final good TFP home
        alpha_f             % final good TFP foreign
        eta                 % elasticity across intermediate varieties
        zeta                % output elasticity of G in final good
        gamma_y_H           % share of foreign intermediates in home final sector
        gamma_y_F           % share of home intermediates in foreign final sector


    % Intermdiate good production
        %Production
        alpha               % labor share in intermediate prod
        Gamma_H             % TFP level home intermediate firms
        Gamma_F
    
        % Carbon / abatement
        sig_H sig_F         % emission intensity
        theta1_H theta2_F   % abatement cost home (theta_1 has no trend)
        theta1_H theta2_F
        tau0_H tau0_F       % carbon taxes

        %other
        kappa_H kappa_F     % price adjustment cost parameters (Rotemberg)
    
    % Trade
    tau_i               % iceberg trade cost

    %Government
        % Monetary policy
        rho                 % interest rate smoothing
        phi_pi              % reaction to inflation
        phi_y               % reaction to output 
    
    % Targets
    y_star      % Natural output (for Taylor rule)  
    pi_star     % Inflation target
    ;

%----------------------------------------------------------------
% 3. Calibration
%----------------------------------------------------------------
sigmaC_H  = 1.5;
sigmaC_F  = 1.2;
sigmaH_H  = 2.0;
sigmaH_F  = 1.9;
beta      = 0.994;
alpha     = 0.7;
hc_H      = 0.7;
hc_F      = 0.6;
psi_B     = 0.007; 
kappa_H   = 100;
kappa_F   = 80;
eta       = 10;
mu_ces    = 2;
zeta      = 0.7;
tau_i     = 0.05;
gamma_y_H = 0.3;
gamma_y_F = 0.3;
gamma_c_H = 0.1;
gamma_c_F = 0.11;
rho       = 0.8;
phi_pi    = 1.5;
phi_y     = 0.05;
n         = 0.4;
pi_ss     = 1.005;
gy_H      = 0.2;
gy_F      = 0.2;
tau0_H    = 0.050;  % home carbon tax ($/unit), higher
tau0_F    = 0.030;  % foreign carbon tax, lower -> CBAM wedge active
sig_H     = 0.2;
sig_F     = 0.2;
y0        = 25;
theta1    = 0.05;
theta2    = 2.6;


%----------------------------------------------------------------
% 5. Model Block
%----------------------------------------------------------------
model;

% ================================================================
% 1. HOUSEHOLDS
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

[name='CPI index, Home']
P_H = (1-gamma_c_H) * p_H^(1-phi) + gamma_c_H * rer^(1-phi);

[name='CPI index, Foreign']
P_F = (1-gamma_c_F) * p_F^(1-phi) + gamma_c_F * (1/rer)^(1-phi);

[name='Relative price dynamics, Home']
p_H / p_H(-1) = pi_H / cpi_H;

[name='Relative price dynamics, Foreign']
p_F / p_F(-1) = pi_F / cpi_F;


% ================================================================
% 2. FINAL GOOD SECTOR
% ================================================================

[name='Final good production, Home']
y_H = alpha_h * G_H^zeta * n^(1-zeta);

[name='Final good production, Foreign']
y_F = alpha_F * G_F^zeta * (1-n)^(1-zeta);

[name='Demand for intermediate aggregate, Home']
PG_H = zeta * p_H * y_H / G_H;

[name='Demand for intermediate aggregate, Foreign']
PG_F = zeta * p_F * y_F / G_F;

[name='Intermediate bundle price index, Home']
PG_H^(1-eta) = (1-gamma_y_H) * p_H^(1-eta) + gamma_y_H * (rer * ((1+tau_i)*p_F + (tau_H - tau_F)*sig_F*(1-mu_F)))^(1-eta);

[name='Intermediate bundle price index, Foreign']
PG_F^(1-eta) = (1-gamma_y_F) * p_F^(1-eta) + gamma_y_F * ((1+tau_i) * p_H / rer)^(1-eta);


% ================================================================
% 3. INTERMEDIATE GOOD SECTOR
% ================================================================

[name='Intermediate production, Home']
y_H = Gamma_H * l_H^(1-alpha) * n_H^alpha;

[name='Intermediate production, Foreign']
y_F = Gamma_F * l_F^(1-alpha) * n_F^alpha;

[name='Optimal abatement, Home']
theta1 * theta2 * mu_H^(theta2-1) = tau_H * sig_H;

[name='Marginal cost, Home']
mc_H = (1/alpha) * w_H * (n_H / y_H) + theta1 * mu_H^theta2 + tau_H * sig_H * (1 - mu_H);

[name='New Keynesian Phillips Curve, Home']
kappa_H * (pi_H - pi_star) * pi_H = (1 - eta) * (p_H / PG_H) + eta * mc_H + beta * kappa_H * (pi_H(+1) - pi_star) * pi_H(+1) * (lb_H(+1) / lb_H) * y_H(+1) / y_H;

[name='Optimal abatement, Foreign']
theta1 * theta2 * mu_F^(theta2-1) = tau_F * sig_F;

[name='Marginal cost, Foreign']
mc_F = (1/alpha) * w_F * (n_F / y_F) + theta1 * mu_F^theta2 + tau_F * sig_F * (1 - mu_F);

[name='New Keynesian Phillips Curve, Foreign']
kappa_F * (pi_F - pi_star) * pi_F = (1 - eta) * (p_F / PG_F) + eta * mc_F + beta * kappa_F * (pi_F(+1) - pi_star) * pi_F(+1) * (lb_F(+1) / lb_F) * y_F(+1) / y_F;


[name='Emissions, Home']
e_H = n * sig_H * (1 - mu_H) * y_H;

[name='Emissions, Foreign']
e_F = (1-n) * sig_F * (1 - mu_F) * y_F;


% ================================================================
% 4. PUBLIC SECTOR & MONETARY POLICY
% ================================================================

[name='Taylor rule, Home']
r_H = r_H(-1)^rho * (STEADY_STATE(r_H) * (cpi_H / STEADY_STATE(cpi_H))^phi_pi * (y_H / STEADY_STATE(y_H))^phi_y)^(1-rho);

[name='Taylor rule, Foreign']
r_F = r_F(-1)^rho * (STEADY_STATE(r_F) * (cpi_F / STEADY_STATE(cpi_F))^phi_pi * (y_F / STEADY_STATE(y_F))^phi_y)^(1-rho);

[name='Government spending, Home']
g_H = gy_H * STEADY_STATE(y_H);

[name='Government spending, Foreign']
g_F = gy_F * STEADY_STATE(y_F);

[name='Carbon tax, Home']
tau_H = tau0_H;

[name='Carbon tax, Foreign']
tau_F = tau0_F;


% ================================================================
% 5. AGGREGATION & MARKET CLEARING
% ================================================================

[name='Resource constraint, Home']
y_H = (1-gamma_y_H) * (p_H / PG_H)^(-eta) * G_H + gamma_y_F * ((1+tau_i) * p_H / rer / PG_F)^(-eta) * G_F + (1-gamma_c_H) * p_H^(-phi) * c_H + gamma_c_F * (p_H / rer)^(-phi) * c_F * (1-n)/n + g_H + theta1 * mu_H^theta2 * y_H + 0.5 * kappa_H * (pi_H - pi_star)^2 * y_H + 0.5 * psi_B * (de * b_H)^2 / (p_H * y_H);

[name='Resource constraint, Foreign']
y_F = (1-gamma_y_F) * (p_F / PG_F)^(-eta) * G_F + gamma_y_H * (rer * ((1+tau_i)*p_F + (tau_H - tau_F)*sig_F*(1-mu_F)) / PG_H)^(-eta) * G_H + (1-gamma_c_F) * p_F^(-phi) * c_F + gamma_c_H * (p_F * rer)^(-phi) * c_H * n/(1-n) + g_F + theta1 * mu_F^theta2 * y_F + 0.5 * kappa_F * (pi_F - pi_star)^2 * y_F + 0.5 * psi_B * (b_F / de)^2 / (p_F * y_F);

[name='NFA identity, Home']
NFA_H = de * b_H;

[name='NFA identity, Foreign']
NFA_F = b_F / de;

[name='NFA accumulation, Home']
NFA_H = r_F(-1) / cpi_H * NFA_H(-1) * de + p_H * gamma_y_F * ((1+tau_i)*p_H/rer/PG_F)^(-eta) * G_F - rer * ((1+tau_i)*p_F + (tau_H-tau_F)*sig_F*(1-mu_F)) * gamma_y_H * (rer*((1+tau_i)*p_F+(tau_H-tau_F)*sig_F*(1-mu_F))/PG_H)^(-eta) * G_H + p_H * gamma_c_F * (p_H/rer)^(-phi) * c_F * (1-n)/n - rer * gamma_c_H * p_H^(-phi) * c_H;

[name='Bond market clearing (Walras law)']
n * b_H + (1-n) * b_F = 0;

[name='UIP / FOC on foreign bonds, Home']
r_H = de(+1) * r_F + psi_B * de * NFA_H / (p_H * y_H) * cpi_H(+1);

[name='Real exchange rate']
rer / rer(-1) = de * cpi_F / cpi_H;

[name='Exports, Home']
ex_H = gamma_y_F * ((1+tau_i)*p_H/rer/PG_F)^(-eta) * G_F + gamma_c_F * (p_H/rer)^(-phi) * c_F * (1-n);

[name='Exports, Foreign']
ex_F = gamma_y_H * (rer*((1+tau_i)*p_F+(tau_H-tau_F)*sig_F*(1-mu_F))/PG_H)^(-eta) * G_H + gamma_c_H * (p_F*rer)^(-phi) * c_H * n;

end;

%----------------------------------------------------------------
% 6. Checks
%----------------------------------------------------------------
resid(1);
check;
