% ================================================================
% CALL_FILE.M  —  Joint calibration (corrected)
% ================================================================

% 1. Initial Dynare run
dynare pf_cbam.mod

% World population reference (millions)
N_total = 1430;
pop_scale = N_total / 1000; % 1.43 billion

% ================================================================
% 2. Interactive target definition
% ================================================================

% --- Total GDP by region (trillion $) — ANNUAL values divided by 4 ---
% Model is quarterly so annual data is converted to quarterly.
% Trade balance (TB/GDP ratio): frequency-invariant, unchanged.
if exist('target_GDP_H_total','var') && ~isempty(target_GDP_H_total)
    t_GDP_H = target_GDP_H_total;
else
    ans_str = input('Enter EU quarterly GDP in T$ (default 4.25 = 17/4): ','s');
    t_GDP_H = parse_target(ans_str, 17/4, 'target_GDP_H_total');
end

if exist('target_GDP_F_total','var') && ~isempty(target_GDP_F_total)
    t_GDP_F = target_GDP_F_total;
else
    ans_str = input('Enter RoW quarterly GDP in T$ (default 13.5 = 54/4): ','s');
    t_GDP_F = parse_target(ans_str, 54/4, 'target_GDP_F_total');
end

% Convert to world per-capita targets
t_y_H = t_GDP_H / pop_scale;
t_y_F = t_GDP_F / pop_scale;

% --- Trade balance ---
if exist('target_tb_H','var') && ~isempty(target_tb_H)
    t_tb_H = target_tb_H;
else
    ans_str = input('Enter target tb_H (default 0.009, i.e. +0.9% of GDP): ','s');
    t_tb_H = parse_target(ans_str, 0.009, 'target_tb_H');
end

% --- Total emissions (billion tons) — ANNUAL values divided by 4 ---
% Same logic as GDP: convert annual to quarterly emissions.
target_emissions = true;
if exist('target_e_H','var') && ~isempty(target_e_H)
    t_e_H = target_e_H;
elseif exist('skip_emissions','var') && skip_emissions
    target_emissions = false;
    t_e_H = NaN;
else
    ans_str = input('Enter target e_H quarterly emissions in billion tons (default 0.625 = 2.5/4): ','s');
    val_e_H = parse_target(ans_str, 2.5/4, 'target_e_H');
    if val_e_H == 0
        target_emissions = false;
        t_e_H = NaN;
    else
        t_e_H = val_e_H / pop_scale; % Normalized to world per-capita
    end
end

if target_emissions
    if exist('target_e_F','var') && ~isempty(target_e_F)
        t_e_F = target_e_F;
    else
        ans_str = input('Enter target e_F quarterly emissions in billion tons (default 2.55 = 10.2/4): ','s');
        val_e_F = parse_target(ans_str, 10.2/4, 'target_e_F');
        if val_e_F == 0
            target_emissions = false;
            t_e_F = NaN;
        else
            t_e_F = val_e_F / pop_scale; % Normalized to world per-capita
        end
    end
else
    t_e_F = NaN;
end

% ================================================================
% 3. Indexes of parameters and variables
% ================================================================
idp_alpha_h   = strmatch('alpha_h',   M_.param_names, 'exact');
idp_alpha_F   = strmatch('alpha_F',   M_.param_names, 'exact');
idp_gamma_y_H = strmatch('gamma_y_H', M_.param_names, 'exact');
idp_gamma_y_F = strmatch('gamma_y_F', M_.param_names, 'exact');
idp_gamma_c_H = strmatch('gamma_c_H', M_.param_names, 'exact');
idp_gamma_c_F = strmatch('gamma_c_F', M_.param_names, 'exact');
idp_sig_H     = strmatch('sig_H',     M_.param_names, 'exact');
idp_sig_F     = strmatch('sig_F',     M_.param_names, 'exact');

idy_y_H  = strmatch('y_H',  M_.endo_names, 'exact');
idy_y_F  = strmatch('y_F',  M_.endo_names, 'exact');
idy_tb_H = strmatch('tb_H', M_.endo_names, 'exact');
idy_e_H  = strmatch('e_H',  M_.endo_names, 'exact');
idy_e_F  = strmatch('e_F',  M_.endo_names, 'exact');

gamma0 = [ M_.params(idp_gamma_y_H) ; M_.params(idp_gamma_y_F) ; M_.params(idp_gamma_c_H) ; M_.params(idp_gamma_c_F) ];
sig0   = [ M_.params(idp_sig_H) ; M_.params(idp_sig_F) ];

lb_gamma = [0.001 ; 0.001 ; 0.001 ; 0.001];
ub_gamma = [0.950 ; 0.950 ; 0.500 ; 0.500];

lb_sig = [0.001 ; 0.001];
ub_sig = [10.0  ; 10.0 ];

% ================================================================
% 4. Coordination loop parameters
% ================================================================
MAX_OUTER  = 40;
TOL_OUTER  = 1e-6;
RELAX      = 0.6;    

opts_fmc = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'TolFun', 1e-9, 'TolCon', 1e-8, 'TolX', 1e-9, 'MaxIterations', 500, 'MaxFunctionEvaluations', 5000);

params_work = M_.params;

fprintf('\n=== CALIBRATION CONJOINTE — 6 PARAMÈTRES (FRÉQUENCE TRIMESTRIELLE) ===\n');
fprintf('Poids pop : l_H = %.4f (EU 450M) | l_F = %.4f (RoW 980M) | N_total = %d M\n', ...
    M_.params(strmatch('l_H', M_.param_names,'exact')), ...
    M_.params(strmatch('l_F', M_.param_names,'exact')), N_total);
fprintf('\nCONVERSION CIBLES TRIMESTRIELLES → PER CAPITA MONDIAL (÷ %.2f bn) :\n', pop_scale);
fprintf('  EU PIB  (trim.) : %5.2f T$ total  →  y_H cible = %.4f\n', t_GDP_H, t_y_H);
fprintf('  RoW PIB (trim.) : %5.2f T$ total  →  y_F cible = %.4f\n', t_GDP_F, t_y_F);
if target_emissions
    fprintf('  EU Emi  (trim.) : %5.4f bn T total →  e_H cible = %.4f\n', val_e_H, t_e_H);
    fprintf('  RoW Emi (trim.) : %5.4f bn T total →  e_F cible = %.4f\n', val_e_F, t_e_F);
end

err_global = Inf;

for outer = 1:MAX_OUTER
    % STEP A: iterative scaling alpha_h / alpha_F
    for inner = 1:40
        try
            [ss_vec, ~, ss_check] = eval([M_.fname '.steadystate(oo_.steady_state, transpose(oo_.exo_steady_state), params_work)']);
        catch
            break;
        end
        if ss_check || any(~isfinite(ss_vec)), break; end

        y_H_cur = real(ss_vec(idy_y_H));
        y_F_cur = real(ss_vec(idy_y_F));

        if abs(y_H_cur - t_y_H)/t_y_H < 1e-9 && abs(y_F_cur - t_y_F)/t_y_F < 1e-9, break; end

        scale_h = min(max(t_y_H / max(abs(y_H_cur), 1e-12), 0.5), 2.0);
        scale_F = min(max(t_y_F / max(abs(y_F_cur), 1e-12), 0.5), 2.0);

        params_work(idp_alpha_h) = params_work(idp_alpha_h) * (scale_h ^ RELAX);
        params_work(idp_alpha_F) = params_work(idp_alpha_F) * (scale_F ^ RELAX);
    end

    % STEP B: fmincon on 4 gamma parameters
    g0 = [ params_work(idp_gamma_y_H) ; params_work(idp_gamma_y_F) ; params_work(idp_gamma_c_H) ; params_work(idp_gamma_c_F) ];
    W = diag([1, 1, 1, 1]);
    obj_gamma = @(g) (g - gamma0)' * W * (g - gamma0);
        N_total = 1430; % World population reference (millions)
    nonlcon_gamma = @(g) deal([], tb_constraint(g, idp_gamma_y_H, idp_gamma_y_F, idp_gamma_c_H, idp_gamma_c_F, t_tb_H, idy_tb_H, params_work, M_, oo_));

    [g_opt, ~, exitflag_b] = fmincon(obj_gamma, g0, [], [], [], [], lb_gamma, ub_gamma, nonlcon_gamma, opts_fmc);

    if exitflag_b > 0 || exitflag_b == -2
        params_work(idp_gamma_y_H) = g_opt(1);
        params_work(idp_gamma_y_F) = g_opt(2);
        params_work(idp_gamma_c_H) = g_opt(3);
        params_work(idp_gamma_c_F) = g_opt(4);
    end

    % STEP C: calibrate sig_H / sig_F
    if target_emissions
        for inner_sig = 1:30
            try
                [ss_sig, ~, chk_sig] = eval([M_.fname '.steadystate(oo_.steady_state, transpose(oo_.exo_steady_state), params_work)']);
            catch
                break;
            end
            if chk_sig || any(~isfinite(ss_sig)), break; end

            e_H_cur = real(ss_sig(idy_e_H));
            e_F_cur = real(ss_sig(idy_e_F));

        t_y_H = t_GDP_H / pop_scale; % Convert to world per-capita targets
        t_y_F = t_GDP_F / pop_scale;
            if err_eH < 1e-8 && err_eF < 1e-8, break; end

            sc_H = min(max(t_e_H / max(abs(e_H_cur), 1e-12), 0.5), 2.0);
            sc_F = min(max(t_e_F / max(abs(e_F_cur), 1e-12), 0.5), 2.0);
            params_work(idp_sig_H) = params_work(idp_sig_H) * (sc_H ^ RELAX);
            params_work(idp_sig_F) = params_work(idp_sig_F) * (sc_F ^ RELAX);
        end
    end

    % Global evaluation
    try
        [ss_vec, ~, ss_check] = eval([M_.fname '.steadystate(oo_.steady_state, transpose(oo_.exo_steady_state), params_work)']);
    catch ME
        continue;
    end
    if ss_check || any(~isfinite(ss_vec))
        continue;
    end

    err_yH = abs(real(ss_vec(idy_y_H)) - t_y_H) / max(abs(t_y_H), 1e-12);
    err_yF = abs(real(ss_vec(idy_y_F)) - t_y_F) / max(abs(t_y_F), 1e-12);
    err_tb = abs(real(ss_vec(idy_tb_H)) - t_tb_H) / max(abs(t_tb_H), 1e-12);

    if target_emissions
        err_eH = abs(real(ss_vec(idy_e_H)) - t_e_H) / max(abs(t_e_H), 1e-12);
        err_eF = abs(real(ss_vec(idy_e_F)) - t_e_F) / max(abs(t_e_F), 1e-12);
        err_global = max([err_yH, err_yF, err_tb, err_eH, err_eF]);
    else
        err_global = max([err_yH, err_yF, err_tb]);
    end

    if err_global < TOL_OUTER
        fprintf('\nConvergence atteinte en %d itérations (err=%.2e).\n', outer, err_global);
        break;
    end
end

if err_global >= TOL_OUTER
    warning('Boucle externe non convergée (err_global=%.2e).', err_global);
end

M_.params(idp_alpha_h)   = params_work(idp_alpha_h);
M_.params(idp_alpha_F)   = params_work(idp_alpha_F);
M_.params(idp_gamma_y_H) = params_work(idp_gamma_y_H);
M_.params(idp_gamma_y_F) = params_work(idp_gamma_y_F);
M_.params(idp_gamma_c_H) = params_work(idp_gamma_c_H);
M_.params(idp_gamma_c_F) = params_work(idp_gamma_c_F);
if target_emissions
    M_.params(idp_sig_H) = params_work(idp_sig_H);
    M_.params(idp_sig_F) = params_work(idp_sig_F);
end

[ys, M_.params, info] = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);

% ================================================================
% 6. Final report (corrected formulas + Dynare format)
% ================================================================
fprintf('\n╔══════════════════════════════════════════════════════════════════════╗\n');
fprintf('║              RÉSULTATS DE CALIBRATION  [TRIMESTRIEL]                 ║\n');
fprintf('╠══════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  SCALING POPULATION                                                  ║\n');
fprintf('║  EU 450M + RoW 980M = 1430M  |  N_total référence = %4d M          ║\n', N_total);
fprintf('║  l_H = %.4f (EU)   l_F = %.4f (RoW)                            ║\n', ...
    M_.params(strmatch('l_H', M_.param_names,'exact')), ...
    M_.params(strmatch('l_F', M_.param_names,'exact')));
fprintf('╠══════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  CIBLES PIB TRIMESTRIELLES (→ per capita mondiale ÷ %.2f bn)        ║\n', pop_scale);
fprintf('║  Annuel→Trim : EU 17T$/4=4.25T$ | RoW 54T$/4=13.5T$                ║\n');
fprintf('║  y_H : EU  %5.2f T$ trim. → cible=%7.4f | obtenu=%7.4f         ║\n', t_GDP_H, t_y_H, ys(idy_y_H));
fprintf('║  y_F : RoW %5.2f T$ trim. → cible=%7.4f | obtenu=%7.4f         ║\n', t_GDP_F, t_y_F, ys(idy_y_F));
fprintf('║  y_H_annual_check : %6.2f T$ (= obtenu × N_total/1000 × 4)       ║\n', ys(idy_y_H) * pop_scale * 4);
fprintf('║  y_F_annual_check : %6.2f T$ (= obtenu × N_total/1000 × 4)       ║\n', ys(idy_y_F) * pop_scale * 4);
fprintf('╠══════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  BALANCE COMMERCIALE (ratio TB/PIB — invariant à la fréquence)      ║\n');
fprintf('║  tb_H : cible = %+.6f  |  obtenu = %+.6f                    ║\n', t_tb_H, ys(idy_tb_H));
fprintf('╠══════════════════════════════════════════════════════════════════════╣\n');
if target_emissions
    fprintf('║  ÉMISSIONS TRIMESTRIELLES (÷ pop mondiale %.2f bn)                 ║\n', pop_scale);
    fprintf('║  Annuel→Trim : EU 2.5Gt/4=0.625Gt | RoW 10.2Gt/4=2.55Gt          ║\n');
    fprintf('║  e_H : cible = %9.4f  |  obtenu = %9.4f                    ║\n', t_e_H, ys(idy_e_H));
    fprintf('║  e_F : cible = %9.4f  |  obtenu = %9.4f                    ║\n', t_e_F, ys(idy_e_F));
    fprintf('║  e_H_annual_check : %5.2f bn T (= obtenu × N_total/1000 × 4)      ║\n', ys(idy_e_H) * pop_scale * 4);
    fprintf('║  e_F_annual_check : %5.2f bn T (= obtenu × N_total/1000 × 4)      ║\n', ys(idy_e_F) * pop_scale * 4);
end
        fprintf('\n=== JOINT CALIBRATION — 6 PARAMETERS (QUARTERLY) ===\n');

% --- COPY-PASTE BLOCK FOR PF_CBAM.MOD ---
fprintf('\n%% ================================================================\n');
        fprintf('\nQUARTERLY → WORLD PER-CAPITA TARGETS (÷ %.2f bn):\n', pop_scale);
fprintf('%% ================================================================\n\n');

fprintf('alpha_h   = %.6f;\n', M_.params(idp_alpha_h));
fprintf('alpha_F   = %.6f;\n', M_.params(idp_alpha_F));
fprintf('gamma_y_H = %.6f;\n', M_.params(idp_gamma_y_H));
fprintf('gamma_y_F = %.6f;\n', M_.params(idp_gamma_y_F));
fprintf('gamma_c_H = %.6f;\n', M_.params(idp_gamma_c_H));
fprintf('gamma_c_F = %.6f;\n', M_.params(idp_gamma_c_F));

if target_emissions
    fprintf('sig_H     = %.6f;\n', M_.params(idp_sig_H));
    fprintf('sig_F     = %.6f;\n', M_.params(idp_sig_F));
end
fprintf('\n%% ================================================================\n\n');

% End of original file (local functions unchanged below)

function ceq = tb_constraint(g, idp_gyH, idp_gyF, idp_gcH, idp_gcF, target_tb, idy_tb, params_in, M_in, oo_in)
    params_in(idp_gyH) = g(1);
    params_in(idp_gyF) = g(2);
    params_in(idp_gcH) = g(3);
    params_in(idp_gcF) = g(4);
    try
        [ss_vec, ~, ss_check] = eval([M_in.fname '.steadystate(oo_in.steady_state, transpose(oo_in.exo_steady_state), params_in)']);
        if ss_check || any(~isfinite(ss_vec))
            ceq = 1e6;
            return;
        end
        ceq = real(ss_vec(idy_tb)) - target_tb;
    catch
        ceq = 1e6;
    end
end

function val = parse_target(str, default_val, name)
    if isempty(str)
        val = default_val;
    else
        val = str2double(str);
        if isnan(val)
            error('Invalid numeric input for %s', name);
        end
    end
end