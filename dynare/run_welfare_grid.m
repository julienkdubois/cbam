% ================================================================
% run_welfare_grid.m
% ----------------------------------------------------------------
% Driver script: welfare of country F as a function of tau_F,
% comparing the no-CBAM baseline to a world with CBAM.
%
% Welfare metric
%   W_F = sum_{t=1}^{T} beta^(t-1) * U(c_F_t, n_F_t)
%   U   = log(c_F) - psi_F * n_F^(1+sigmaH_F) / (1+sigmaH_F)
% ================================================================

% clear; close all; clc;

% ================================================================
% 0.  CONFIGURATION
% ================================================================

tau_F_grid  = linspace(0.00, 0.15, 31);
tau_H_ss    = 0.1;
T           = 200;
beta        = 0.994;
sigmaH_F    = 2.0;
mod_name    = 'welfare_cbam';

% ================================================================
% 0b.  LOAD ESTIMATED PARAMETERS (interactive)
% ================================================================

estimated_params = struct();
use_estimated    = false;

answer = input('Utiliser les paramètres estimés de cbam_estimation ? (y/n) : ', 's');
if strcmpi(strtrim(answer), 'y')
    [estimated_params, mode_file] = load_estimated_parameters();
    if ~isempty(fieldnames(estimated_params))
        use_estimated = true;
        fprintf('  → %d paramètre(s) chargé(s) depuis :\n    %s\n', ...
                numel(fieldnames(estimated_params)), mode_file);
        disp('  → Paramètres récupérés :');
        fn = fieldnames(estimated_params);
        for k = 1:numel(fn)
            fprintf('      %-20s = %.6g\n', fn{k}, estimated_params.(fn{k}));
        end

        % Override beta / sigmaH_F if estimated
        if isfield(estimated_params, 'beta')
            beta     = estimated_params.beta;
            fprintf('  → beta mis à jour : %.6g\n', beta);
        end
        if isfield(estimated_params, 'sigmah_f')
            sigmaH_F = estimated_params.sigmah_f;
            fprintf('  → sigmaH_F mis à jour : %.6g\n', sigmaH_F);
        end
    else
        fprintf('  ⚠ Aucun paramètre trouvé dans le mode file. Calibration de base utilisée.\n');
    end
else
    fprintf('  → Calibration de base conservée.\n');
end

% Build (or select) the Dynare .mod to use
if use_estimated
    base_mod    = [mod_name '.mod'];
    mod_name    = [mod_name '__estimated'];
    variant_mod = [mod_name '.mod'];
    n_overrides = build_estimated_variant_mod(base_mod, variant_mod, estimated_params, mode_file);
    fprintf('  → %d override(s) injecté(s) dans %s\n', n_overrides, variant_mod);
end

% ================================================================
% 1.  PRE-ALLOCATE RESULTS
% ================================================================

n_grid       = numel(tau_F_grid);
W_no_cbam    = nan(n_grid, 1);
W_cbam       = nan(n_grid, 1);
lambda_CEV   = nan(n_grid, 1);
disc_weights = beta .^ (0:T-1)';

% ================================================================
% 2.  MAIN LOOP
% ================================================================

fprintf('\n=== Welfare grid: %d points, cbam = tau_H_ss - tau_F_val (dynamic) ===\n', n_grid);
if use_estimated
    fprintf('    (paramètres estimés actifs)\n');
end

for i = 1:n_grid
    tau_F_val  = tau_F_grid(i);
    cbam_level = max(tau_H_ss - tau_F_val, 0);   % CBAM = tau_H - tau_F (>= 0)
    fprintf('[%2d/%2d] tau_F = %.4f  cbam = %.4f ... ', i, n_grid, tau_F_val, cbam_level);

    % ---- steady-state without CBAM --------------------------------
    try
        sig_H=0.2; sig_F=0.2; theta1_H=0.05; theta2_H=2.6;
        theta1_F=0.05; theta2_F=2.6; eta=2.0; gamma_y_H=0.3; gamma_y_F=0.3;
        tau_i=0.05; l_H=0.4; l_F=0.6; Gamma_H=1.0; Gamma_F=1.0;
        alpha=0.7; zeta=0.33; alpha_h=1.0; alpha_F=1.0;
        p_H_ss=1; p_F_ss=1; sigmaC_H=1.5; sigmaC_F=1.5; sigmaH_H=2.0;
        phi=2.0; gamma_c_H=0.1; gamma_c_F=0.1; pi_star=1;

        % Apply estimated structural params to the SS call when available
        if use_estimated
            if isfield(estimated_params,'eta'),        eta        = estimated_params.eta;        end
            if isfield(estimated_params,'phi'),        phi        = estimated_params.phi;        end
            if isfield(estimated_params,'sigmac_h'),   sigmaC_H   = estimated_params.sigmac_h;   end
            if isfield(estimated_params,'sigmac_f'),   sigmaC_F   = estimated_params.sigmac_f;   end
            if isfield(estimated_params,'sigmah_h'),   sigmaH_H   = estimated_params.sigmah_h;   end
            if isfield(estimated_params,'gamma_c_h'),  gamma_c_H  = estimated_params.gamma_c_h;  end
            if isfield(estimated_params,'gamma_c_f'),  gamma_c_F  = estimated_params.gamma_c_f;  end
            if isfield(estimated_params,'gamma_y_h'),  gamma_y_H  = estimated_params.gamma_y_h;  end
            if isfield(estimated_params,'gamma_y_f'),  gamma_y_F  = estimated_params.gamma_y_f;  end
            if isfield(estimated_params,'sig_h'),      sig_H      = estimated_params.sig_h;      end
            if isfield(estimated_params,'sig_f'),      sig_F      = estimated_params.sig_f;      end
            if isfield(estimated_params,'theta1_h'),   theta1_H   = estimated_params.theta1_h;   end
            if isfield(estimated_params,'theta1_f'),   theta1_F   = estimated_params.theta1_f;   end
            if isfield(estimated_params,'theta2_h'),   theta2_H   = estimated_params.theta2_h;   end
            if isfield(estimated_params,'theta2_f'),   theta2_F   = estimated_params.theta2_f;   end
            if isfield(estimated_params,'alpha'),      alpha      = estimated_params.alpha;      end
            if isfield(estimated_params,'zeta'),       zeta       = estimated_params.zeta;       end
        end

        [~,~,~,~,~,~,~, nH_i, nF_i, cH_i, cF_i, wH_i, wF_i] = ...
            ss_pf(0, tau_H_ss, tau_F_val, sig_F, sig_H, ...
                  theta1_H, theta2_H, theta1_F, theta2_F, eta, ...
                  gamma_y_H, gamma_y_F, tau_i, l_H, l_F, ...
                  Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, ...
                  p_H_ss, p_F_ss, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, ...
                  phi, gamma_c_H, gamma_c_F, pi_star, beta, 1, 1);
        psi_H_cal = (cH_i^(-sigmaC_H) * wH_i) / (nH_i^sigmaH_H);
        psi_F_cal = (cF_i^(-sigmaC_F) * wF_i) / (nF_i^sigmaH_F);

        [~,~,~,~,~,~,~, ~, nF_0, ~, cF_0, ~, ~] = ...
            ss_pf(0, tau_H_ss, tau_F_val, sig_F, sig_H, ...
                  theta1_H, theta2_H, theta1_F, theta2_F, eta, ...
                  gamma_y_H, gamma_y_F, tau_i, l_H, l_F, ...
                  Gamma_H, Gamma_F, alpha, zeta, alpha_h, alpha_F, ...
                  p_H_ss, p_F_ss, sigmaC_H, sigmaC_F, sigmaH_H, sigmaH_F, ...
                  phi, gamma_c_H, gamma_c_F, pi_star, beta, psi_H_cal, psi_F_cal);

        U_ss_0       = log(cF_0) - psi_F_cal * nF_0^(1+sigmaH_F) / (1+sigmaH_F);
        W_no_cbam(i) = U_ss_0 / (1 - beta);
    catch ME
        fprintf('SS FAILED (%s)\n', ME.message);
        continue
    end

    % ---- PF simulation with CBAM ----------------------------------
    try
        tau_F_int = round(tau_F_val  * 1e8);
        cbam_int  = round(cbam_level * 1e8);

        dynare_cmd = sprintf('dynare %s -Dtau_F_val=%d -Dcbam_val=%d noclearall nolog nostrict', ...
                             mod_name, tau_F_int, cbam_int);
        clear M_ oo_ options_ var_list_;
        evalc(dynare_cmd);

        idx_cF = find(strcmp(M_.endo_names, 'c_F'));
        idx_nF = find(strcmp(M_.endo_names, 'n_F'));
        if isempty(idx_cF) || isempty(idx_nF)
            error('Variable c_F or n_F not found in M_.endo_names');
        end

        cF_path = oo_.endo_simul(idx_cF, 2:T+1)';
        nF_path = oo_.endo_simul(idx_nF, 2:T+1)';
        U_path  = log(cF_path) - psi_F_cal .* nF_path.^(1+sigmaH_F) ./ (1+sigmaH_F);

        W_cbam(i)     = disc_weights' * U_path;
        lambda_CEV(i) = exp((W_no_cbam(i) - W_cbam(i)) * (1 - beta));

        fprintf('PF ok | W_cbam=%.4f\n', W_cbam(i));
    catch ME
        fprintf('PF FAILED (%s)\n', ME.message);
    end
end

% ================================================================
% 3.  SAVE RESULTS
% ================================================================

results.tau_F_grid    = tau_F_grid;
results.W_no_cbam     = W_no_cbam;
results.W_cbam        = W_cbam;
results.lambda_CEV    = lambda_CEV;
results.beta          = beta;
results.tau_H_ss      = tau_H_ss;
results.cbam_grid     = max(tau_H_ss - tau_F_grid, 0);  % dynamic cbam per grid point
results.use_estimated = use_estimated;
if use_estimated
    results.estimated_params = estimated_params;
    results.mode_file        = mode_file;
end

save('welfare_cbam_results.mat', 'results');
fprintf('\nResults saved to welfare_cbam_results.mat\n');

% ================================================================
% 4.  PLOT
% ================================================================

fig1 = figure('Name','Welfare CBAM','NumberTitle','off','Color','w');
plot(tau_F_grid, W_cbam, 'r-o', 'LineWidth',1.8, 'MarkerSize',4);
xlabel('\tau_F','FontSize',12);
ylabel('W_{CBAM}','FontSize',12);
ttl = 'Welfare under CBAM as a function of \tau_F';
if use_estimated, ttl = [ttl ' [estimated params]']; end
title(ttl,'FontSize',13);
grid on; box on;

% ================================================================
% LOCAL FUNCTIONS
% ================================================================

function [estimated_params, mode_file] = load_estimated_parameters()
% Search for cbam_estimation mode file in standard locations.
estimated_params = struct();
mode_file        = '';

candidates = {
    fullfile(pwd,                               'cbam_estimation', 'Output', 'cbam_estimation_mode.mat')
    fullfile(pwd,                    'dynare',  'cbam_estimation', 'Output', 'cbam_estimation_mode.mat')
    fullfile(fileparts(mfilename('fullpath')),  'cbam_estimation', 'Output', 'cbam_estimation_mode.mat')
};

for i = 1:numel(candidates)
    if exist(candidates{i}, 'file')
        mode_file = candidates{i};
        break
    end
end

if isempty(mode_file)
    fprintf('  ⚠ Mode file introuvable. Chemins testés :\n');
    for i = 1:numel(candidates)
        fprintf('      %s\n', candidates{i});
    end
    return
end

S = load(mode_file);
if ~isfield(S, 'xparam1') || ~isfield(S, 'parameter_names')
    fprintf('  ⚠ Mode file incomplet (xparam1 ou parameter_names manquant).\n');
    return
end

xparam1 = S.xparam1(:);
names   = S.parameter_names;
if ischar(names),    names = cellstr(names); end
if isnumeric(names), names = cellstr(names); end

n = min(numel(xparam1), numel(names));
for k = 1:n
    raw = names{k};
    if isstring(raw), raw = char(raw); end
    key = lower(strtrim(raw));
    key = strrep(key, ' ', '');
    key = strrep(key, '-', '_');
    key = strrep(key, '.', '_');
    if ~isempty(key) && isvarname(key)
        estimated_params.(key) = xparam1(k);
    end
end
end

% ----------------------------------------------------------------
function n_overrides = build_estimated_variant_mod(base_mod, variant_mod, estimated_params, mode_file)
% Inject estimated structural parameters into a copy of the .mod file.

txt = fileread(base_mod);
model_param_names = parse_model_parameter_names(txt);

excluded = {'tau_h_ss','tau_f_ss','p_h_ss','p_f_ss'};
override_lines = {};

for i = 1:numel(model_param_names)
    name = model_param_names{i};
    key  = lower(name);
    if ismember(key, excluded), continue; end
    if isfield(estimated_params, key)
        v = estimated_params.(key);
        if isfinite(v) && isscalar(v) && v > 0
            override_lines{end+1} = sprintf('%s = %.16g;', name, v); %#ok<AGROW>
        end
    end
end

n_overrides = numel(override_lines);
if n_overrides == 0
    writelines(txt, variant_mod);
    return
end

header = sprintf('%% ---- Estimated parameter overrides from %s ----\n', mode_file);
override_block = [header sprintf('%s\n', override_lines{:})];

anchor = '% ---- tau_F_ss is set to the grid value passed by the driver ----';
idx    = strfind(txt, anchor);
if ~isempty(idx)
    new_txt = [txt(1:idx(1)-1) override_block txt(idx(1):end)];
else
    model_anchor = 'model;';
    idx2 = strfind(txt, model_anchor);
    if isempty(idx2)
        new_txt = [txt newline override_block];
    else
        new_txt = [txt(1:idx2(1)-1) override_block txt(idx2(1):end)];
    end
end

fid = fopen(variant_mod, 'w');
if fid < 0, error('Cannot create variant mod file: %s', variant_mod); end
fwrite(fid, new_txt);
fclose(fid);
end

% ----------------------------------------------------------------
function names = parse_model_parameter_names(mod_text)
names     = {};
start_idx = strfind(mod_text, 'parameters');
if isempty(start_idx), return; end

tail    = mod_text(start_idx(1):end);
end_idx = strfind(tail, ';');
if isempty(end_idx), return; end

block  = tail(1:end_idx(1));
block  = regexprep(block, 'parameters', '', 'once');
block  = strrep(block, ';', ' ');
tokens = regexp(block, '[A-Za-z][A-Za-z0-9_]*', 'match');

if isempty(tokens), return; end

seen = containers.Map('KeyType','char','ValueType','logical');
for i = 1:numel(tokens)
    t = tokens{i};
    k = lower(t);
    if ~isKey(seen, k)
        seen(k) = true;
        names{end+1} = t; %#ok<AGROW>
    end
end
end