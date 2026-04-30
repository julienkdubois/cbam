function compare_pf(start_dates)
% compare_pf  Runs pf_cbam with different CBAM start dates or compares nominal rigidities.
% Usage: compare_pf([1, 13, 30]) % Compares multiple start dates using sticky prices
%        compare_pf()           % Compares sticky vs flexible prices for the default start date

global oo_ M_ options_;

clearvars -except start_dates oo_ M_ options_;
close all;
clc;

base_mod = 'pf_cbam.mod';
flex_mod = 'pf_cbam_transition_flex.mod';

% ================================================================
% 0.  LOAD ESTIMATED PARAMETERS (interactive)
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
    else
        fprintf('  ⚠ Aucun paramètre trouvé dans le mode file. Calibration de base utilisée.\n');
    end
else
    fprintf('  → Calibration de base conservée.\n');
end

if use_estimated
    patched_base = 'pf_cbam__estimated.mod';
    patch_mod_with_estimates(base_mod, patched_base, estimated_params, mode_file);
    run_mod = patched_base;
    fprintf('  → Variant patché : %s\n', patched_base);
else
    run_mod = base_mod;
end

% ================================================================
% Save figures option
% ================================================================
save_figs = false;
answer = input('Sauvegarder les graphiques ? (y/n) : ', 's');
if strcmpi(strtrim(answer), 'y')
    save_figs = true;
    save_dir = fullfile(pwd, 'pf_results_figs');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    fprintf('  → Les graphiques seront sauvegardés dans : %s\n', save_dir);
end

% ================================================================
% Diagnostics & Layout settings
% ================================================================

vars_to_compare = {
    'tau_H','tau_F', 'mu_H','mu_F', 'e_H','e_F', 'y_H','y_F', ...
    'c_H','c_F', 'n_H','n_F', 'w_H','w_F', 'p_int_H','p_int_F', ...
    'PG_H','PG_F', 'cpi_H','cpi_F', 'rer', 'ex_H','ex_F', ...
    'NFA_H','NFA_F', 'mc_H','mc_F', 'r_H','r_F', 'pi_H','pi_F', 'T_H','T_F'
};

plot_groups = {
    struct('title','Consumption and Labor',  'vars',{{'c_H','c_F','n_H','n_F','w_H','w_F'}},           'nrows',2,'ncols',3), ...
    struct('title','Prices and Inflation',   'vars',{{'p_int_H','p_int_F','PG_H','PG_F','cpi_H','cpi_F','pi_H','pi_F'}}, 'nrows',2,'ncols',4), ...
    struct('title','Output and Emissions',   'vars',{{'y_H','y_F','e_H','e_F','mu_H','mu_F','mc_H','mc_F'}},  'nrows',2,'ncols',4), ...
    struct('title','Policy and Rates',       'vars',{{'tau_H','tau_F','T_H','T_F','r_H','r_F'}},        'nrows',2,'ncols',3), ...
    struct('title','External Block',         'vars',{{'rer','ex_H','ex_F','NFA_H','NFA_F'}},            'nrows',2,'ncols',3),...
    struct('title','Intermediates',         'vars',{{'y_int_H','y_int_F'}},            'nrows',2,'ncols',1)...
};

if nargin < 1 || isempty(start_dates)
    % ================================================================
    % 1. Baseline logic: Sticky vs Flexible (Default Start Date)
    % ================================================================
    dynare(run_mod, 'noclearall', 'nolog');
    
    path_sticky       = oo_.endo_simul;
    endo_names_sticky = cellstr(M_.endo_names);
    T_sticky          = size(path_sticky, 2);

    txt   = fileread(run_mod);
    txt_lines = splitlines(txt);
    for k = 1:numel(txt_lines)
        l = strtrim(txt_lines{k});
        if ~isempty(regexp(l, '^kappa_H\s*=\s*', 'once'))
            indent   = regexp(txt_lines{k}, '^\s*', 'match', 'once');
            txt_lines{k} = [indent 'kappa_H   = 0;'];
        elseif ~isempty(regexp(l, '^kappa_F\s*=\s*', 'once'))
            indent   = regexp(txt_lines{k}, '^\s*', 'match', 'once');
            txt_lines{k} = [indent 'kappa_F   = 0;'];
        end
    end
    txt_flex = strjoin(txt_lines, newline);
    fid = fopen(flex_mod, 'w');
    if fid < 0, error('Could not create %s', flex_mod); end
    fwrite(fid, txt_flex);
    fclose(fid);

    dynare(flex_mod, 'noclearall', 'nolog');
    
    path_flex       = oo_.endo_simul;
    endo_names_flex = cellstr(M_.endo_names);
    T_flex          = size(path_flex, 2);

    cbam_start_period = 13;
    match = regexp(txt, 'periods\s+(\d+):\d+', 'tokens', 'once');
    if ~isempty(match)
        cbam_start_period = str2double(match{1});
    end

    T_plot = min([T_sticky, T_flex, 40]);
    time   = 1:T_plot;

    legend_sticky = 'Nominal rigidities (kappa > 0)';
    legend_flex   = 'Flexible prices (kappa = 0)';
    if use_estimated
        legend_sticky = [legend_sticky ' [estimated]'];
        legend_flex   = [legend_flex   ' [estimated]'];
    end

    for g = 1:numel(plot_groups)
        G  = plot_groups{g};
        fig_h = figure('Name', ['PF compare - ' G.title]);
        tl = tiledlayout(G.nrows, G.ncols, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, G.title);

        for i = 1:numel(G.vars)
            v     = G.vars{i};
            idx_s = find(strcmp(strtrim(endo_names_sticky), v), 1);
            idx_f = find(strcmp(strtrim(endo_names_flex),   v), 1);
            if isempty(idx_s) || isempty(idx_f), continue; end
            
            nexttile;
            plot(time, path_sticky(idx_s, 1:T_plot), 'LineWidth', 1.6);
            hold on;
            plot(time, path_flex(idx_f,   1:T_plot), '--', 'LineWidth', 1.6);
            hold off;
            grid on;
            title(v, 'Interpreter', 'none');
            xlabel('Period');
        end

        lg = legend(legend_sticky, legend_flex, 'Location', 'southoutside', 'Orientation', 'horizontal');
        lg.Layout.Tile = 'south';
        
        if save_figs
            fig_name = sprintf('pf_compare_%s.png', regexprep(G.title, '\s+', '_'));
            fig_path = fullfile(save_dir, fig_name);
            saveas(fig_h, fig_path);
            fprintf('  → Sauvegardé : %s\n', fig_name);
        end
    end

else
    % ================================================================
    % 2. Multiple Start Dates logic: Loop across provided dates
    % ================================================================
    results = cell(length(start_dates), 1);
    legends = cell(length(start_dates), 1);
    
    txt = fileread(run_mod);
    
    for i = 1:length(start_dates)
        sd = max(1, start_dates(i)); % Avoid 0, minimum period is 1
        
        new_txt = regexprep(txt, 'periods\s+\d+:\d+', sprintf('periods %d:1000', sd));
        
        temp_mod = sprintf('pf_cbam_temp_%d.mod', sd);
        fid = fopen(temp_mod, 'w');
        fwrite(fid, new_txt);
        fclose(fid);
        
        dynare(temp_mod, 'noclearall', 'nolog');
        
        results{i}.path = oo_.endo_simul;
        results{i}.endo_names = cellstr(M_.endo_names);
        legends{i} = sprintf('Start date: %d', sd);
        
        delete(temp_mod);
    end
    
    T_plot = 40;
    for i = 1:length(start_dates)
        T_plot = min(T_plot, size(results{i}.path, 2));
    end
    time = 1:T_plot;
    
    colors = lines(length(start_dates));
    title_suffix = '';
    if use_estimated
        title_suffix = ' [estimated parameters]';
    end
    
    for g = 1:numel(plot_groups)
        G = plot_groups{g};
        fig_h = figure('Name', ['PF compare dates - ' G.title]);
        tl = tiledlayout(G.nrows, G.ncols, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, [G.title title_suffix]);

        for j = 1:numel(G.vars)
            v = G.vars{j};
            
            if isempty(find(strcmp(strtrim(results{1}.endo_names), v), 1))
                continue;
            end
            
            nexttile;
            hold on;
            h_lines = gobjects(length(start_dates), 1);
            for i = 1:length(start_dates)
                idx = find(strcmp(strtrim(results{i}.endo_names), v), 1);
                lw = max(1, 1.5 + 1 * (length(start_dates) - i));
                h_lines(i) = plot(time, results{i}.path(idx, 1:T_plot), 'LineWidth', lw, 'Color', colors(i,:));
            end
            
            yl = ylim;
            for i = 1:length(start_dates)
                sd = max(1, start_dates(i));
                if sd <= T_plot
                    p_patch = patch([sd T_plot T_plot sd], [yl(1) yl(1) yl(2) yl(2)], colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.1);
                    uistack(p_patch, 'bottom');
                end
            end
            
            hold off;
            grid on;
            title(v, 'Interpreter', 'none');
            xlabel('Period');
        end
        
        lg = legend(h_lines, legends, 'Location', 'southoutside', 'Orientation', 'horizontal');
        lg.Layout.Tile = 'south';
        
        if save_figs
            fig_name = sprintf('pf_compare_dates_%s.png', regexprep(G.title, '\s+', '_'));
            fig_path = fullfile(save_dir, fig_name);
            saveas(fig_h, fig_path);
            fprintf('  → Sauvegardé : %s\n', fig_name);
        end
    end
end

end

% ================================================================
% LOCAL FUNCTIONS
% ================================================================

function [estimated_params, mode_file] = load_estimated_parameters()
estimated_params = struct();
mode_file        = '';

candidates = {
    fullfile(pwd,                       'cbam_estimation', 'Output', 'cbam_estimation_mode.mat')
    fullfile(pwd,            'dynare',  'cbam_estimation', 'Output', 'cbam_estimation_mode.mat')
    fullfile(fileparts(mfilename('fullpath')), 'cbam_estimation', 'Output', 'cbam_estimation_mode.mat')
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
if ischar(names),   names = cellstr(names); end
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

function patch_mod_with_estimates(base_mod, out_mod, estimated_params, mode_file)
txt = fileread(base_mod);
excluded = {'tau_h_ss','tau_f_ss','p_h_ss','p_f_ss'};
override_lines = {};
model_param_names = parse_param_names(txt);
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

if isempty(override_lines)
    writelines(txt, out_mod);
    return
end

header = sprintf('%% ---- Estimated param overrides from %s ----\n', mode_file);
block  = [header sprintf('%s\n', override_lines{:})];

anchor = 'model;';
idx    = strfind(txt, anchor);
if isempty(idx)
    new_txt = [txt newline block];
else
    new_txt = [txt(1:idx(1)-1) block txt(idx(1):end)];
end

fid = fopen(out_mod, 'w');
if fid < 0, error('Cannot create %s', out_mod); end
fwrite(fid, new_txt);
fclose(fid);
end

function names = parse_param_names(mod_text)
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