% Compare perfect foresight paths with and without nominal rigidities
% Usage (from MATLAB in dynare folder):
%   compare_pf_nominal_rigidities

clearvars;
close all;
clc;

base_mod = 'pf_cbam.mod';
flex_mod = 'pf_cbam_transition_flex.mod';

% 1) Baseline run (nominal rigidities as calibrated in the .mod file)
dynare(base_mod, 'noclearall', 'nolog');

path_sticky = oo_.endo_simul;
endo_names_sticky = cellstr(M_.endo_names);
T_sticky = size(path_sticky, 2);

% 2) Build a flexible-price counterfactual .mod with kappa_H = kappa_F = 0
txt = fileread(base_mod);
lines = splitlines(txt);
for k = 1:numel(lines)
    l = strtrim(lines{k});
    if ~isempty(regexp(l, '^kappa_H\s*=\s*', 'once'))
        indent = regexp(lines{k}, '^\s*', 'match', 'once');
        lines{k} = [indent 'kappa_H   = 0;'];
    elseif ~isempty(regexp(l, '^kappa_F\s*=\s*', 'once'))
        indent = regexp(lines{k}, '^\s*', 'match', 'once');
        lines{k} = [indent 'kappa_F   = 0;'];
    end
end
txt_flex = strjoin(lines, newline);
fid = fopen(flex_mod, 'w');
if fid < 0
    error('Could not create %s', flex_mod);
end
fwrite(fid, txt_flex);
fclose(fid);

% 3) Flexible-price run from scratch
dynare(flex_mod, 'noclearall', 'nolog');

path_flex = oo_.endo_simul;
endo_names_flex = cellstr(M_.endo_names);
T_flex = size(path_flex, 2);

% 4) Variables kept for diagnostics (table at the end)
vars_to_compare = {
    'tau_H','tau_F', ...
    'mu_H','mu_F', ...
    'e_H','e_F', ...
    'y_H','y_F', ...
    'c_H','c_F', ...
    'n_H','n_F', ...
    'w_H','w_F', ...
    'p_int_H','p_int_F', ...
    'PG_H','PG_F', ...
    'cpi_H','cpi_F', ...
    'rer', ...
    'ex_H','ex_F', ...
    'NFA_H','NFA_F', ...
    'mc_H','mc_F', ...
    'r_H','r_F', ...
    'pi_H','pi_F', ...
    'T_H','T_F'
};

% 5) Grouped layout plots (selected variables only)
T_plot = min(T_sticky, T_flex);
time = 1:T_plot;

plot_groups = {
    struct('title','Consumption and Labor','vars',{{'c_H','c_F','n_H','n_F','w_H','w_F'}},'nrows',2,'ncols',3), ...
    struct('title','Prices and Inflation','vars',{{'p_int_H','p_int_F','PG_H','PG_F','cpi_H','cpi_F','pi_H','pi_F'}},'nrows',2,'ncols',4), ...
    struct('title','Output and Emissions','vars',{{'y_H','y_F','e_H','e_F','mu_H','mu_F','mc_H','mc_F'}},'nrows',2,'ncols',4), ...
    struct('title','Policy and Rates','vars',{{'tau_H','tau_F','T_H','T_F','r_H','r_F'}},'nrows',2,'ncols',3), ...
    struct('title','External Block','vars',{{'rer','ex_H','ex_F','NFA_H','NFA_F'}},'nrows',2,'ncols',3)
};

for g = 1:numel(plot_groups)
    G = plot_groups{g};
    figure('Name', ['PF compare - ' G.title]);
    tl = tiledlayout(G.nrows, G.ncols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, G.title);

    for i = 1:numel(G.vars)
        v = G.vars{i};
        idx_s = find(strcmp(strtrim(endo_names_sticky), v), 1);
        idx_f = find(strcmp(strtrim(endo_names_flex), v), 1);
        if isempty(idx_s) || isempty(idx_f)
            continue;
        end

        nexttile;
        plot(time, path_sticky(idx_s, 1:T_plot), 'LineWidth', 1.6);
        hold on;
        plot(time, path_flex(idx_f, 1:T_plot), '--', 'LineWidth', 1.6);
        hold off;
        grid on;
        title(v, 'Interpreter', 'none');
        xlabel('Period');
    end

    lg = legend('Nominal rigidities (kappa > 0)', 'Flexible prices (kappa = 0)', 'Location', 'southoutside', 'Orientation', 'horizontal');
    lg.Layout.Tile = 'south';
end

% Optional quick check on final-period gaps
disp('Last-period difference (sticky - flexible):');
for i = 1:numel(vars_to_compare)
    v = vars_to_compare{i};
    idx_s = find(strcmp(strtrim(endo_names_sticky), v), 1);
    idx_f = find(strcmp(strtrim(endo_names_flex), v), 1);
    if ~isempty(idx_s) && ~isempty(idx_f)
        d = path_sticky(idx_s, T_plot) - path_flex(idx_f, T_plot);
        fprintf('%12s : %+ .6e\n', v, d);
    end
end
