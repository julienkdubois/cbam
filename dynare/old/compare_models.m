clear; close all; clc;

% ================================================================
%  USER SETTINGS — edit this block only
% ================================================================

% --- Shocks ---
SHOCKS = { ...
    'eta_t_H', 'Carbon tax shock — Home (tau_H up)'; ...
    'eta_t_F', 'Carbon tax shock — Foreign (tau_F up)'; ...
    'eta_z_H', 'TFP shock — Home'; ...
    'eta_p_H', 'Cost-push shock — Home'; ...
    'eta_r_H', 'Monetary policy shock — Home'; ...
    'eta_e',   'UIP shock' ...
};

% --- Plot style ---
C_CBAM   = [0.13 0.41 0.82];    % blue  — with CBAM
C_NOCBAM = [0.84 0.18 0.13];    % red   — no CBAM
LW       = 2.0;
FONT     = 8.5;
IRF_T    = 20;    % must match irf= in stoch_simul

% ================================================================
%  STEP 1 — CALIBRATION
% ================================================================
print_sep('STEP 1 — CALIBRATION AUDIT');

P1 = read_params('model.mod');
P2 = read_params('model_no_cbam.mod');
all_pnames = union(fieldnames(P1), fieldnames(P2));

fprintf('\n%-22s  %12s  %12s  %s\n','Parameter','With CBAM','No CBAM','Note');
fprintf('%s\n', repmat('-',1,66));
n_diff = 0;
for k = 1:numel(all_pnames)
    nm = all_pnames{k};
    v1 = get_or_nan(P1,nm);  v2 = get_or_nan(P2,nm);
    if     isnan(v1),            note = '<-- only in No-CBAM'; n_diff=n_diff+1;
    elseif isnan(v2),            note = '<-- only in CBAM';    n_diff=n_diff+1;
    elseif abs(v1-v2)>1e-12,    note = '*** DIFFER ***';      n_diff=n_diff+1;
    else,                        note = '';
    end
    if ~isempty(note)
        fprintf('  %-22s  %12.6f  %12.6f  %s\n', nm, v1, v2, note);
    end
end
if n_diff==0, fprintf('  >> All parameters identical.\n'); end

fprintf('\n  => Using calibration directly from .mod files.\n\n');

% ================================================================
%  STEP 2 — RUN DYNARE
% ================================================================
print_sep('STEP 2 — RUNNING DYNARE');

fprintf('\n  [1/2] Running model.mod (WITH CBAM)...\n\n');
dynare model noclearall nograph nointeractive;
[SS1, IRF1, NAMES1, M1_] = extract_dynare_results();

fprintf('\n  [2/2] Running model_no_cbam.mod (NO CBAM)...\n\n');
dynare model_no_cbam noclearall nograph nointeractive;
[SS2, IRF2, NAMES2, M2_] = extract_dynare_results();

% ================================================================
%  STEP 3 — IRF FIGURES  (one PDF per shock)
% ================================================================
print_sep('STEP 3 — IRF FIGURES');
fprintf('\n');

for s = 1:size(SHOCKS,1)
    shock  = SHOCKS{s,1};
    slabel = SHOCKS{s,2};

    % ---- Build IRF matrix for this shock ----
    % Each call returns a struct with all the series we want to plot,
    % already scaled to % deviation from SS (x100).
    D1 = build_plot_data(IRF1, SS1, NAMES1, shock, IRF_T);
    D2 = build_plot_data(IRF2, SS2, NAMES2, shock, IRF_T);
    
    % Alias the active shock variable so it maps to 'e_sh' for plotting
    shk_var = regexprep(shock, '^eta_', 'e_');
    if isfield(D1, shk_var), D1.e_sh = D1.(shk_var); end
    if isfield(D2, shk_var), D2.e_sh = D2.(shk_var); end

    if isempty(fieldnames(D1)) && isempty(fieldnames(D2))
        fprintf('  [skip] No IRF data found for shock: %s\n', shock);
        continue;
    end

    % ---- Layout and plotting ----

    fig = figure('Name', slabel, 'Color','w', 'Position',[30 30 1380 900]);

    % Title
    annotation(fig,'textbox',[0 0.955 1 0.042], ...
        'String', ['IRF — ' slabel], ...
        'FontSize',14,'FontWeight','bold', ...
        'EdgeColor','none','HorizontalAlignment','center');

    % --- Define panels as {field, label, row, col} in a 6x5 grid ---
    panels = { ...
        % Row 1
        'e_sh',      'Shock (e_.)',      1, 1; ...
        'rer',       'Real exch. rate',  1, 2; ...
        'NFA_H',     'NFA Home',         1, 3; ...

        % Row 2
        'pi_H',      'Inflation H',      2, 1; ...
        'pi_F',      'Inflation F',      2, 2; ...
        'r_H',       'Interest rate H',  2, 3; ...
        'r_F',       'Interest rate F',  2, 4; ...
        
        % Row 3
        'y_H',       'Output H',         3, 1; ...
        'y_int_H',   'Int. Output H',    3, 2; ...
    
        'y_F',       'Output F',         3, 3; ...
        'y_int_F',   'Int. Output F',    3, 4; ...
        
          
        % Row 4
        'c_H',       'Consump. H',       4, 1; ...
        
        'c_H_h',     'Cons. H (Home)',   4, 2; ...
        'c_H_f',     'Cons. H (Foreign)',4, 3; ...
        
        % Row 5
        'c_F',       'Consump. F',       5, 1; ...

        'c_F_f',     'Cons. F (Foreign)',5, 2; ...
        'c_F_h',     'Cons. F (Home)',   5, 3; ...

        'c_tot',     'Total Consump.',   5, 4; ...
        
        % Row 6
        'mu_H',      'Abatement H',      6, 1; ...
        'mu_F',      'Abatement F',      6, 2; ...
        'e_H',       'Emissions H',      6, 3; ...
        'e_F',       'Emissions F',      6, 4; ...
        'e_tot',     'Total Emissions',  6, 5; ...
    };

    NROWS = 6; NCOLS = 5;
    mg_l=0.055; mg_r=0.02; mg_b=0.07; mg_t=0.09;
    gap_w = 0.012; gap_h = 0.045;
    pw = (1-mg_l-mg_r - (NCOLS-1)*gap_w)/NCOLS;
    ph = (1-mg_b-mg_t - (NROWS-1)*gap_h)/NROWS;

    legend_done = false;
    for p = 1:size(panels,1)
        field = panels{p,1};
        lbl   = panels{p,2};
        pr    = panels{p,3};   % row (1=top)
        pc    = panels{p,4};   % col

        % Position: [left, bottom, width, height]
        left   = mg_l + (pc-1)*(pw+gap_w);
        bottom = 1 - mg_t - pr*ph - (pr-1)*gap_h;
        ax = axes(fig,'Position',[max(0,left) max(0,bottom) pw ph]); %#ok<LAXES>
        hold(ax,'on');

        h1=[]; h2=[];
        if isfield(D1, field) && ~all(isnan(D1.(field)))
            h1 = plot(ax, 0:IRF_T-1, D1.(field), ...
                      'Color',C_CBAM,'LineWidth',LW,'DisplayName','With CBAM');
        end
        if isfield(D2, field) && ~all(isnan(D2.(field)))
            h2 = plot(ax, 0:IRF_T-1, D2.(field), ...
                      'Color',C_NOCBAM,'LineWidth',LW,'LineStyle','--', ...
                      'DisplayName','No CBAM');
        end

        yline(ax,0,'Color',[0.55 0.55 0.55],'LineWidth',0.7, ...
              'HandleVisibility','off');

        title(ax,lbl,'FontSize',FONT,'FontWeight','bold');
        if pr==NROWS, xlabel(ax,'Time','FontSize',FONT-1); end
        ylabel(ax,'% dev. SS','FontSize',FONT-1);
        xlim(ax,[0 IRF_T-1]);
        grid(ax,'on'); box(ax,'on');
        set(ax,'GridAlpha',0.2,'TickDir','out','FontSize',FONT-1, ...
               'XTick',0:4:IRF_T-1);

        % One legend, top-left panel only
        if ~legend_done && (~isempty(h1)||~isempty(h2))
            hh = [h1,h2]; hh = hh(~isempty(hh));
            if numel(hh)>=1
                lg = legend(ax,'Location','best','FontSize',FONT-1,'Box','off');
                legend_done = true;
            end
        end
    end

    % Row labels (annotations on left margin)
    row_titles = {'Macro H','Macro F','Cons/RER','Rates & Prices','Trade','Emissions'};
    for r = 1:NROWS
        bottom = 1 - mg_t - r*ph - (r-1)*gap_h;
        annotation(fig,'textbox', ...
            [0, max(0, bottom), max(0, mg_l-0.005), min(1, ph)], ...
            'String', row_titles{r}, ...
            'FontSize',7.5,'FontWeight','bold','Color',[0.4 0.4 0.4], ...
            'Rotation',90,'EdgeColor','none', ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end

    fname = sprintf('IRF_%s.pdf', shock);
    exportgraphics(fig, fname, 'ContentType','vector');
    fprintf('  Saved: %s\n', fname);
end

print_sep('DONE');
fprintf('\n  Figures saved to: %s\n\n', pwd);


% ================================================================
%  FUNCTIONS
% ================================================================

function D = build_plot_data(IRF, SS, NAMES, shock, T)
    D = struct();
    nm = @(s) strtrim(s);
    idx = @(s) find(strcmp(strtrim(NAMES), s),1);


    function y = get_irf(vn)
        key = [vn '_' shock];
        if isfield(IRF, key)
            raw = IRF.(key)(:);
            ss_val = SS(idx(vn));
            if abs(ss_val) > 1e-10
                y = raw / ss_val * 100;   % % deviation from SS
            else
                y = raw * 100;            % level deviation (SS≈0)
            end
            if numel(y) < T, y(end+1:T) = 0; end
            y = y(1:T);
        else
            y = nan(T,1);
        end
    end

    std_vars = {'y_H','y_F','y_int_H','y_int_F','c_H','c_F','pi_H','pi_F', ...
                'r_H','r_F','rer','e_H','e_F','mu_H','mu_F', ...
                'c_H_h','c_H_f','c_F_h','c_F_f','NFA_H', ...
                'e_z_H', 'e_z_F', 'e_p_H', 'e_p_F', 'e_r_H', 'e_r_F', ...
                'e_x_H', 'e_x_F', 'e_t_H', 'e_t_F', 'e_e'};
    for k = 1:numel(std_vars)
        D.(std_vars{k}) = get_irf(std_vars{k});
    end

    % --- Aggregate variables ---
    if isfield(D, 'e_H') && isfield(D, 'e_F')
        D.e_tot = D.e_H + D.e_F;
    else
        D.e_tot = nan(T,1);
    end
    if isfield(D, 'c_H') && isfield(D, 'c_F')
        D.c_tot = D.c_H + D.c_F;
    else
        D.c_tot = nan(T,1);
    end
end

function [SS, IRF, NAMES, M_out] = extract_dynare_results()
% Extracts the results from global variables after running Dynare
    global M_ oo_
    SS    = oo_.steady_state;
    IRF   = oo_.irfs;
    NAMES = M_.endo_names;
    M_out = M_;
end

function P = read_params(modfile)
    txt = fileread(modfile);
    txt = regexprep(txt,'%[^\n]*','');
    txt = regexprep(txt,'model\b.*?\bend\b','','ignorecase');
    pat = '(?m)^\s*([A-Za-z_]\w*)\s*=\s*([+-]?\s*\d+\.?\d*(?:[eE][+-]?\d+)?)\s*;';
    toks = regexp(txt, pat, 'tokens');
    P = struct();
    for k = 1:numel(toks)
        nm = matlab.lang.makeValidName(strtrim(toks{k}{1}));
        P.(nm) = str2double(strrep(toks{k}{2},' ',''));
    end
end

function v = get_or_nan(S,f)
    if isfield(S,f), v=S.(f); else, v=nan; end
end
function v = get_or_nan_vec(vec,idx)
    if isempty(idx)||isnan(idx), v=nan; else, v=vec(idx); end
end
function print_sep(t)
    n=68; bar=repmat('=',1,n);
    pad=repmat(' ',1,max(0,floor((n-2-numel(t))/2)));
    fprintf('\n%s\n%s %s %s\n%s\n',bar,pad,t,pad,bar);
end
