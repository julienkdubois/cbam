% 1. Exécution initiale de Dynare
dynare cbam_estimation.mod 

% 2. Définition interactive des cibles
% Use existing workspace variables if present, else prompt the user with sensible defaults
if exist('target_y_H','var') && ~isempty(target_y_H)
    t_y_H = target_y_H;
else
    ans_str = input('Enter target y_H (default 15): ','s');
    if isempty(ans_str)
        t_y_H = 15;
    else
        t_y_H = str2double(ans_str);
        if isnan(t_y_H)
            error('Invalid numeric input for target_y_H');
        end
    end
end
if exist('target_y_F','var') && ~isempty(target_y_F)
    t_y_F = target_y_F;
else
    ans_str = input('Enter target y_F (default 80): ','s');
    if isempty(ans_str)
        t_y_F = 80;
    else
        t_y_F = str2double(ans_str);
        if isnan(t_y_F)
            error('Invalid numeric input for target_y_F');
        end
    end
end
my_targets_ss = [t_y_H; t_y_F];

% 3. Paramètres de contrôle: 2 paramètres pour 2 cibles
idp_control = [strmatch('alpha_h', M_.param_names, 'exact'), ...
               strmatch('alpha_F', M_.param_names, 'exact')];

% 4. Indices des cibles
idy_target = [strmatch('y_H', M_.endo_names, 'exact'), ...
              strmatch('y_F', M_.endo_names, 'exact')];
idy_emissions = [strmatch('e_H', M_.endo_names, 'exact'), ...
                 strmatch('e_F', M_.endo_names, 'exact')];

% 5. Guesses initiaux réalistes autour de la calibration de départ
x0 = [1.0; 1.0];

% 6. Calibration directe par mise à l'échelle itérative (alpha_h, alpha_F)
maxiter = 50;
tol_abs = 1e-6;
relax = 0.7; % under-relaxation to stabilise
params_vec = M_.params; % working copy
err = Inf(2,1);
for it = 1:maxiter
    params_vec(idp_control) = params_vec(idp_control); % current candidates
    try
        [ss_vec, params_vec, ss_check] = eval([ M_.fname '.steadystate(oo_.steady_state, transpose(oo_.exo_steady_state), params_vec)']);
    catch
        warning('steady-state call failed at iter %d', it);
        break;
    end
    if ss_check || any(isnan(ss_vec))
        warning('steady-state not solved (check=%d) at iter %d', ss_check, it);
        break;
    end
    y_vals = real(ss_vec(idy_target));
    err = my_targets_ss(:) - y_vals;
    if all(abs(err) < tol_abs)
        fprintf('Scaling converged in %d iterations.\n', it);
        break;
    end
    % Compute scale factors and damp them
    y_safe = y_vals;
    y_safe(abs(y_safe) < 1e-12) = 1e-12;
    scale = (my_targets_ss(:) ./ y_safe);
    % constrain scale to avoid wild jumps
    scale = min(max(scale, 0.1), 10);
    params_vec(idp_control) = params_vec(idp_control) .* (scale .^ relax);
end
% Final injection
x_opt = params_vec(idp_control);
resnorm = sum(err.^2);

% 7. Injection et vérification finale via Dynare
M_.params(idp_control) = x_opt;
[ys, M_.params, info] = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);

fprintf('\n--- VERIFICATION DES CIBLES ---\n');
target_names = {'y_H', 'y_F'};
for i=1:2
    fprintf('%s: Cible = %.2f | Obtenu = %.4f\n', target_names{i}, my_targets_ss(i), ys(idy_target(i)));
end

fprintf('\n--- EMISSIONS ---\n');
fprintf('e_H = %.4f\n', ys(idy_emissions(1)));
fprintf('e_F = %.4f\n', ys(idy_emissions(2)));

fprintf('\n--- PARAMETRES CALIBRES ---\n');
for i = 1:numel(idp_control)
    fprintf('%s = %.6f\n', M_.param_names{idp_control(i)}, M_.params(idp_control(i)));
end