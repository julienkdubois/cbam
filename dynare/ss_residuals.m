function residuals = ss_residuals(x, idp_control, idy_target, my_targets_ss, M_, oo_)
    M_.params(idp_control) = x;
    try
        [ss, M_.params, check] = eval([ M_.fname '.steadystate(oo_.steady_state, transpose(oo_.exo_steady_state), M_.params)']);
        % Si Dynare échoue (check > 0) ou produit des NaNs, renvoyer une pénalité
        if check || any(isnan(ss))
            residuals = ones(size(idy_target)) * 1e12; 
        else
            % real() est indispensable pour lsqnonlin
            residuals = real(ss(idy_target) - my_targets_ss(:)); 
        end
    catch
        residuals = ones(size(idy_target)) * 1e12; 
    end
end