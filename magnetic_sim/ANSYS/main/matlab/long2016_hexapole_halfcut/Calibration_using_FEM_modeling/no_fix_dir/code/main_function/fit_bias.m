function [ell, e_hat, J] = fit_bias(P, Bstack, Pc_base, ell0)
%FIT_BIAS  Fit the 18-param bias model {ell, e_hat(17)} by lsqnonlin (g_j profiled out).
%   [ell, e_hat, J] = FIT_BIAS(P, Bstack, Pc_base, ell0)
%   Initial guess per no_fix_l.pdf step 5b: ell = ell0 (= ell_design), e_hat = 0.
%   Uses lsqnonlin on bias_resid (avoids the R2025b fminunc optim:fminusub bug).
%   Units: fit in SI (P [m], ell [m]; optimizer packs ell in mm = well-scaled vs bias e~O(0.01)).
%     Bstack [mT], e_hat dimensionless. 呼叫端把回傳 ell(m) ×1e6 成 µm（Unit Sheet）。
%   Returns ell [m], e_hat (17x1), J (cost = sum of squared residuals, mT^2).
    opts = optimoptions('lsqnonlin','Display','off', ...
        'MaxFunctionEvaluations',1e5,'MaxIterations',4e3, ...
        'FunctionTolerance',1e-20,'StepTolerance',1e-12);
    x0 = [ell0*1e3; zeros(17,1)];                       % ell 打包成 mm（well-scaled），bias e = 0
    xf = lsqnonlin(@(x) bias_resid(x, P, Bstack, Pc_base), x0, [], [], opts);
    ell   = xf(1)*1e-3;                                 % mm -> m
    e_hat = xf(2:18);
    J     = sum(bias_resid(xf, P, Bstack, Pc_base).^2);
end
