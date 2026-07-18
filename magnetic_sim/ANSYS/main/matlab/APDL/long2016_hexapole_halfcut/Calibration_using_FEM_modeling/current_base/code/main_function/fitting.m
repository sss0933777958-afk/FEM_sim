function [ell, e_hat, Pc, J] = fitting(P, Bstack, Pc_base, ell0, USE_BIAS)
%FITTING  Variable-projection fit of the point-charge model geometry.
%   [ell, e_hat, Pc, J] = FITTING(P, Bstack, Pc_base, ell0, USE_BIAS)
%   Iterates the nonlinear parameters {ell (+ e_hat(17) if USE_BIAS)}; the 6 per-
%   simulation charge magnitudes g_j are profiled out by least squares inside the
%   residual (variable projection).  Returns the fitted ell [m], e_hat (17x1; zeros
%   when ~USE_BIAS), the assembled charge lattice Pc = Pc_base + E(e_hat) (3x6), and
%   J (cost = sum of squared residuals, mT^2).
%     USE_BIAS=true  -> {ell, e_hat(17)}  (off-axis 18-param bias model)
%     USE_BIAS=false -> {ell}             (on-axis fix model, e_hat frozen 0, 1-D)
%   Units: fit in SI (P [m], ell [m]; optimizer packs ell in mm = well-scaled vs
%   bias e~O(0.01)); Bstack [mT]; e_hat dimensionless.
%   (Merges the old fit_varpro + bias_resid + make_Pc into one file, 2026-07-16.)
    if nargin < 5, USE_BIAS = true; end
    opts = optimoptions('lsqnonlin','Display','off', ...
        'MaxFunctionEvaluations',1e5,'MaxIterations',4e3, ...
        'FunctionTolerance',1e-20,'StepTolerance',1e-12);
    if USE_BIAS
        x0 = [ell0*1e3; zeros(17,1)];                       % ell(mm) + bias e(17)
        xf = lsqnonlin(@(x) resid(x, P, Bstack, Pc_base), x0, [], [], opts);
        ell   = xf(1)*1e-3;                                 % mm -> m
        e_hat = xf(2:18);
    else
        x0 = ell0*1e3;                                      % only ell (mm); e frozen 0
        xf = lsqnonlin(@(x) resid([x; zeros(17,1)], P, Bstack, Pc_base), x0, [], [], opts);
        ell   = xf(1)*1e-3;
        e_hat = zeros(17,1);
    end
    Pc = make_Pc(e_hat, Pc_base);
    J  = sum(resid([ell*1e3; e_hat], P, Bstack, Pc_base).^2);
end

% ---- lsqnonlin residual: g_j profiled out per simulation (variable projection) ----
function r = resid(x, P, Bstack, Pc_base)
    A = build_S_matrix(x(1)*1e-3, make_Pc(x(2:18), Pc_base), P);   % 3Np x 6
    M = A.' * A;                                                   % 6 x 6
    r = [];
    for j = 1:size(Bstack,2)
        gj = M \ (A.' * Bstack(:,j));
        r  = [r; A*gj - Bstack(:,j)]; %#ok<AGROW>
    end
end

% ---- assemble charge lattice Pc = Pc_base + E(e17), with e6z constraint ----
function Pc = make_Pc(e17, Pc_base)
    E = zeros(3, 6);
    E(:,1) = e17(1:3);    E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);    E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);     E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);   % e6z constraint (paper step 2c)
    Pc = Pc_base + E;
end
