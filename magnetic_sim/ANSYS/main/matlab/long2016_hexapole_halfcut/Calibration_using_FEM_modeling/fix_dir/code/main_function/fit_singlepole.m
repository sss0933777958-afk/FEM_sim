function [ell, gI, J, res, Bmod] = fit_singlepole(P, Bfem, dhat, I)
%FIT_SINGLEPOLE  No-bias single-pole point-charge fit: fit {ell, ghat_I} (2 scalars).
%   [ell, gI, J, res, Bmod] = FIT_SINGLEPOLE(P, Bfem, dhat, I)
%   Model (single effective monopole, no K_bar, scalar current I):
%       B(p) = ghat_I * I * (p/ell - dhat) / |p/ell - dhat|^3       [mT]
%   with the single charge on-axis at ell*dhat (no bias). Frame-invariant → run in the
%   actuator frame (dhat = e1). Solved by lsqnonlin (well-scaled: ell packed in mm).
%   Inputs: P Nx3 [m], Bfem Nx3 [mT] (all-source, actuator frame), dhat 3x1 unit, I [A].
%   Returns: ell [m], gI [mT/A], J = cost = sum(res.^2) [mT^2] (per-node actuator-frame
%            residual B_model-B_FEM), res (stacked residual), Bmod.
    x0   = [0.5; 1];                  % [ell_mm ; gI]
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'FunctionTolerance',1e-18, 'StepTolerance',1e-12, ...
        'MaxFunctionEvaluations',1e4, 'MaxIterations',2e3);
    xf = lsqnonlin(@(x) sp_res(x, P, Bfem, dhat, I), x0, [], [], opts);
    ell = xf(1)*1e-3;  gI = xf(2);
    [res, Bmod] = sp_res(xf, P, Bfem, dhat, I);
    J = sum(res.^2);                  % cost（lsqnonlin 目標）= Σ|B_model-B_FEM|²
end

function [r, Bmod] = sp_res(x, P, Bfem, dhat, I)
    ell = x(1)*1e-3;  gI = x(2);
    D    = P/ell - dhat.';               % Nx3 (broadcast dhat row)
    nrm  = vecnorm(D, 2, 2);             % Nx1
    Bmod = gI * I * D ./ (nrm.^3);       % Nx3  [mT]
    r    = Bmod - Bfem;  r = r(:);
end
