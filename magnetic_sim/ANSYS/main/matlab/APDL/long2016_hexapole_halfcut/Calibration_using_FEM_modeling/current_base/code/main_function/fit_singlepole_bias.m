function [ell, ey, ez, gI, J, res, Bmod] = fit_singlepole_bias(P, Bfem, I)
%FIT_SINGLEPOLE_BIAS  Single-pole point-charge fit WITH transverse bias: {ell, ey, ez, ghat_I}.
%   [ell, ey, ez, gI, J, res, Bmod] = FIT_SINGLEPOLE_BIAS(P, Bfem, I)
%   Bias analog of fit_singlepole (current_base): the charge is let off-axis by a transverse
%   offset (e_y, e_z) while the axial component stays gauge-fixed at 1 (breaks the
%   ell<->gI scale degeneracy, mirrors the hexapole no_fix_l 18-param bias):
%       pc_hat = [1 ; e_y ; e_z]                          (dimensionless charge position)
%       B(p)   = ghat_I * I * (p/ell - pc_hat) / |p/ell - pc_hat|^3    [mT]
%   ell captures the axial charge depth (e.g. tipcut), (e_y,e_z) captures shape asymmetry
%   (e.g. halfcut z-cut -> e_z != 0). Frame-invariant -> run in the +x mesh frame
%   (pole axis = +x, WP at origin). Solved by lsqnonlin (ell packed in mm; e ~ O(0.01)).
%   Inputs: P Nx3 [m], Bfem Nx3 [mT] (all-source, +x frame), I [A].
%   Returns: ell [m], e_y, e_z (dimensionless), gI [mT/A],
%            J = sum(res.^2) [mT^2], res (stacked residual), Bmod Nx3 [mT].
    x0   = [0.5; 0; 0; 1];            % [ell_mm ; e_y ; e_z ; gI]
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'FunctionTolerance',1e-18, 'StepTolerance',1e-12, ...
        'MaxFunctionEvaluations',2e4, 'MaxIterations',4e3);
    xf = lsqnonlin(@(x) spb_res(x, P, Bfem, I), x0, [], [], opts);
    ell = xf(1)*1e-3;  ey = xf(2);  ez = xf(3);  gI = xf(4);
    [res, Bmod] = spb_res(xf, P, Bfem, I);
    J = sum(res.^2);                 % cost（lsqnonlin 目標）= Σ|B_model-B_FEM|²
end

function [r, Bmod] = spb_res(x, P, Bfem, I)
    ell = x(1)*1e-3;  pchat = [1, x(2), x(3)];     % [1, e_y, e_z]
    D    = P/ell - pchat;                          % Nx3 (broadcast)
    nrm  = vecnorm(D, 2, 2);                        % Nx1
    Bmod = x(4) * I * D ./ (nrm.^3);               % Nx3  [mT]
    r    = Bmod - Bfem;  r = r(:);
end
