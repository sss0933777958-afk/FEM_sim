function [ell, ghat_I_B, K_bar, J] = fit_KI_fixl(coil, dhat, I)
%FIT_KI_FIXL  Fit {K_bar (6x6, (1,1)=5/6), ell, ^Bg_I} to the FEM field (fix-l model).
%   [ell, ghat_I_B, K_bar, J] = FIT_KI_FIXL(coil, dhat, I)
%   Document fix-l point-charge model: the 6 charges are locked on the pole axes at
%   ell*dhat (no bias); free parameters are ell, ^Bg_I (ghat_I_B) and the 35 off-fixed K_bar entries.
%   Solved by lsqnonlin on charge_residual. Returns ell [m], ghat_I_B, K_bar (6x6), J (cost).
%   Initial guess per the document: K_bar0 = eye(6)-ones(6)/6, ell0 = 0.5 mm, ghat_I_B0 = 1.
%   Units: fit in SI (coil.p [m], ell [m]; well-scaled optimizer packs ell in mm); coil.bfem [mT]
%     -> ghat_I_B [mT/A]. 呼叫端把回傳 ell(m) ×1e6 成 µm（Unit Sheet）。
%   (Named *_fixl to avoid clashing with hexapole-long2016\analysis\fit_charge_model.m.)
    K_bar0     = eye(6) - ones(6)/6;   % diag 5/6, off-diag -1/6
    ell0       = 0.5e-3;               % 0.5 mm [m]
    ghat_I_B0  = 1;
    freemask = true(6);  freemask(1,1) = false;     % K_bar(1,1) fixed = 5/6
    opts = optimoptions('lsqnonlin','Display','off', ...
        'MaxFunctionEvaluations',1e5,'MaxIterations',4e3, ...
        'FunctionTolerance',1e-20,'StepTolerance',1e-12);
    x0 = [ell0*1e3; ghat_I_B0; K_bar0(freemask)];   % ell 打包成 mm（well-scaled）
    xf = lsqnonlin(@(x) charge_residual(x, coil, dhat, I, freemask), x0, [], [], opts);
    [ell, ghat_I_B, K_bar] = unpack_params(xf, freemask);
    J = sum(charge_residual(xf, coil, dhat, I, freemask).^2);
end
