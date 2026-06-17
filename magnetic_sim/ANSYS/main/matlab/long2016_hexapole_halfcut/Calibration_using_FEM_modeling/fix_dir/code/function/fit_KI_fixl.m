function [ell, gB, Khat, J] = fit_KI_fixl(coil, dhat, I)
%FIT_KI_FIXL  Fit {Khat_I^FEM (6x6, (1,1)=5/6), ell, gB} to the FEM field (fix-l model).
%   [ell, gB, Khat, J] = FIT_KI_FIXL(coil, dhat, I)
%   Document fix-l point-charge model: the 6 charges are locked on the pole axes at
%   ell*dhat (no bias); free parameters are ell, gB and the 35 off-fixed Khat entries.
%   Solved by lsqnonlin on charge_residual. Returns ell [m], gB, Khat (6x6), J (cost).
%   Initial guess per the document: Khat0 = eye(6)-ones(6)/6, ell0 = 0.5 mm, gB0 = 1.
%   (Named *_fixl to avoid clashing with hexapole-long2016\analysis\fit_charge_model.m.)
    Khat0 = eye(6) - ones(6)/6;        % diag 5/6, off-diag -1/6
    ell0  = 0.5e-3;                    % 0.5 mm
    gB0   = 1;
    freemask = true(6);  freemask(1,1) = false;     % Khat(1,1) fixed = 5/6
    opts = optimoptions('lsqnonlin','Display','off', ...
        'MaxFunctionEvaluations',1e5,'MaxIterations',4e3, ...
        'FunctionTolerance',1e-20,'StepTolerance',1e-12);
    x0 = [ell0*1e3; gB0; Khat0(freemask)];
    xf = lsqnonlin(@(x) charge_residual(x, coil, dhat, I, freemask), x0, [], [], opts);
    [ell, gB, Khat] = unpack_params(xf, freemask);
    J = sum(charge_residual(xf, coil, dhat, I, freemask).^2);
end
