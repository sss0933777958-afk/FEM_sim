function [l_hat, e_hat, H_hat, J, out] = fitting_force(f_m, P, Pc_base, l0, H0, u, USE_BIAS, NITER)
%FITTING_FORCE  Minimise J(l, e, H) = sum over (i,j,k) of |f_m - f_model|^2.
%
%   [l_hat, e_hat, H_hat, J, out] = FITTING_FORCE(f_m, P, Pc_base, l0, H0, u, USE_BIAS)
%     f_m      : 3 x Np x N_I  measured force [pN]
%     P        : Np x 3        sample positions [um], actuator frame
%     Pc_base  : 3 x 6         ideal charge grid
%     l0       : initial effective length [um]     (l_design)
%     H0       : 6 x 6 initial force transfer matrix  (eye(6) by decision)
%     u        : 6 x N_I excitation matrix, column j is u_j; default eye(6).
%                Current base: Fmap [A]. Voltage base: V [mV]. Taken from the
%                flux calibration by the caller, not assumed here.
%     USE_BIAS : true  -> solve l + e(17) + H(36) = 54 unknowns  (default)
%                false -> freeze e = 0, solve l + H(36) = 37 unknowns
%     NITER    : optional fixed iteration count. Given, the tolerances are
%                switched off and the optimiser runs exactly that many steps,
%                so two runs can be compared over the same horizontal axis.
%                Omitted, it stops on its own convergence criteria.
%
%     l_hat : [um]     e_hat : 17x1      H_hat : 6x6      J : residual sum of squares
%     out   : nvar / USE_BIAS / npts / NI / resnorm / rms / elapsed / Jhist
%             Jhist is the cost per iteration (column 1 = iteration index,
%             column 2 = J), recorded through an OutputFcn so the convergence
%             can be plotted afterwards.
%
%   Unlike the flux-side fitting, the charges CANNOT be projected out here: the
%   force is quadratic in H, so H stays in the optimisation vector. That is why
%   this solves 54 unknowns where the flux fit solves 18.
%
%   The residual caches L across the H-only perturbations of the finite-
%   difference Jacobian: L depends on (l, e) alone, which is 18 of the 54
%   variables, so the cache removes roughly two thirds of the build_L calls.
%
%   BOUNDS.  l is bounded below at 1 um. Without it the optimiser can and does
%   step l straight through zero (observed on zhi_peng single, where H0 = I is
%   ~5x below the true H and the first steps are large), which makes the
%   dimensionless position P/l blow up and build_L reject the argument.
%
%   Packing: x = [ l/1e3 ; e(17) ; H(:) (36) ].  l is divided by 1e3 (um -> mm) so
%   the optimiser sees ~0.87 instead of ~870, putting it on the same scale as e
%   and H. Same trick as the flux fit.
%
%   See also BUILD_L, CURRENT_BASE_MODEL, SOLVE_FORCE_CURRENT.

    if nargin < 6 || isempty(u),   u   = eye(6); end
    if nargin < 7 || isempty(USE_BIAS), USE_BIAS = true;  end
    if nargin < 8, NITER = []; end
    validateattributes(P, {'numeric'}, {'2d','ncols',3,'real','finite'}, mfilename, 'P');
    validateattributes(H0, {'numeric'}, {'size',[6 6],'real','finite'}, mfilename, 'H0');
    validateattributes(l0, {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'l0');
    Np = size(P,1);   NI = size(u,2);
    assert(isequal(size(f_m), [3 Np NI]), 'fitting_force:fmSize', ...
           'f_m must be 3 x %d x %d (got %s)', Np, NI, mat2str(size(f_m)));

    JH = zeros(0, 2);                                 % [iteration, J] per step
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'MaxFunctionEvaluations',1e6, 'MaxIterations',4e3, ...
        'FunctionTolerance',1e-20, 'StepTolerance',1e-12, ...
        'OutputFcn', @track);
    if ~isempty(NITER)
        opts = optimoptions(opts, 'MaxIterations',NITER, ...
            'FunctionTolerance',0, 'StepTolerance',0, 'OptimalityTolerance',0);
    end

    % L cache, shared by every residual evaluation in this call
    Lc = [];   gec = [];

    t0 = tic;
    if USE_BIAS
        x0 = [l0/1e3; zeros(17,1); H0(:)];
        lb = -inf(size(x0));   lb(1) = 1e-3;          % l >= 1 um
        xf = lsqnonlin(@resid, x0, lb, [], opts);
        e_hat = xf(2:18);
    else
        x0 = [l0/1e3; H0(:)];
        lb = -inf(size(x0));   lb(1) = 1e-3;          % l >= 1 um
        xf = lsqnonlin(@(x) resid([x(1); zeros(17,1); x(2:37)]), x0, lb, [], opts);
        xf = [xf(1); zeros(17,1); xf(2:37)];
        e_hat = zeros(17,1);
    end
    l_hat = xf(1)*1e3;
    H_hat = reshape(xf(19:54), 6, 6);
    r     = resid(xf);
    J     = sum(r.^2);

    [~, keep] = unique(JH(:,1), 'first');             % 'init' and the first 'iter' both report 0
    JH = JH(sort(keep), :);
    out = struct('nvar', numel(x0), 'USE_BIAS', USE_BIAS, 'npts', Np, 'NI', NI, ...
                 'resnorm', J, 'rms', sqrt(J/numel(r)), 'elapsed', toc(t0), 'Jhist', JH);

    % ---- convergence history -------------------------------------------
    function stop = track(~, ov, state)
        stop = false;
        if any(strcmp(state, {'init','iter'}))
            JH(end+1, :) = [ov.iteration, ov.resnorm];
        end
    end

    % ---- residual, with L cached over the (l,e) block --------------------
    function rr = resid(x)
        ge = x(1:18);
        if isempty(gec) || ~isequal(ge, gec)
            Lc  = build_L(P, x(1)*1e3, x(2:18), Pc_base);
            gec = ge;
        end
        f_model = base_model(Lc, reshape(x(19:54), 6, 6), u);
        rr = f_m(:) - f_model(:);
    end
end
