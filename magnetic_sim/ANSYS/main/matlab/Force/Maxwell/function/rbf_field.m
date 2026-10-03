function [b, grad, info] = rbf_field(R, kern)
%RBF_FIELD  Radial-basis interpolant of the six FEM field exports.
% =========================================================================
%   One call builds the interpolant and hands back the two evaluators, so the
%   expensive part -- one Phi, one factorisation, 18 right-hand sides -- happens
%   once no matter how many points are asked for afterwards.
%
%       [b, grad]       = rbf_field(R, kern)              one data range
%       [bs, gs, infos] = rbf_field({R1 R2 ...}, kern)    sweep the data range
%
%   INPUT 1  資料範圍  R
%     scalar   [um]   the nodes with |p| <= R become the RBF centres
%     [r1 r2]  [um]   an annulus instead, r1 <= |p| <= r2
%     cell            SWEEP: one model per entry, each entry being either of the
%                     above.  A numeric R keeps its meaning exactly, so the cell
%                     is what tells a sweep apart from an annulus -- {300 500} is
%                     two models, [300 500] is one annulus.
%
%   INPUT 2  核函數  kern   (struct; a bare scalar is read as rho with lam = 0)
%     .rho     [um]   shape parameter
%     .lam            ridge:  (Phi + lam*I) W = B.   lam = 0 is exact interpolation
%     .name           'gauss' (default) | 'mq' | 'imq'
%     .phi            optional handle phi(r2)            -- overrides .name
%     .dphi_r         optional handle phi'(r)/r of r2    -- required with .phi
%     .nb             block size for the build and factorisation (default 1024)
%     .solver         'chol' | 'lu'.  Default picks itself; see MEMORY below.
%
%       gauss   phi = exp(-r^2/rho^2)        phi'(r)/r = -(2/rho^2) phi
%       mq      phi = sqrt(r^2 + rho^2)      phi'(r)/r = 1/sqrt(r^2 + rho^2)
%       imq     phi = 1/sqrt(r^2 + rho^2)    phi'(r)/r = -(r^2 + rho^2)^(-3/2)
%
%     Every kernel is parametrised by phi'(r)/r rather than phi'(r), which is
%     what removes the 0/0 at r = 0 -- a query point may sit exactly on a node.
%     The kernel is held FIXED across a sweep; that is the point of sweeping R.
%
%   OUTPUT 1  b(p)     -- a HANDLE, so it reads the way it is written
%     b(Pq)               -> Nq x 3      [mT]   Pq in [um], ACTUATOR frame
%     b(Pq, I)            -> Nq x 3             for the 6-coil current vector I [A]
%
%   OUTPUT 2  那點的梯度  -- also a handle
%     grad(Pq)            -> Nq x 3      [mT^2/um]   column j = d(b.b)/dx_j
%     grad(Pq, I)         -> Nq x 3                  for current vector I
%     grad(Pq, I, 'jac')  -> Nq x 3 x 3  [mT/um]     J(:,c,j) = db_c/dx_j
%     grad(Pq, I, 'd2')   -> Nq x 3      [mT^2/um^2]  column j = d2(b.b)/dx_j^2
%     grad(Pq, I, 'd3')   -> Nq x 3      [mT^2/um^3]  column j = d3(b.b)/dx_j^3
%
%   'd2' is the analytic second derivative, not a difference of the sampled
%   first one, so it does not depend on how densely Pq is laid out:
%       d2phi/dxj^2   = dphi_r + 2*(xj - xjn)^2 * ddphi,  ddphi = d(phi'(r)/r)/d(r2)
%       d2(b.b)/dxj^2 = 2*sum_c [ (dbc/dxj)^2 + bc * d2bc/dxj^2 ]
%       d3phi/dxj^3   = 6*(xj-xjn)*ddphi + 4*(xj-xjn)^3*dddphi
%       d3(b.b)/dxj^3 = 2*sum_c [ 3*(dbc/dxj)*(d2bc/dxj^2) + bc * d3bc/dxj^3 ]
%   'd3' is the sharper smoothness test.  d2 carries an unconditionally positive
%   term, 2*sum_c (dbc/dxj)^2, which at R = 250 um is 40 % of it and acts as a
%   floor that keeps d2 from ever changing sign -- so d2 > 0 is just a restatement
%   of F being monotone and cannot see a ripple.  Neither term of d3 is sign
%   definite, and a ripple of wavelength L is amplified by ~2*pi/L each time the
%   curve is differentiated, so d3 is where a wave in F shows up.
%   Its sign is the smoothness test used in this study: on the P1 axis over
%   |x_a| <= 150 um the true F_xa rises with increasing slope throughout, so
%   b.b is convex there and ANY point with d2 < 0 is an artefact of the fit,
%   not physics.  That is a fact about this segment and this excitation, not a
%   general rule -- elsewhere an inflection can be real.
%
%   In sweep mode both come back as cells, so bs{k}(Pq) is entry k's field and
%   infos(k) its record; an entry that failed has [] for its handles and
%   infos(k).ok false with the message in infos(k).err, which keeps a long sweep
%   from losing the levels that did succeed.
%
%   I defaults to [1 0 0 0 0 0], i.e. coil 1 (pole P1) at 1 A, the excitation
%   every figure in this study used.  The FIELD superposes because mu_r = 280 is
%   constant, so W(I) = sum_k I_k W_k and the six weight sets collapse to one
%   once the current vector is given.  The FORCE does not superpose: it is
%   quadratic in b.
%
%   The project's force follows from the gradient in one line,
%       dU = grad(Pq);   F_j = 0.5 * mgB * UF * dU(:,j)
%   with mgB = 0.0451 and UF = 1e3, giving F in pN.
%
%   OUTPUT 3  info  -- the built model itself, plus its arithmetic checks.
%       P B6 W6 Np R kern map R_act    the model
%       wmax                           max|W|
%       resid                          max|Phi W - B| [mT], from the EXACT identity
%                                      (Phi + lam I) W = B  =>  Phi W - B = -lam W,
%                                      so resid = lam*wmax and no second product
%                                      of Phi is needed
%       pivmin pivmax                  smallest / largest Cholesky pivot
%       cndest                         (pivmax/pivmin)^2, a conditioning proxy
%       solver secs ok err             which path ran, wall time, and how it went
%   pivmin / cndest replace the old rcnd, which needed a second full copy of Phi.
%
%   ── MEMORY, and why the solve is written the way it is ───────────────────
%   Phi is Np x Np dense and Np grows as R^3 on the 20 um export grid:
%       R   250     300      350      400      450      500   um
%       Np  8225    14082    22455    33460    47693    65362
%       Phi 0.5     1.6      4.0      9.0      18.2     34.2  GB
%   The straightforward path costs TWICE that, because `Phi + lam*eye(Np)` makes
%   a second dense matrix and `A\B` copies A again before factorising -- 68 GB at
%   R = 500 on a 64 GB machine.  Three changes bring the peak back to one Phi:
%     1. Phi is written column block by column block into a preallocated array,
%        so no Np x Np temporary is ever formed (the Gram form s2 + s2' - 2*P*P'
%        makes at least two).
%     2. lam goes onto the diagonal in place; the diagonal is exactly phi(0)+lam.
%     3. A blocked Cholesky overwrites Phi's lower triangle, and the two
%        triangular solves use linsolve's LT/TRANSA options, which read only that
%        triangle -- tril(Phi) would allocate another full copy.  The trailing
%        update runs in column panels so its temporary stays Np x nb.
%   The stale upper triangle left behind by step 3 is never read.  A sweep frees
%   Phi before starting the next entry, so its peak is that of its largest entry,
%   not the sum; what survives per entry is only P and W6, a few MB each.
%
%   Step 3 needs Phi POSITIVE DEFINITE.  Gaussian and inverse multiquadric are
%   strictly positive definite, so they take that path.  Plain multiquadric is
%   only CONDITIONALLY positive definite -- its Phi carries a negative eigenvalue
%   and chol would fail -- so 'mq' and any custom kernel fall back to the LU
%   path, which still builds Phi in blocks but pays the one extra copy in the
%   solve.  Override with kern.solver if you know better.
%
%   Frame and units are the study's convention throughout: positions in um and
%   fields in mT, both rotated into the actuator frame by R_act and translated
%   so the six pole tips are centred on the origin.  The six .fld files are read
%   once, at the largest radius the call asks for, and cached for later calls.
%
%   ⚠ The outer shells of the larger radii are where the field explodes.  Per
%   50 um shell the maximum |B| runs 11.2 mT at the centre, 40.6 at 350-400,
%   72.1 at 400-450, 260.9 at 450-500, then 2722 at 500-550 -- that last one is
%   iron, so R_norm = 500 um is exactly where the steel starts.  A single global
%   rho tuned for the gentle centre cannot also resolve a tip whose own radius is
%   40 um, so expect the outer nodes to demand large weights as R grows.
%
%   Example
%     kern = struct('rho',120, 'lam',1e-4);
%     [b, grad] = rbf_field(250, kern);
%     xa = linspace(-150,150,1001).';   Pq = [xa, zeros(1001,2)];
%     U  = sum(b(Pq).^2, 2);                        % b.b        [mT^2]
%     dU = grad(Pq);   F = 0.5*0.0451*1e3*dU(:,1);  % F_xa       [pN]
%
%     [bs, gs, infos] = rbf_field(num2cell(300:50:500), kern);   % the range sweep
%     Us = cell2mat(cellfun(@(f) sum(f(Pq).^2,2), bs, 'uni',0)); % 1001 x 5
%
%   [ADDED 2026-09-11] merged from temp_code/scripts/{rbf_w.m, rbf_w6.m,
%   rbf_bb.m}: the weight solve, the evaluation and the analytic gradient were
%   three separate scripts with the Gaussian hard-coded in each.
%   [MODIFIED 2026-09-13] memory-lean build and factorisation, plus the cell
%   sweep, so the data range can be swept past R = 400 um on this machine.
% =========================================================================

    % ---- kernel (resolved once; a sweep holds it fixed) ------------------
    if isnumeric(kern), kern = struct('rho', kern, 'lam', 0); end
    if ~isfield(kern,'lam')  || isempty(kern.lam),  kern.lam  = 0;       end
    if ~isfield(kern,'name') || isempty(kern.name), kern.name = 'gauss'; end
    if ~isfield(kern,'nb')   || isempty(kern.nb),   kern.nb   = 1024;    end
    assert(isfield(kern,'rho') && isscalar(kern.rho) && kern.rho > 0, ...
           'kern.rho [um] is required and must be positive');
    assert(isscalar(kern.lam) && kern.lam >= 0, 'kern.lam must be a non-negative scalar');

    if isfield(kern,'phi') && ~isempty(kern.phi)
        assert(isfield(kern,'dphi_r') && ~isempty(kern.dphi_r), ...
               'a custom kern.phi also needs kern.dphi_r = phi''(r)/r as a function of r^2');
        phi = kern.phi;   dphi_r = kern.dphi_r;   kern.name = 'custom';   pd = false;
        if isfield(kern,'ddphi') && ~isempty(kern.ddphi)
            ddphi = kern.ddphi;
        else
            ddphi = @(r2) error(['a custom kernel needs kern.ddphi = ' ...
                                 'd(phi''(r)/r)/d(r^2) before ''d2'' can be used']);
        end
        if isfield(kern,'dddphi') && ~isempty(kern.dddphi)
            dddphi = kern.dddphi;
        else
            dddphi = @(r2) error(['a custom kernel needs kern.dddphi = ' ...
                                  'd2(phi''(r)/r)/d(r^2)^2 before ''d3'' can be used']);
        end
    else
        c2 = kern.rho^2;
        switch lower(kern.name)
        case 'gauss'
            phi    = @(r2) exp(-r2 / c2);
            dphi_r = @(r2) -(2/c2) * exp(-r2 / c2);
            ddphi  = @(r2)  (2/c2^2) * exp(-r2 / c2);
            dddphi = @(r2) -(2/c2^3) * exp(-r2 / c2);    pd = true;
        case 'mq'
            phi    = @(r2) sqrt(r2 + c2);
            dphi_r = @(r2) 1 ./ sqrt(r2 + c2);
            ddphi  = @(r2) -0.5  * (r2 + c2).^(-1.5);
            dddphi = @(r2)  0.75 * (r2 + c2).^(-2.5);    pd = false;   % only conditionally PD
        case 'imq'
            phi    = @(r2) 1 ./ sqrt(r2 + c2);
            dphi_r = @(r2) -(r2 + c2).^(-1.5);
            ddphi  = @(r2)  1.5  * (r2 + c2).^(-2.5);
            dddphi = @(r2) -3.75 * (r2 + c2).^(-3.5);    pd = true;
        otherwise
            error('unknown kernel ''%s'' -- use gauss | mq | imq, or supply phi + dphi_r', ...
                  kern.name);
        end
    end
    if ~isfield(kern,'solver') || isempty(kern.solver)
        if pd, kern.solver = 'chol'; else, kern.solver = 'lu'; end
    end
    assert(any(strcmpi(kern.solver, {'chol','lu'})), 'kern.solver must be ''chol'' or ''lu''');

    % ---- the data ranges to build ---------------------------------------
    sweep = iscell(R);
    if sweep, RL = R(:).';   else,   RL = {R};   end
    assert(~isempty(RL), 'R must name at least one data range');
    RMAX = 0;
    for k = 1:numel(RL)
        Rk = RL{k};
        assert(isnumeric(Rk) && ~isempty(Rk) && numel(Rk) <= 2 && all(Rk > 0), ...
               'entry %d: each data range is a positive scalar or [r1 r2], in um', k);
        RMAX = max(RMAX, max(Rk(:)));
    end

    % ---- nodes (cached across calls; the .fld read is the slow part) -----
    persistent C
    if isempty(C) || C.R < RMAX
        C = read_exports_(max(RMAX, 250));
    end

    % ---- one model per entry; Phi is freed before the next --------------
    n    = numel(RL);
    b    = cell(1,n);   grad = cell(1,n);
    info = repmat(blank_info_(RL{1}, kern), 1, n);
    for k = 1:n
        if sweep
            fprintf('rbf_field: entry %d of %d, R = %s um ...\n', k, n, mat2str(RL{k}));
        end
        try
            [b{k}, grad{k}, info(k)] = build_one_(C, RL{k}, kern, phi, dphi_r, ddphi, dddphi);
        catch ME
            % [MODIFIED 2026-09-14] a SWEEP keeps going and records the failure, so one
            % bad level does not throw away the others.  A single-R call must RETHROW:
            % returning empty handles made the caller fail later with an unrelated
            % "Index in position 1 is invalid" instead of the real cause.
            if ~sweep, rethrow(ME); end
            info(k)     = blank_info_(RL{k}, kern);
            info(k).err = ME.message;
            warning('rbf_field:entryFailed', 'R = %s um failed: %s', ...
                    mat2str(RL{k}), ME.message);
        end
        if sweep
            s = info(k);
            if s.ok
                fprintf(['    Np %d | %.1f s | max|W| %.3e | resid %.3e mT | cndest %.3e' newline], ...
                        s.Np, s.secs, s.wmax, s.resid, s.cndest);
            else
                fprintf(['    FAILED' newline]);
            end
        end
    end
    if ~sweep
        b = b{1};   grad = grad{1};   info = info(1);
    end
end

% =========================================================================
function [b, grad, info] = build_one_(C, R, kern, phi, dphi_r, ddphi, dddphi)
% Select the nodes, build Phi in place, factorise and solve the 18 right-hand
% sides, and wrap the result in the two evaluators.
    tic;
    rn = vecnorm(C.P, 2, 2);
    if isscalar(R), sel = rn <= R;
    else,           sel = rn >= min(R) & rn <= max(R);
    end
    P  = C.P(sel,:);   B6 = C.B6(sel,:,:);   Np = size(P,1);
    assert(Np > 0, 'no nodes fall inside the requested data range');
    nb = min(kern.nb, Np);

    % ---- Phi, one column block at a time, straight into place ------------
    Phi = zeros(Np, Np);
    for j0 = 1:nb:Np
        j  = j0:min(j0+nb-1, Np);
        d2 = (P(:,1)-P(j,1).').^2 + (P(:,2)-P(j,2).').^2 + (P(:,3)-P(j,3).').^2;
        Phi(:,j) = phi(d2);
        clear d2
    end
    Phi(1:Np+1:end) = phi(0) + kern.lam;                 % ridge, in place

    % ---- factorise and solve --------------------------------------------
    pivmin = NaN;   pivmax = NaN;
    if strcmpi(kern.solver, 'chol')
        for k0 = 1:nb:Np                                 % blocked, in place
            kb = min(nb, Np-k0+1);
            j  = k0:k0+kb-1;
            Phi(j,j) = chol(Phi(j,j), 'lower');
            dg = diag(Phi(j,j));
            pivmin = min([pivmin; dg]);   pivmax = max([pivmax; dg]);
            if k0+kb <= Np
                r = k0+kb:Np;
                Phi(r,j) = Phi(r,j) / Phi(j,j).';
                for c0 = k0+kb:nb:Np                     % panel the trailing update so
                    c  = c0:min(c0+nb-1, Np);            % the temporary stays Np x nb
                    rr = c0:Np;
                    Phi(rr,c) = Phi(rr,c) - Phi(rr,j) * Phi(c,j).';
                end
            end
        end
        oL = struct('LT',true);   oT = struct('LT',true,'TRANSA',true);
        W  = linsolve(Phi, linsolve(Phi, reshape(B6, Np, 18), oL), oT);
    else
        W  = Phi \ reshape(B6, Np, 18);                  % one extra copy of Phi
    end
    clear Phi
    W6 = reshape(W, Np, 3, 6);   clear W

    b    = @(varargin) eval_b_(P, W6, phi,         varargin{:});
    grad = @(varargin) eval_g_(P, W6, phi, dphi_r, ddphi, dddphi, varargin{:});

    info = blank_info_(R, kern);
    info.P    = P;    info.B6     = B6;     info.W6   = W6;    info.Np = Np;
    info.map  = C.map;                      info.R_act = C.R_act;
    info.wmax = max(abs(W6), [], 'all');    info.resid = kern.lam * info.wmax;
    info.pivmin = pivmin;  info.pivmax = pivmax;  info.cndest = (pivmax/pivmin)^2;
    info.secs = toc;       info.ok = true;
end

function s = blank_info_(R, kern)
% One constructor for both the success and the failure record, so a sweep's
% struct array always has the same fields in the same order.
    s = struct('P',[], 'B6',[], 'W6',[], 'Np',0, 'R',R, 'kern',kern, ...
               'map',[], 'R_act',[], 'wmax',NaN, 'resid',NaN, ...
               'pivmin',NaN, 'pivmax',NaN, 'cndest',NaN, ...
               'solver',kern.solver, 'secs',NaN, 'ok',false, 'err','');
end

% =========================================================================
function bq = eval_b_(P, W6, phi, Pq, I)
    if nargin < 5, I = []; end
    W  = wsum_(W6, I);
    bq = phi(dist2_(Pq, P)) * W;
end

function out = eval_g_(P, W6, phi, dphi_r, ddphi, dddphi, Pq, I, what)
    if nargin < 8, I    = [];   end
    if nargin < 9, what = 'du'; end
    assert(any(strcmpi(what, {'du','jac','d2','d3'})), ...
           'the third argument must be ''du'' (default), ''jac'', ''d2'' or ''d3''');
    W  = wsum_(W6, I);
    D2 = dist2_(Pq, P);
    K  = dphi_r(D2);                                  % phi'(r)/r, finite at r = 0
    Nq = size(Pq,1);
    J  = zeros(Nq, 3, 3);                             % J(:,c,j) = db_c/dx_j
    for j = 1:3
        J(:,:,j) = (K .* (Pq(:,j) - P(:,j).')) * W;
    end
    if strcmpi(what, 'jac'), out = J;  return;  end
    bq  = phi(D2) * W;
    out = zeros(Nq, 3);
    if strcmpi(what, 'du')                            % d(b.b)/dx_j = 2 sum_c b_c db_c/dx_j
        for j = 1:3
            out(:,j) = 2 * sum(bq .* J(:,:,j), 2);
        end
        return
    end
    KK = ddphi(D2);                                   % d(phi'(r)/r)/d(r^2)
    if strcmpi(what, 'd2')
        for j = 1:3                                   % d2(b.b)/dx_j^2
            Uj  = Pq(:,j) - P(:,j).';
            d2b = (K + 2*(Uj.^2).*KK) * W;            % d2 b_c / dx_j^2
            out(:,j) = 2 * sum(J(:,:,j).^2 + bq .* d2b, 2);
        end
        return
    end
    KKK = dddphi(D2);                                 % d2(phi'(r)/r)/d(r^2)^2
    for j = 1:3                                       % d3(b.b)/dx_j^3
        Uj  = Pq(:,j) - P(:,j).';
        d2b = (K + 2*(Uj.^2).*KK) * W;
        d3b = (6*Uj.*KK + 4*(Uj.^3).*KKK) * W;        % d3 b_c / dx_j^3
        out(:,j) = 2 * sum(3*J(:,:,j).*d2b + bq .* d3b, 2);
    end
end

function W = wsum_(W6, I)
% Collapse the six weight sets onto one current vector.  Legal because the FIELD
% is linear in the currents (mu_r constant); the force is not.
    if isempty(I), I = [1;0;0;0;0;0]; end             % coil 1 (P1) at 1 A
    I = I(:);
    assert(numel(I) == 6, 'the current vector I needs six entries [A], one per coil');
    W = reshape(reshape(W6, [], 6) * I, size(W6,1), 3);
end

function D2 = dist2_(Pq, P)
% Componentwise, not the Gram form: at a query point sitting on a node the Gram
% form loses the cancellation and phi comes out slightly wrong.
    assert(size(Pq,2) == 3, 'the query points Pq must be Nq x 3, in um');
    D2 = (Pq(:,1) - P(:,1).').^2 + (Pq(:,2) - P(:,2).').^2 + (Pq(:,3) - P(:,3).').^2;
end

% =========================================================================
function C = read_exports_(R)
% Read the six .fld exports once, translate to the tip-sphere centre and rotate
% into the actuator frame.  The field is kept RAW, i.e. the physical field for
% +1 A in the deck's own current direction.  The all-source sign flip of
% charge-model-source-convention.md is a presentation choice for the charge model
% and must NOT be applied here, because superposing b = sum_k I_k b_k needs the
% real per-ampere fields.
    HERE = fileparts(mfilename('fullpath'));                  % .../Maxwell/utils
    FMX  = fileparts(HERE);                                   % .../matlab/Force/Maxwell
    MAIN = fileparts(fileparts(fileparts(FMX)));              % .../main
    CAL  = fullfile(MAIN, 'matlab', 'Flux', 'Maxwell');
    addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'), fullfile(CAL,'common_path'));

    cfg = model_config('long2016_hexapole_halfcut', 'tip40um');
    ZC  = cfg.SPH_OFST;
    % The [1 3 6 5 2 4] map is the APDL deck's build order and does NOT apply to
    % Maxwell: this config declares identity, i.e. export B_p<k> IS pole P<k>
    % (pole-coil-numbering.md -- the map is per-model, never copied between branches).
    map = cfg.apdl_to_paper_idx;   if isempty(map), map = 1:6; end

    tip   = [cfg.pole_tip_x; cfg.pole_tip_y; cfg.pole_tip_z_wp];
    dhat  = tip ./ vecnorm(tip);
    R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';
    assert(abs(det(R_act)-1) < 1e-9 && norm(R_act.'*R_act - eye(3)) < 1e-9, ...
           'R_act is not a rotation -- the pole axes or the magic angle are wrong');

    P = [];  B6 = [];  Np = 0;
    for k = 1:6
        FB = fullfile(cfg.fld_dir, 'baseline', sprintf('B_p%d.fld', k));
        if ~isfile(FB), FB = fullfile(cfg.fld_dir, sprintf('B_p%d.fld', k)); end
        assert(isfile(FB), 'missing export %s', FB);
        A   = readmatrix(FB, 'FileType','text', 'NumHeaderLines',2);
        Pk  = (R_act * ([A(:,1), A(:,2), A(:,3)-ZC*1e3] * 1e3).').';   % um, actuator
        Bk  = (R_act * (A(:,4:6) * 1e3).').';                          % mT, actuator
        sel = vecnorm(Pk,2,2) <= R;
        if k == 1
            P  = Pk(sel,:);   Np = size(P,1);   B6 = zeros(Np,3,6);
        else
            assert(max(abs(Pk(sel,:) - P), [], 'all') < 1e-6, ...
                   'coil %d sits on a different export grid', k);
        end
        B6(:,:,k) = Bk(sel,:);
        clear A Pk Bk
    end
    C = struct('P',P, 'B6',B6, 'R',R, 'map',map, 'R_act',R_act, ...
               'dact',R_act*dhat, 'ZC',ZC);
    fprintf('rbf_field: read %d nodes (r <= %g um) from the six exports\n', Np, R);
end
