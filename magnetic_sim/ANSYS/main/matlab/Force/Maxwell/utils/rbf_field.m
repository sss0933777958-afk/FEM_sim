function [b, grad, info] = rbf_field(R, kern)
%RBF_FIELD  Radial-basis interpolant of the six FEM field exports.
% =========================================================================
%   One call builds the interpolant and hands back the two evaluators, so the
%   expensive part -- one Phi, one factorisation, 18 right-hand sides -- happens
%   once no matter how many points are asked for afterwards.
%
%       [b, grad] = rbf_field(R, kern)
%
%   INPUT 1  資料範圍  R
%     scalar   [um]   the nodes with |p| <= R become the RBF centres
%     [r1 r2]  [um]   an annulus instead, r1 <= |p| <= r2
%
%   INPUT 2  核函數  kern   (struct; a bare scalar is read as rho with lam = 0)
%     .rho     [um]   shape parameter
%     .lam            ridge:  (Phi + lam*I) W = B.   lam = 0 is exact interpolation
%     .name           'gauss' (default) | 'mq' | 'imq'
%     .phi            optional handle phi(r2)            -- overrides .name
%     .dphi_r         optional handle phi'(r)/r of r2    -- required with .phi
%
%       gauss   phi = exp(-r^2/rho^2)        phi'(r)/r = -(2/rho^2) phi
%       mq      phi = sqrt(r^2 + rho^2)      phi'(r)/r = 1/sqrt(r^2 + rho^2)
%       imq     phi = 1/sqrt(r^2 + rho^2)    phi'(r)/r = -(r^2 + rho^2)^(-3/2)
%
%     Every kernel is parametrised by phi'(r)/r rather than phi'(r), which is
%     what removes the 0/0 at r = 0 -- a query point may sit exactly on a node.
%
%   OUTPUT 1  b(p)     -- a HANDLE, so it reads the way it is written
%     b(Pq)               -> Nq x 3      [mT]   Pq in [um], ACTUATOR frame
%     b(Pq, I)            -> Nq x 3             for the 6-coil current vector I [A]
%
%   OUTPUT 2  那點的梯度  -- also a handle
%     grad(Pq)            -> Nq x 3      [mT^2/um]   column j = d(b.b)/dx_j
%     grad(Pq, I)         -> Nq x 3                  for current vector I
%     grad(Pq, I, 'jac')  -> Nq x 3 x 3  [mT/um]     J(:,c,j) = db_c/dx_j
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
%     Optional: asking for it costs one extra Phi*W and one rcond, so a caller
%     that only wants the two evaluators pays nothing.  Use it when you need
%     something the evaluators do not cover.
%       P      Np x 3       node positions [um], actuator frame
%       B6     Np x 3 x 6   the raw FEM field at those nodes [mT], per coil
%       W6     Np x 3 x 6   the solved weights, per coil
%       Np                  how many nodes the range actually selected
%       resid               max|Phi W - B| over all six coils [mT]
%       wmax                max|W|
%       rcnd                rcond(Phi + lam I): how close Phi is to singular
%       map    1 x 6        export -> paper pole; identity on the Maxwell branch
%       R_act  3 x 3        the rotation used, to get back to the measure frame
%       R, kern             echoed so a saved .mat is self-describing
%     resid must equal lam*wmax exactly -- that identity is the arithmetic check
%     on the whole solve, and it is how the nodal residual was verified to be
%     -lam*W rather than anything data-dependent.  B6 and P are what you score
%     the model against at its own nodes; rcnd is what tells you the Gaussian
%     has gone numerically singular, which it does around rho = 120 um.
%
%   Frame and units are the study's convention throughout: positions in um and
%   fields in mT, both rotated into the actuator frame by R_act and translated
%   so the six pole tips are centred on the origin.
%
%   Example
%     [b, grad] = rbf_field(250, struct('rho',120, 'lam',1e-4));
%     xa = linspace(-150,150,1001).';   Pq = [xa, zeros(1001,2)];
%     U  = sum(b(Pq).^2, 2);                        % b.b        [mT^2]
%     dU = grad(Pq);   F = 0.5*0.0451*1e3*dU(:,1);  % F_xa       [pN]
%
%   [ADDED 2026-09-11] merged from temp_code/scripts/{rbf_w.m, rbf_w6.m,
%   rbf_bb.m}: the weight solve, the evaluation and the analytic gradient were
%   three separate scripts with the Gaussian hard-coded in each.
% =========================================================================

    % ---- kernel ---------------------------------------------------------
    if isnumeric(kern), kern = struct('rho', kern, 'lam', 0); end
    if ~isfield(kern,'lam')  || isempty(kern.lam),  kern.lam  = 0;       end
    if ~isfield(kern,'name') || isempty(kern.name), kern.name = 'gauss'; end
    assert(isfield(kern,'rho') && isscalar(kern.rho) && kern.rho > 0, ...
           'kern.rho [um] is required and must be positive');
    assert(isscalar(kern.lam) && kern.lam >= 0, 'kern.lam must be a non-negative scalar');

    if isfield(kern,'phi') && ~isempty(kern.phi)
        assert(isfield(kern,'dphi_r') && ~isempty(kern.dphi_r), ...
               'a custom kern.phi also needs kern.dphi_r = phi''(r)/r as a function of r^2');
        phi = kern.phi;   dphi_r = kern.dphi_r;   kern.name = 'custom';
    else
        c2 = kern.rho^2;
        switch lower(kern.name)
        case 'gauss'
            phi    = @(r2) exp(-r2 / c2);
            dphi_r = @(r2) -(2/c2) * exp(-r2 / c2);
        case 'mq'
            phi    = @(r2) sqrt(r2 + c2);
            dphi_r = @(r2) 1 ./ sqrt(r2 + c2);
        case 'imq'
            phi    = @(r2) 1 ./ sqrt(r2 + c2);
            dphi_r = @(r2) -(r2 + c2).^(-1.5);
        otherwise
            error('unknown kernel ''%s'' -- use gauss | mq | imq, or supply phi + dphi_r', ...
                  kern.name);
        end
    end

    % ---- nodes (cached across calls; the .fld read is the slow part) -----
    persistent C
    assert(isnumeric(R) && ~isempty(R) && numel(R) <= 2 && all(R > 0), ...
           'R must be a positive scalar or a two-element [r1 r2], in um');
    RMAX = max(R(:));
    if isempty(C) || C.R < RMAX
        C = read_exports_(max(RMAX, 250));
    end

    rn = vecnorm(C.P, 2, 2);
    if isscalar(R), sel = rn <= R;
    else,           sel = rn >= min(R) & rn <= max(R);
    end
    P  = C.P(sel,:);   B6 = C.B6(sel,:,:);   Np = size(P,1);
    assert(Np > 0, 'no nodes fall inside the requested data range');

    % ---- one Phi, one factorisation, 18 right-hand sides ----------------
    s2  = sum(P.^2, 2);
    D2n = max(s2 + s2.' - 2*(P*P.'), 0);   D2n(1:Np+1:end) = 0;   % exact zero diagonal
    Phi = phi(D2n);   clear D2n s2
    W6  = reshape((Phi + kern.lam*eye(Np)) \ reshape(B6, Np, 18), Np, 3, 6);

    % ---- evaluators ------------------------------------------------------
    b    = @(varargin) eval_b_(P, W6, phi,         varargin{:});
    grad = @(varargin) eval_g_(P, W6, phi, dphi_r, varargin{:});

    if nargout < 3, return; end
    resid = 0;
    for k = 1:6
        resid = max(resid, max(abs(Phi*W6(:,:,k) - B6(:,:,k)), [], 'all'));
    end
    info = struct('P',P, 'B6',B6, 'W6',W6, 'Np',Np, 'R',R, 'kern',kern, ...
                  'map',C.map, 'R_act',C.R_act, 'resid',resid, ...
                  'wmax',max(abs(W6),[],'all'), 'rcnd',rcond(Phi + kern.lam*eye(Np)));
end

% =========================================================================
function bq = eval_b_(P, W6, phi, Pq, I)
    if nargin < 5, I = []; end
    W  = wsum_(W6, I);
    bq = phi(dist2_(Pq, P)) * W;
end

function out = eval_g_(P, W6, phi, dphi_r, Pq, I, what)
    if nargin < 6, I    = [];   end
    if nargin < 7, what = 'du'; end
    W  = wsum_(W6, I);
    D2 = dist2_(Pq, P);
    K  = dphi_r(D2);                                  % phi'(r)/r, finite at r = 0
    Nq = size(Pq,1);
    J  = zeros(Nq, 3, 3);                             % J(:,c,j) = db_c/dx_j
    for j = 1:3
        J(:,:,j) = (K .* (Pq(:,j) - P(:,j).')) * W;
    end
    if strcmpi(what, 'jac'), out = J;  return;  end
    assert(strcmpi(what,'du'), 'the third argument must be ''du'' (default) or ''jac''');
    bq  = phi(D2) * W;                                % d(b.b)/dx_j = 2 sum_c b_c db_c/dx_j
    out = zeros(Nq, 3);
    for j = 1:3
        out(:,j) = 2 * sum(bq .* J(:,:,j), 2);
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

function C = read_exports_(R)
% Read the six .fld exports once, translate to the tip-sphere centre and rotate
% into the actuator frame.  Lifted from rbf_w6.m; the field is kept RAW, i.e. the
% physical field for +1 A in the deck's own current direction.  The all-source
% sign flip of charge-model-source-convention.md is a presentation choice for the
% charge model and must NOT be applied here, because superposing
% b = sum_k I_k b_k needs the real per-ampere fields.
    HERE = fileparts(mfilename('fullpath'));                  % .../Maxwell/utils
    FMX  = fileparts(HERE);                                   % .../matlab/Force/Maxwell
    MAIN = fileparts(fileparts(fileparts(FMX)));              % .../main
    CAL  = fullfile(MAIN, 'matlab', 'Flux', 'Maxwell');
    addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'), fullfile(CAL,'common_path'));

    cfg = model_config('long2016_hexapole_halfcut', 'tip40um');
    ZC  = cfg.SPH_OFST;
    % [MODIFIED 2026-09-11] the [1 3 6 5 2 4] map is the APDL deck's build order and does
    % NOT apply to Maxwell: this config declares identity, i.e. export B_p<k> IS pole P<k>
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
    end
    C = struct('P',P, 'B6',B6, 'R',R, 'map',map, 'R_act',R_act, ...
               'dact',R_act*dhat, 'ZC',ZC);
    fprintf('rbf_field: cached %d nodes (r <= %g um) from the six exports\n', Np, R);
end
