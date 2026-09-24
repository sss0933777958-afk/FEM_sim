function [b, grad, info] = sph_field(R, FEM, opts)
%SPH_FIELD  Solid-harmonic (spherical-harmonic) model of a source-free FEM field.
% =========================================================================
%   The working region carries no current, so the field derives from a scalar
%   potential that satisfies Laplace's equation:
%
%       Phi(p) = sum_{k=1..K} c_k p_k(x,y,z),        b = grad(Phi)
%
%   with p_k the REGULAR SOLID HARMONICS -- r^l P_l^m(cos th) {cos m ph, sin m ph}
%   written out as homogeneous polynomials of degree l in x, y, z.  Counting
%   2l+1 terms per degree from l = 1 to L,
%
%       K = (L+1)^2 - 1
%
%   (l = 0 is the constant potential; its gradient is zero, so it is left out).
%   Only the c_k are unknown.  They come from a least-squares fit to the FEM
%   field: every node contributes three rows, b_c = sum_k c_k dp_k/dx_c, so the
%   system is 3*Np x K and heavily overdetermined.
%
%   WHY THIS BASIS.  Every p_k is harmonic by construction, so div b = 0 and
%   curl b = 0 hold identically and the model cannot invent a ripple: its
%   derivatives are polynomials, obtained exactly rather than by differencing
%   or by tuning a smoothing parameter.  Unlike an RBF there is no shape
%   parameter, no ridge, and no node-cloud edge to collapse against -- the only
%   choice is L.  The expansion is valid inside a source-free ball about the
%   origin; accuracy near the sources (pole tips) is what sets the L needed.
%
%   MODEL-AGNOSTIC.  Nothing here knows which magnet it is looking at: the
%   caller passes the data, so the same function serves every model.
%
%       [b, grad, info] = SPH_FIELD(R, FEM)
%       [b, grad, info] = SPH_FIELD(R, FEM, opts)
%
%   INPUT 1  R
%     scalar   [um]   fit the nodes with |p| <= R
%     [r1 r2]  [um]   an annulus instead, r1 <= |p| <= r2
%
%   INPUT 2  FEM   the data, as a struct
%     .P     Np x 3        node positions [um], in the frame the answers are wanted in
%     .B     Np x 3 x NC   field at those nodes [mT], one page per excitation
%                          (.B6 is accepted as an alias)
%     The frame, the units and the excitation order are the caller's business;
%     they come back untouched in b and grad.
%
%   INPUT 3  opts  (optional struct)
%     .L        maximum degree (default 8, i.e. K = 80)
%     .Rn       normalising radius [um] for the polynomials (default: max |p|
%               of the fitted nodes).  Coordinates enter as p/Rn, which keeps
%               the r^l columns from spanning orders of magnitude.
%     .quiet    true to suppress the one-line summary (default false)
%
%   OUTPUT 1  b(Pq[, I])      -> Nq x 3  field at Pq [same units as FEM.B]
%   OUTPUT 2  grad(Pq[, I[, what]])  -- by default grad(B) itself, i.e. the three
%     derivative vectors written directly from the second derivatives of p_k:
%
%         dB/dx = sum_k c_k [d_xx p_k; d_yx p_k; d_zx p_k]
%         dB/dy = sum_k c_k [d_xy p_k; d_yy p_k; d_zy p_k]
%         dB/dz = sum_k c_k [d_xz p_k; d_yz p_k; d_zz p_k]
%
%       what = 'jac' (default) -> Nq x 3 x 3.  G(:,:,j) is dB/dx_j, so
%                                 G(:,c,j) = dB_c/dx_j = sum_k c_k d_{c j} p_k.
%                                 It is symmetric in (c,j) -- both equal
%                                 d^2 Phi/dx_c dx_j -- which is curl B = 0.
%              'du'            -> Nq x 3      column j = d(b.b)/dx_j   (the force)
%              'd2'            -> Nq x 3      column j = d2(b.b)/dx_j^2
%              'd3'            -> Nq x 3      column j = d3(b.b)/dx_j^3
%     I is the excitation vector (NC x 1); it defaults to the first excitation.
%     The FIELD superposes, so one set of coefficients per excitation is enough
%     and c(I) = sum_k I_k c_k; the FORCE does not, being quadratic in b.
%     WARNING  the signatures match RBF_FIELD but the DEFAULTS do not: rbf_field
%     returns d(b.b)/dx_j when `what` is omitted, this one returns the Jacobian.
%     Always pass `what` explicitly when comparing the two, e.g.
%         dU = grad(Pq, I, 'du');   F_j = 0.5 * mgB * UF * dU(:,j)      [pN]
%
%   OUTPUT 3  info
%       c        K x NC    the fitted coefficients
%       L K R Rn Np NC     what was fitted
%       resid    Np x NC   || b_model - B_FEM || at the fitted nodes [mT]
%       NMAE     1 x NC    sum||resid|| / sum||B|| * 100, per excitation
%       NMAE_all scalar    the same over every excitation
%       cnd      condition number of the (column-scaled) design matrix
%       basis    the polynomial table, and secs / L / K for the record
%
%   HOW THE POLYNOMIALS ARE BUILT.  The complex combination
%       A_l^m = r^l P_l^m(cos th) exp(i m ph)
%   is a polynomial in x, y, z, and the standard Legendre recurrence
%       (l-m) P_l^m = (2l-1) u P_{l-1}^m - (l+m-1) P_{l-2}^m
%   becomes, after multiplying by r^l exp(i m ph) and using u*r = z,
%       (l-m) A_l^m = (2l-1) z A_{l-1}^m - (l+m-1) r^2 A_{l-2}^m
%   seeded by A_m^m = (2m-1)!! (x + i y)^m and A_{m+1}^m = (2m+1) z A_m^m.
%   Every step is an exact polynomial operation (multiply by z, multiply by
%   r^2 = x^2+y^2+z^2), so the monomial table is built with no quadrature and
%   no special functions, and it is exact on the axis and at the origin -- the
%   places where the spherical form itself is singular.  Re and Im of A_l^m are
%   the cos(m ph) and sin(m ph) members.
%
%   See also RBF_FIELD.

    t0 = tic;
    if nargin < 3 || isempty(opts), opts = struct(); end
    L     = getd(opts, 'L', 8);
    quiet = getd(opts, 'quiet', false);

    validateattributes(L, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'opts.L');
    validateattributes(R, {'numeric'}, {'vector','real','positive','finite'}, mfilename, 'R');
    assert(numel(R) <= 2, 'sph_field:R', 'R is a scalar radius or [r1 r2]');

    % ---- data ------------------------------------------------------------------
    assert(isstruct(FEM) && isfield(FEM,'P'), 'sph_field:FEM', ...
           'FEM must be a struct with .P (Np x 3) and .B (Np x 3 x NC)');
    if     isfield(FEM,'B'),  Ball = FEM.B;
    elseif isfield(FEM,'B6'), Ball = FEM.B6;
    else,  error('sph_field:FEM', 'FEM needs a .B (or .B6) field');
    end
    Pall = FEM.P;
    validateattributes(Pall, {'numeric'}, {'2d','ncols',3,'real','finite'}, mfilename, 'FEM.P');
    assert(size(Ball,1) == size(Pall,1) && size(Ball,2) == 3, 'sph_field:size', ...
           'FEM.B must be Np x 3 x NC to match FEM.P');

    rr = vecnorm(Pall, 2, 2);
    if numel(R) == 1, keep = rr <= R;  else, keep = rr >= R(1) & rr <= R(2); end
    assert(any(keep), 'sph_field:empty', 'no node inside R');
    P  = Pall(keep,:);   Bd = Ball(keep,:,:);
    Np = size(P,1);      NC = size(Bd,3);

    Rn = getd(opts, 'Rn', max(vecnorm(P,2,2)));
    X  = P / Rn;                                   % normalised coordinates

    % ---- basis and its first derivatives at the nodes ---------------------------
    [basis, lab] = build_basis_(L);                % lab(k,:) = [l m t], t = 1 cos, 0 sin
    K     = numel(basis);
    assert(K == (L+1)^2 - 1, 'sph_field:K', 'basis count %d is not (L+1)^2-1', K);

    D = zeros(Np, K, 3);
    for j = 1:3
        a = zeros(1,3);  a(j) = 1;
        D(:,:,j) = eval_basis_(basis, X, a) / Rn;  % dp/dx_j in physical units
    end

    % ---- least squares  J*C = B  (all NC excitations share one factorisation) ---
    %   J  3Np x K   column k, row (node,component c) = dp_k/dx_c
    %   C  K   x NC  the unknown coefficients, one column per excitation
    %   B  3Np x NC  the FEM field, stacked in the same row order
    J = [D(:,:,1); D(:,:,2); D(:,:,3)];
    B = [reshape(Bd(:,1,:), Np, NC); reshape(Bd(:,2,:), Np, NC); reshape(Bd(:,3,:), Np, NC)];
    sc = vecnorm(J, 2, 1);   sc(sc == 0) = 1;      % column scaling, undone below
    Js = J ./ sc;
    C  = (Js \ B) ./ sc.';
    % [ADDED 2026-09-17] cond() is an SVD of a 3Np x K matrix and is purely
    % diagnostic -- nothing downstream reads info.cnd.  Its cost is ~2*(3Np)*K^2,
    % so on an uncapped degree sweep over every node it overtakes the fit itself:
    % at Np = 47,693 and L = 50 (K = 2600) it is minutes per degree.  opts.condmax
    % skips it above a given K and reports NaN.  Default Inf = always compute, so
    % every existing caller behaves exactly as before.
    if K <= getd(opts, 'condmax', Inf), cnd = cond(Js);  else, cnd = NaN;  end

    % ---- residual  r = J*C - B  -------------------------------------------------
    r       = J*C - B;                             % 3Np x NC [same units as B]
    rms_exc = sqrt(mean(r.^2, 1));                 % per excitation
    rms     = sqrt(mean(r.^2, 'all'));             % all six together

    Bm = zeros(Np, 3, NC);
    for j = 1:3, Bm(:,j,:) = reshape(D(:,:,j) * C, Np, 1, NC); end
    resid = squeeze(vecnorm(Bm - Bd, 2, 2));       % Np x NC, vector norm per node
    nB    = squeeze(vecnorm(Bd, 2, 2));
    NMAE  = sum(resid, 1) ./ sum(nB, 1) * 100;

    % ---- handles ----------------------------------------------------------------
    b    = @(varargin) eval_b_(basis, Rn, C, varargin{:});
    grad = @(varargin) eval_grad_(basis, Rn, C, varargin{:});

    info = struct('C', C, 'lab', lab, 'L', L, 'K', K, 'R', R, 'Rn', Rn, 'Np', Np, 'NC', NC, ...
                  'rms', rms, 'rms_exc', rms_exc, ...
                  'resid', resid, 'NMAE', NMAE, ...
                  'NMAE_all', sum(resid,'all')/sum(nB,'all')*100, ...
                  'cnd', cnd, 'basis', {basis}, 'secs', toc(t0));
    if ~quiet
        fprintf(['sph_field: L = %d (K = %d), %d nodes, %d excitations | ' ...
                 'rms(J*C-B) = %.5g | NMAE %.4f %% | cond %.2e | %.1f s' newline], ...
                L, K, Np, NC, rms, info.NMAE_all, cnd, info.secs);
    end
end

% ================================ evaluation =====================================
function B = eval_b_(basis, Rn, C, Pq, I)
% b_c = dPhi/dx_c = sum_k C(k) dp_k/dx_c
    if nargin < 5, I = []; end
    cv = coeff_(C, I);
    X  = Pq / Rn;
    B  = zeros(size(Pq,1), 3);
    for j = 1:3
        a = zeros(1,3);  a(j) = 1;
        B(:,j) = eval_basis_(basis, X, a) * cv / Rn;
    end
end

function G = eval_grad_(basis, Rn, C, Pq, I, what)
% grad(B) by default, and the derivatives of U = b.b on request.  b_c = dPhi/dx_c,
% so the n-th derivative of b_c along axis j is the (n+1)-th derivative of Phi;
% in particular dB_c/dx_j = sum_k c_k d_{cj} p_k, the form written in the header.
    if nargin < 5, I    = [];    end
    if nargin < 6, what = 'jac'; end
    cv = coeff_(C, I);
    X  = Pq / Rn;   Nq = size(Pq,1);

    dB = @(c_, j_, n_) eval_basis_(basis, X, unit_(c_) + n_*unit_(j_)) * cv / Rn^(n_+1);

    switch lower(what)
        case 'jac'                                   % J(:,c,j) = db_c/dx_j
            G = zeros(Nq,3,3);
            for cc = 1:3, for j = 1:3, G(:,cc,j) = dB(cc, j, 1); end, end
        case {'du','d2','d3'}
            nmax = 1;  if strcmpi(what,'d2'), nmax = 2; elseif strcmpi(what,'d3'), nmax = 3; end
            G = zeros(Nq,3);
            for j = 1:3
                b0 = zeros(Nq,3);  b1 = zeros(Nq,3);  b2 = zeros(Nq,3);  b3 = zeros(Nq,3);
                for cc = 1:3
                    b0(:,cc) = dB(cc, j, 0);
                    b1(:,cc) = dB(cc, j, 1);
                    if nmax >= 2, b2(:,cc) = dB(cc, j, 2); end
                    if nmax >= 3, b3(:,cc) = dB(cc, j, 3); end
                end
                switch nmax
                    case 1, G(:,j) = 2*sum(b0.*b1, 2);
                    case 2, G(:,j) = 2*sum(b1.^2 + b0.*b2, 2);
                    case 3, G(:,j) = 2*sum(3*b1.*b2 + b0.*b3, 2);
                end
            end
        otherwise
            error('sph_field:what', 'what must be ''du'', ''jac'', ''d2'' or ''d3''');
    end
end

function e = unit_(j), e = zeros(1,3);  e(j) = 1;  end

function cv = coeff_(C, I)
% One column of coefficients: the excitation picks a linear combination of the
% columns of C, which is legitimate because the FIELD superposes (the force,
% being quadratic in b, does not).
    NC = size(C,2);
    if isempty(I)
        cv = C(:,1);
    else
        I = I(:);
        assert(numel(I) == NC, 'sph_field:I', 'I needs %d entries, one per excitation', NC);
        cv = C * I;
    end
end

% ================================ the basis ======================================
function V = eval_basis_(basis, X, a)
% Np x K matrix of d^a p_k / dx^a at the (normalised) points X.
%   a = [ax ay az] derivative orders.  Monomials differentiate by the falling
%   factorial, and a term dies once the order exceeds its exponent.
    Np = size(X,1);   K = numel(basis);
    dmax = max(cellfun(@(bk) max(bk.e(:)), basis));            % highest exponent present
    PW = cell(1,3);
    for d = 1:3
        PW{d} = ones(Np, dmax+1);
        for q = 2:dmax+1, PW{d}(:,q) = PW{d}(:,q-1) .* X(:,d); end   % column q = x^(q-1)
    end
    V = zeros(Np, K);
    for k = 1:K
        e = basis{k}.e;   cf = basis{k}.c;
        v = zeros(Np,1);
        for t = 1:size(e,1)
            ex = e(t,:);
            if any(ex < a), continue, end                       % differentiated away
            f = cf(t);
            for d = 1:3
                f = f * prod(ex(d) - (0:a(d)-1));               % falling factorial
            end
            if f == 0, continue, end
            ex = ex - a;
            v = v + f * (PW{1}(:,ex(1)+1) .* PW{2}(:,ex(2)+1) .* PW{3}(:,ex(3)+1));
        end
        V(:,k) = v;
    end
end

function [basis, lab] = build_basis_(L)
% The regular solid harmonics up to degree L as monomial tables, in the order
%   l = 1..L, and within each l:  m = 0, then (cos, sin) for m = 1..l.
% A_l^m = r^l P_l^m(cos th) exp(i m ph) is built by the Legendre recurrence
% rewritten for polynomials (see the header), so only exact integer arithmetic
% on exponents is involved.
    basis = cell(1, (L+1)^2 - 1);
    lab   = zeros((L+1)^2 - 1, 3);                      % [l m t], t = 1 cos m ph, 0 sin m ph
    n = 0;
    for l = 1:L
        for m = 0:l
            [Re, Im] = solid_(l, m);
            n = n + 1;   basis{n} = Re;   lab(n,:) = [l m 1];   % m = 0 gives the real one only
            if m > 0, n = n + 1;  basis{n} = Im;  lab(n,:) = [l m 0];  end
        end
    end
    % reorder so that every degree keeps its 2l+1 members together (already true)
    assert(n == numel(basis), 'sph_field:count', 'built %d of %d basis members', n, numel(basis));
end

function [Re, Im] = solid_(l, m)
% Re / Im parts of A_l^m, as monomial tables.
    % seed: A_m^m = (2m-1)!! (x + i y)^m
    df = 1;  for q = 1:2:(2*m-1), df = df * q;  end            % (2m-1)!!
    eR = zeros(0,3);  cR = zeros(0,1);   eI = zeros(0,3);  cI = zeros(0,1);
    for t = 0:m
        cf = df * nchoosek(m,t);
        switch mod(t,4)
            case 0, eR(end+1,:) = [m-t t 0];  cR(end+1,1) =  cf;   %#ok<AGROW>
            case 1, eI(end+1,:) = [m-t t 0];  cI(end+1,1) =  cf;   %#ok<AGROW>
            case 2, eR(end+1,:) = [m-t t 0];  cR(end+1,1) = -cf;   %#ok<AGROW>
            case 3, eI(end+1,:) = [m-t t 0];  cI(end+1,1) = -cf;   %#ok<AGROW>
        end
    end
    Amm = {struct('e',eR,'c',cR), struct('e',eI,'c',cI)};      % {Re, Im} of A_m^m

    if l == m
        Re = tidy_(Amm{1});  Im = tidy_(Amm{2});  return
    end

    prev2 = Amm;                                               % A_{m}^m
    prev1 = {mulz_(Amm{1}, 2*m+1), mulz_(Amm{2}, 2*m+1)};      % A_{m+1}^m
    if l == m+1
        Re = tidy_(prev1{1});  Im = tidy_(prev1{2});  return
    end
    for ll = m+2:l
        cur = cell(1,2);
        for part = 1:2
            t1  = mulz_(prev1{part}, 2*ll-1);
            t2  = mulr2_(prev2{part}, -(ll+m-1));
            cur{part} = scale_(add_(t1, t2), 1/(ll-m));
        end
        prev2 = prev1;   prev1 = cur;
    end
    Re = tidy_(prev1{1});  Im = tidy_(prev1{2});
end

% ---- polynomial helpers (monomial tables: .e exponents, .c coefficients) --------
function q = mulz_(p, s),  q = struct('e', p.e + [0 0 1], 'c', s * p.c);  end

function q = mulr2_(p, s)                                     % multiply by r^2
    q = struct('e', [p.e + [2 0 0]; p.e + [0 2 0]; p.e + [0 0 2]], 'c', s * [p.c; p.c; p.c]);
end

function q = scale_(p, s), q = struct('e', p.e, 'c', s * p.c);  end

function q = add_(p1, p2), q = tidy_(struct('e', [p1.e; p2.e], 'c', [p1.c; p2.c]));  end

function q = tidy_(p)                                          % merge like terms
    if isempty(p.c), q = struct('e', zeros(0,3), 'c', zeros(0,1));  return, end
    [e, ~, ix] = unique(p.e, 'rows');
    c = accumarray(ix, p.c);
    k = abs(c) > 0;
    q = struct('e', e(k,:), 'c', c(k));
    if isempty(q.c), q = struct('e', zeros(1,3), 'c', 0);  end
end

function v = getd(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f);  else, v = d;  end
end
