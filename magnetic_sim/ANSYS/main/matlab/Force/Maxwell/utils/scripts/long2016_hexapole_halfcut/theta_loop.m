%% theta_loop.m -- one ring around theta at fixed phi, at the origin: does the walk come back to the same solution?
%  Problem as branch.m (branch_jump.pdf, 6-dim Phihat, Long's model: R_norm 550 um, e = 0):
%      min Phihat' Phihat   s.t.   Phihat' L_k Phihat = r_hat_k,   r_hat = [cos(phi) cos(theta); cos(phi) sin(theta); sin(phi)]
%  Start: the solution of the branch file BRF at (theta, phi) = (TH0, PH0); then theta += DA, every step fmincon (sqp)
%  from the previous solution (sign aligned to it), all the way round to theta = TH0 + 360.
%  Parameters: BRF ('branch_oct_h15_a45_lg_ring_single'), PH0 (0), TH0 (6), DA (6) [deg]
%  Output: utils/data/long2016_hexapole_halfcut/theta_loop.mat (TH0 = 6, PH0 = 0), else theta_loop_t<TH0>_p<PH0>.mat
%    th_s (1 x N, TH0 : DA : TH0+360), X (6 x N, the walked solutions), cost (1 x N, ||Phihat||^2), cs_step (|cos| to the
%    previous step), cs_br (|cos| to the BRF solution at the same theta, N-1 values: the last theta equals the first),
%    cs_back / cs_back_signed (|cos| / cos of the end to the start), ef (fmincon exit flags), BRF, PH0, TH0, DA
if ~exist('BRF','var'), BRF = 'branch_oct_h15_a45_lg_ring_single'; end
if ~exist('PH0','var'), PH0 = 0; end
if ~exist('TH0','var'), TH0 = 6; end
if ~exist('DA','var'), DA = 6; end
% TP: probe instead of the ring -- start from the solution of theta_loop.mat (the walk from TH0 = 6) at theta = TP(1),
%     then walk through the theta values TP (any spacing, e.g. up to and across 90 deg), each from the previous solution;
%     record cost, |cos| to the previous, the eigenvalues of N = Z' M Z, lambda and cond([L_x Phihat, L_y Phihat, L_z Phihat])
%     -> theta_loop_probe.mat (tp, Xp, costp, csp, eNp (3 x n), lamp (3 x n), condA, ef)
if ~exist('TP','var'), TP = []; end
clearvars -except BRF PH0 TH0 DA TP;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
DAT = fullfile(CALROOT, 'utils', 'data', MODEL);
cal = load_flux_calib(MODEL, GEOM, 'calib_current_maxwell_axshN43_R150_single', 'current');   % Pc_base only
L4  = build_L([0 0 0], 550, zeros(17,1), cal.Pc_base);  L = reshape(L4(1,:,:,:), 6, 6, 3);
for k = 1:3, L(:,:,k) = (L(:,:,k) + L(:,:,k).')/2; end

if ~isempty(TP)                                                % ---- probe along theta (see TP above) ----
    W = load(fullfile(DAT, 'theta_loop.mat'), 'X', 'th_s', 'PH0');
    xs = W.X(:, abs(W.th_s - TP(1)) < 1e-9);  assert(~isempty(xs), 'theta_loop:tp', 'TP(1) must be a theta of theta_loop.mat');
    opt = optimoptions('fmincon', 'Algorithm','sqp', 'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true, ...
                       'OptimalityTolerance',1e-12, 'ConstraintTolerance',1e-13, 'StepTolerance',1e-15, 'MaxIterations',3000, 'Display','off');
    tp = TP(:).';  n = numel(tp);  Xp = zeros(6, n);  [costp, csp, condA] = deal(nan(1, n));  [eNp, lamp] = deal(nan(3, n));  ef = ones(1, n);
    for q = 1:n
        r = [cosd(W.PH0)*cosd(tp(q)); cosd(W.PH0)*sind(tp(q)); sind(W.PH0)];
        if q > 1, [y, ~, ef(q)] = fmincon(@obj_, xs, [], [], [], [], [], [], @(x) con_(x, L, r), opt);  if y.'*xs < 0, y = -y; end
        else, y = xs; end
        A = [L(:,:,1)*y, L(:,:,2)*y, L(:,:,3)*y];  lam = A \ y;
        M = eye(6) - lam(1)*L(:,:,1) - lam(2)*L(:,:,2) - lam(3)*L(:,:,3);  M = (M + M.')/2;
        Z = null(A.');  N = Z.'*M*Z;
        Xp(:,q) = y;  costp(q) = y.'*y;  csp(q) = abs(y.'*xs)/(norm(y)*norm(xs));  eNp(:,q) = sort(eig((N + N.')/2));
        lamp(:,q) = lam;  condA(q) = cond(A);  xs = y;
        fprintf('%10.4f  cost %.6f  |cos| %.5f  eig N [%+.4f %+.4f %+.4f]  lambda [%+.4f %+.4f %+.4f]  cond %.3g\n', tp(q), costp(q), csp(q), eNp(:,q), lam, condA(q));
    end
    OUT = fullfile(DAT, 'theta_loop_probe.mat');
    save(OUT, 'tp', 'Xp', 'costp', 'csp', 'eNp', 'lamp', 'condA', 'ef');
    fprintf('-> %s\n', OUT);
    return
end
R  = load(fullfile(DAT, [BRF '.mat']), 'Phihat', 'ph', 'th', 'P_act');
i0 = find(vecnorm(R.P_act, 2, 2) < 1e-9, 1);  a0 = find(abs(R.ph - PH0) < 1e-9, 1);
x0 = R.Phihat(:, a0, abs(R.th - mod(TH0, 360)) < 1e-9, i0);
opt = optimoptions('fmincon', 'Algorithm','sqp', 'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true, ...
                   'OptimalityTolerance',1e-10, 'ConstraintTolerance',1e-12, 'StepTolerance',1e-14, 'MaxIterations',1000, 'Display','off');
th_s = TH0 : DA : TH0 + 360;  N = numel(th_s);
X = zeros(6, N);  X(:,1) = x0;  cost = zeros(1, N);  cost(1) = x0.'*x0;  cs_step = nan(1, N);  ef = ones(1, N);
for n = 2:N
    r = [cosd(PH0)*cosd(th_s(n)); cosd(PH0)*sind(th_s(n)); sind(PH0)];
    [x, ~, ef(n)] = fmincon(@obj_, X(:,n-1), [], [], [], [], [], [], @(x) con_(x, L, r), opt);
    if x.'*X(:,n-1) < 0, x = -x; end
    X(:,n) = x;  cost(n) = x.'*x;  cs_step(n) = abs(x.'*X(:,n-1))/(norm(x)*norm(X(:,n-1)));
end
cs_br = nan(1, N-1);
for n = 1:N-1
    y = R.Phihat(:, a0, abs(R.th - mod(th_s(n), 360)) < 1e-9, i0);  cs_br(n) = abs(X(:,n).'*y)/(norm(X(:,n))*norm(y));
end
xe = X(:,end);  cs_back = abs(xe.'*x0)/(norm(xe)*norm(x0));  cs_back_signed = xe.'*x0/(norm(xe)*norm(x0));
OUT = fullfile(DAT, 'theta_loop.mat');
if TH0 ~= 6 || PH0 ~= 0, OUT = fullfile(DAT, sprintf('theta_loop_t%g_p%g.mat', TH0, PH0)); end   % other starts: own file
save(OUT, 'th_s', 'X', 'cost', 'cs_step', 'cs_br', 'cs_back', 'cs_back_signed', 'ef', 'BRF', 'PH0', 'TH0', 'DA');
[m, k] = min(cs_step(2:end));
fprintf('phi = %g deg, start theta = %g deg (%s), %d steps of %g deg\n', PH0, TH0, BRF, N-1, DA);
fprintf('back at the start: |cos| %.6f (signed %+.6f), cost start %.6f, end %.6f\n', cs_back, cs_back_signed, cost(1), cost(end));
fprintf('step |cos|: min %.4f (theta %g -> %g deg), steps < 0.95: %d, fmincon not converged: %d\n', m, th_s(k), th_s(k+1), nnz(cs_step(2:end) < 0.95), nnz(ef(2:end) <= 0));
fprintf('-> %s\n', OUT);

function [v, g] = obj_(x)
v = x.'*x;  g = 2*x;
end

function [c, ceq, gc, gceq] = con_(x, L, r)
c = [];  gc = [];
ceq  = [x.'*L(:,:,1)*x - r(1); x.'*L(:,:,2)*x - r(2); x.'*L(:,:,3)*x - r(3)];
gceq = 2*[L(:,:,1)*x, L(:,:,2)*x, L(:,:,3)*x];
end
