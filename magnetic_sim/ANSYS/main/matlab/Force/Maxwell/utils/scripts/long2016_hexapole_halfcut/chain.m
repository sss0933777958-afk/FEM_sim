%% chain.m -- continuation test at the centre pbar = 0: does chaining the Lagrange solution from one direction to the
%  next remove the jumps of I_opt_hat between neighbouring directions?
%  Problem (as nn_label.m / rho_solve.m):  I_opt_hat = argmin ||I_hat||^2  s.t.  I_hat' Q_k I_hat = r_k,  k = 1..3,
%  Q_k = F H_I' L_k(pbar = 0) F H_I,  r_hat = [cos(phi) cos(theta), cos(phi) sin(theta), sin(phi)].
%  Directions j = 1..N_j in the ls_fit.m order (phi fastest, then theta; phi = -90 first, i.e. towards P6).
%    j = 1   : Lagrange conditions by fmincon multistart (NS random starts), then the Gurobi global check
%              (as rho_solve.m stage 2; the polished Gurobi point replaces it if smaller)
%    j >= 2  : fmincon from the solution of direction j-1 only (no multistart, no Gurobi)
%  Compared with the labels nn_label_sgn.mat at pbar = 0 (multistart + Gurobi per direction).
%  Output: utils/data/long2016_hexapole_halfcut/chain.mat
%    phi_j, theta_j  N_j x 1        [deg]
%    I_chain         6 x N_j        [A]    chained solution
%    I_lab           6 x N_j        [A]    label at pbar = 0
%    ok              N_j x 1               fmincon converged and feasible
%    cos_step        N_j-1 x 1             cos(I_chain(:,j), I_chain(:,j+1))   (signed)
%    ratio           N_j x 1               ||I_chain||^2 / ||I_lab||^2   (1 = global minimum reached)
%    GSTAT           Gurobi status at j = 1
if ~exist('NS','var'), NS = 40; end
if ~exist('TLIM','var'), TLIM = 600; end
clearvars -except NS TLIM;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
DAT = fullfile(CALROOT, 'utils', 'data', MODEL);
S   = load(fullfile(DAT, 'nn_label_sgn.mat'), 'P_act', 'ph', 'th', 'IOPT', 'l_hat', 'MATNAME');

% ---- Q_k at pbar = 0 (force model as nn_label.m) ----
mgB = 0.0451;  UF = 1e3;                                       % M g_B [A um^2/mT]; unit factor -> pN
cal = load_flux_calib(MODEL, GEOM, S.MATNAME, 'current');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;                     % F g_I [pN/A^2]
FH  = sqrt(gF) * cal.Mbar_m;                                   % F H_I [sqrt(pN)/A]
L   = build_L([0 0 0], cal.l_m, cal.e_m, cal.Pc_base);         % 1 x 6 x 6 x 3
Q   = zeros(6, 6, 3);
for k = 1:3, Q_k = FH.' * squeeze(L(1,:,:,k)) * FH;  Q(:,:,k) = (Q_k + Q_k.')/2; end

% ---- directions (ls_fit.m order) and labels at pbar = 0 ----
[PH, TH] = ndgrid(S.ph, S.th);
U = TH < 360 & (abs(PH) < 90 | TH == 0);
phi_j = PH(U);  theta_j = TH(U);  N_j = numel(phi_j);
p0 = find(all(S.P_act == 0, 2));
I_lab = reshape(S.IOPT(:,:,:,p0), 6, []);  I_lab = I_lab(:, U(:));

TOLC = 1e-6;
opt = optimoptions('fmincon', 'Algorithm','sqp', 'SpecifyObjectiveGradient',true, ...
                   'SpecifyConstraintGradient',true, 'OptimalityTolerance',1e-10, 'ConstraintTolerance',1e-12, ...
                   'StepTolerance',1e-14, 'MaxIterations',1000, 'Display','off');
r_hat = @(j) [cosd(phi_j(j))*cosd(theta_j(j)); cosd(phi_j(j))*sind(theta_j(j)); sind(phi_j(j))];
I_chain = nan(6, N_j);  ok = false(N_j, 1);

% ---- j = 1: multistart Lagrange + Gurobi check ----
rng(0);
r = r_hat(1);  R0 = randn(6, NS);  X0 = [[0 0 0 0 0 1].', 0.35 * R0 ./ vecnorm(R0)];
best = inf;
for s = 1:size(X0, 2)
    [I, ~, ef] = fmincon(@obj_, X0(:,s), [],[],[],[],[],[], @(x) con_(x, Q, r), opt);
    if ef > 0 && norm(fq_(I, Q) - r) < TOLC && I.'*I < best, best = I.'*I;  I_chain(:,1) = I; end
end
assert(isfinite(best), 'chain:j1', 'no feasible start at j = 1');
I0 = I_chain(:,1);  B = norm(I0);
gp = struct('NonConvex',2, 'MIPGap',1e-9, 'MIPGapAbs',1e-12, 'FeasibilityTol',1e-9, 'OptimalityTol',1e-9, ...
            'NumericFocus',2, 'TimeLimit',TLIM, 'Seed',0, 'OutputFlag',0);
m = struct('Q',speye(6), 'obj',zeros(6,1), 'A',sparse(0,6), 'rhs',zeros(0,1), 'sense','', ...
           'lb',-B*ones(6,1), 'ub',B*ones(6,1), 'start',I0);
for k = 1:3, m.quadcon(k) = struct('Qc',sparse(Q(:,:,k)), 'q',zeros(6,1), 'rhs',r(k), 'sense','='); end
m.quadcon(4) = struct('Qc',speye(6), 'q',zeros(6,1), 'rhs',B^2, 'sense','<');      % ball bound
res = gurobi(m, gp);  GSTAT = string(res.status);
if isfield(res, 'x')
    [I, ~, ef] = fmincon(@obj_, res.x, [],[],[],[],[],[], @(x) con_(x, Q, r), opt);   % polish
    if ef > 0 && norm(fq_(I, Q) - r) < TOLC && I.'*I < I0.'*I0, I_chain(:,1) = I; end
end
ok(1) = true;

% ---- j >= 2: fmincon from the previous direction's solution ----
for j = 2:N_j
    r = r_hat(j);
    [I, ~, ef] = fmincon(@obj_, I_chain(:,j-1), [],[],[],[],[],[], @(x) con_(x, Q, r), opt);
    ok(j) = ef > 0 && norm(fq_(I, Q) - r) < TOLC;
    I_chain(:,j) = I;
end

% ---- jumps and optimality ----
cos_step = (sum(I_chain(:,1:end-1) .* I_chain(:,2:end), 1) ./ (vecnorm(I_chain(:,1:end-1)) .* vecnorm(I_chain(:,2:end)))).';
ratio    = (vecnorm(I_chain).^2 ./ vecnorm(I_lab).^2).';

save(fullfile(DAT, 'chain.mat'), 'phi_j', 'theta_j', 'I_chain', 'I_lab', 'ok', 'cos_step', 'ratio', 'GSTAT', 'NS', 'TLIM');
fprintf('j = 1 (%g, %g): Gurobi %s | converged %d / %d | ratio: min %.6f, max %.4f, > 1.001 at %d directions\n', ...
        phi_j(1), theta_j(1), GSTAT, nnz(ok), N_j, min(ratio), max(ratio), nnz(ratio > 1.001));
fprintf('cos_step: min %.3f, < 0.9 at %d of %d steps -> %s\n', min(cos_step), nnz(cos_step < 0.9), N_j - 1, ...
        fullfile(DAT, 'chain.mat'));

% ---------------------------------------------------------------------------------------
function [v, g] = obj_(x)
% ||I_hat||^2 and its gradient
    v = x.'*x;  g = 2*x;
end

function [c, ceq, gc, gceq] = con_(x, Q, r)
% equality constraints I_hat' Q_k I_hat - r_k = 0 (k = 1..3) and their gradients (6 x 3)
    c = [];  gc = [];
    ceq  = [x.'*Q(:,:,1)*x - r(1); x.'*Q(:,:,2)*x - r(2); x.'*Q(:,:,3)*x - r(3)];
    gceq = 2*[Q(:,:,1)*x, Q(:,:,2)*x, Q(:,:,3)*x];
end

function f = fq_(x, Q)
% force fhat_k = I_hat' Q_k I_hat
    f = [x.'*Q(:,:,1)*x; x.'*Q(:,:,2)*x; x.'*Q(:,:,3)*x];
end
