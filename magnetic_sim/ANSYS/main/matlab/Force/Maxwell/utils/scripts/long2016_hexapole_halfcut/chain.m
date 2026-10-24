%% chain.m -- position continuation of I_opt_hat inside the LS fitting cube, one direction at a time.
%  Problem (as nn_label.m / rho_solve.m):  I_opt_hat = argmin ||I_hat||^2  s.t.  I_hat' Q_k I_hat = r_k,  k = 1..3,
%  Q_k = F H_I' L_k(pbar) F H_I,  r_hat = [cos(phi) cos(theta), cos(phi) sin(theta), sin(phi)].
%  Positions: the 30 um grid nodes of the cube |x|, |y|, |z| <= A (A = 45 um: 90 um cube, all octants, N_p = 27).
%  For every direction j (ls_fit.m order):
%    centre pbar = 0 : global minimum = the label of nn_label_sgn.mat (multistart + Gurobi)
%    step 1  along x : (+-1,0,0), (+-2,0,0), ...   start = the node one step closer to the centre along x
%    step 2  along y : (a,+-1,0), (a,+-2,0), ...   start = (a, b -+ 1, 0)
%    step 3  along z : (a,b,+-1), (a,b,+-2), ...   start = (a, b, c -+ 1)
%    i.e. the start of every node is its parent: one step back along z if c ~= 0, else along y, else along x;
%    fmincon (Lagrange conditions) from the parent's solution only, no multistart, no Gurobi.
%  LBL = 'd15': positions and labels of nn_label_d15.mat (15 um cube grid, 343 points, lattice step 15 um)
%  instead of nn_label_sgn.mat / nn_grid.mat -> chain_d15.mat.
%  Output: utils/data/long2016_hexapole_halfcut/chain.mat
%    phi_j, theta_j  N_j x 1          [deg]
%    pbar            N_p x 3                 positions (normalised), solve order
%    par             N_p x 1                 parent of each position (0 for the centre)
%    I_chain         6 x N_j x N_p    [A]    continued solution
%    I_lab           6 x N_j x N_p    [A]    labels (global minimum per position)
%    ok              N_j x N_p               fmincon converged and feasible
%    cos_par         N_j x N_p               cos(I_chain, I_chain at the parent)   (signed; 1 at the centre)
%    ratio           N_j x N_p               ||I_chain||^2 / ||I_lab||^2   (1 = global minimum reached)
if ~exist('A','var'), A = 45; end
if ~exist('LBL','var'), LBL = ''; end
clearvars -except A LBL;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
DAT = fullfile(CALROOT, 'utils', 'data', MODEL);
if isempty(LBL)
    S = load(fullfile(DAT, 'nn_label_sgn.mat'), 'P_act', 'ph', 'th', 'IOPT', 'l_hat', 'MATNAME');
    G = load(fullfile(DAT, 'nn_grid.mat'), 'H');  OUT = 'chain.mat';
else
    S = load(fullfile(DAT, ['nn_label_' LBL '.mat']), 'P_act', 'ph', 'th', 'IOPT', 'l_hat', 'MATNAME', 'DP');
    G = struct('H', S.DP);  OUT = ['chain_' LBL '.mat'];
end
l_hat = S.l_hat;

% ---- positions in solve order and their parents ----
n  = find(all(abs(S.P_act) <= A + 1e-9, 2));
K  = round(S.P_act(n,:) / G.H);                                % lattice index (a, b, c)
st = ones(numel(n), 1);  st(K(:,2) ~= 0) = 2;  st(K(:,3) ~= 0) = 3;  st(all(K == 0, 2)) = 0;   % 0 centre, 1 x-axis, 2 z = 0 plane, 3 rest
d  = abs(K(sub2ind(size(K), (1:numel(n)).', max(st, 1))));     % distance along the axis being extended
[~, o] = sortrows([st, d, K]);  n = n(o);  K = K(o,:);  st = st(o);
N_p = numel(n);  par = zeros(N_p, 1);
for p = 2:N_p
    Kp = K(p,:);  ax = st(p);                                  % step back along the axis of this node's step
    Kp(ax) = Kp(ax) - sign(Kp(ax));
    par(p) = find(all(K == Kp, 2));
    assert(par(p) < p, 'chain:order', 'parent solved after the node');
end
pbar = S.P_act(n,:) / l_hat;

% ---- Q_k at every position (force model as nn_label.m) ----
mgB = 0.0451;  UF = 1e3;                                       % M g_B [A um^2/mT]; unit factor -> pN
cal = load_flux_calib(MODEL, GEOM, S.MATNAME, 'current');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;                     % F g_I [pN/A^2]
FH  = sqrt(gF) * cal.Mbar_m;                                   % F H_I [sqrt(pN)/A]
L   = build_L(pbar * l_hat, cal.l_m, cal.e_m, cal.Pc_base);    % N_p x 6 x 6 x 3
Q   = zeros(6, 6, 3, N_p);
for p = 1:N_p
    for k = 1:3, Q_k = FH.' * squeeze(L(p,:,:,k)) * FH;  Q(:,:,k,p) = (Q_k + Q_k.')/2; end
end

% ---- directions (ls_fit.m order) and labels ----
[PH, TH] = ndgrid(S.ph, S.th);
U = TH < 360 & (abs(PH) < 90 | TH == 0);
phi_j = PH(U);  theta_j = TH(U);  N_j = numel(phi_j);
I_lab = reshape(S.IOPT(:,:,:,n), 6, [], N_p);  I_lab = I_lab(:, U(:), :);
assert(all(K(1,:) == 0), 'chain:centre', 'first node must be the centre');

TOLC = 1e-6;
opt = optimoptions('fmincon', 'Algorithm','sqp', 'SpecifyObjectiveGradient',true, ...
                   'SpecifyConstraintGradient',true, 'OptimalityTolerance',1e-10, 'ConstraintTolerance',1e-12, ...
                   'StepTolerance',1e-14, 'MaxIterations',1000, 'Display','off');
I_chain = nan(6, N_j, N_p);  ok = false(N_j, N_p);
I_chain(:,:,1) = I_lab(:,:,1);  ok(:,1) = true;
t0 = tic;
for j = 1:N_j
    r = [cosd(phi_j(j))*cosd(theta_j(j)); cosd(phi_j(j))*sind(theta_j(j)); sind(phi_j(j))];
    for p = 2:N_p
        Qp = Q(:,:,:,p);
        [I, ~, ef] = fmincon(@obj_, I_chain(:,j,par(p)), [],[],[],[],[],[], @(x) con_(x, Qp, r), opt);
        ok(j,p) = ef > 0 && norm(fq_(I, Qp) - r) < TOLC;
        I_chain(:,j,p) = I;
    end
    if mod(j, 250) == 0, fprintf('  %d / %d directions, %.0f s\n', j, N_j, toc(t0)); end
end

% ---- continuity and optimality ----
Ip = I_chain(:,:,max(par, 1));
cos_par = squeeze(sum(I_chain .* Ip, 1) ./ (vecnorm(I_chain, 2, 1) .* vecnorm(Ip, 2, 1)));
ratio   = squeeze(vecnorm(I_chain, 2, 1).^2 ./ vecnorm(I_lab, 2, 1).^2);

save(fullfile(DAT, OUT), 'phi_j', 'theta_j', 'pbar', 'par', 'I_chain', 'I_lab', 'ok', 'cos_par', 'ratio', 'A', 'N_p', 'N_j');
fprintf('N_p = %d, N_j = %d | converged %d / %d | cos_par < 0.9: %d | ratio > 1.001: %d (max %.3f) -> %s\n', ...
        N_p, N_j, nnz(ok), numel(ok), nnz(cos_par < 0.9), nnz(ratio > 1.001), max(ratio(:)), fullfile(DAT, OUT));

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
