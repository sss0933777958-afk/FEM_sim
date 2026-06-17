%% fit_B6x_6C.m — [B-6x] with per-pole C_i (6 amplitudes, not 1)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
%
%  Alternating superposition: I_apdl = [+1,-1,+1,-1,-1,+1]
%  w = I_diss = [-1,+1,-1,+1,+1,-1], sum=0
%  18 position params (fminsearch) + 6 per-pole C_i (VarPro linear LS)

cnst = mt_constants();
cube_half = 50e-6;

fprintf('=== [B-6x] with per-pole C_i ===\n\n');

% Alternating current
I_apdl = [1; -1; 1; -1; -1; 1];
I_diss = -I_apdl;
w = I_diss;  % sum=0 → K_I drops out
fprintf('w = [%+d, %+d, %+d, %+d, %+d, %+d]\n\n', w);

%% 1. Build alternating superposition FEM data
d1 = import_ansys_data(fullfile('..','..','results','coil1','filleted'), 'wp', 'coil1');
mask = abs(d1.x)<cube_half & abs(d1.y)<cube_half & abs(d1.z)<cube_half;
N = sum(mask);
px = d1.x(mask); py = d1.y(mask); pz = d1.z(mask);
bx_f = I_apdl(1)*d1.bx(mask); by_f = I_apdl(1)*d1.by(mask); bz_f = I_apdl(1)*d1.bz(mask);
for k = 2:6
    dk = import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'), ...
        'wp', sprintf('coil%d',k));
    bx_f = bx_f + I_apdl(k)*dk.bx(mask);
    by_f = by_f + I_apdl(k)*dk.by(mask);
    bz_f = bz_f + I_apdl(k)*dk.bz(mask);
end
b_fem = [bx_f; by_f; bz_f];
fprintf('Data: %d nodes in +/-%.0f um, |B| mean=%.3f mT\n\n', N, cube_half*1e6, ...
    mean(sqrt(bx_f.^2+by_f.^2+bz_f.^2))*1e3);

%% 2. Load initial positions from [J]
J = load('../../data/J_ideal_fit.mat');
pos_J = J.pos;

% Also build sphere positions from KI_fit
KI_data = load('../../data/KI_fit.mat');
pos_sphere = KI_data.pos;

alpha = cnst.alpha;
d_hat = zeros(3,6);
for i = 1:6
    theta = cnst.pole_angles(i)*pi/180;
    zs = sign(cnst.pole_tip_z_wp(i));
    d_hat(:,i) = [cos(theta)*sin(alpha); sin(theta)*sin(alpha); zs*cos(alpha)];
end

%% 3. Sphere baseline (per-pole C_i)
[cost0, C0] = cost_perC(pos_sphere(:), w, px, py, pz, b_fem, N);
fprintf('Sphere baseline: cost=%.6e\n', cost0);
fprintf('  C_i: '); fprintf('%.4e ', C0); fprintf('\n\n');

%% 4. Optimize
cost_fn = @(x) cost_perC(x, w, px, py, pz, b_fem, N);

opts = optimset('TolX', 1e-10, 'TolFun', 1e-22, ...
    'MaxIter', 50000, 'MaxFunEvals', 50000, 'Display', 'off');

x0_list = {};
x0_names = {};
x0_list{end+1} = pos_sphere(:);    x0_names{end+1} = 'sphere';
x0_list{end+1} = pos_J(:);         x0_names{end+1} = '[J] positions';
rng(42); x0_list{end+1} = pos_J(:) + 20e-6*randn(18,1); x0_names{end+1} = '[J]+20um';

n_trials = length(x0_list);
trials = struct([]);

for t = 1:n_trials
    fprintf('----- Trial %d: %s -----\n', t, x0_names{t});
    tic;
    [x_opt, fval, ~, output] = fminsearch(cost_fn, x0_list{t}, opts);
    elapsed = toc;
    [~, Ci] = cost_perC(x_opt, w, px, py, pz, b_fem, N);
    trials(t).name = x0_names{t};
    trials(t).pos = reshape(x_opt, 3, 6);
    trials(t).cost = fval;
    trials(t).C_i = Ci;
    trials(t).iters = output.iterations;
    trials(t).fevals = output.funcCount;
    trials(t).elapsed = elapsed;
    fprintf('  Cost=%.6e (iters %d, fevals %d, %.0fs)\n', fval, output.iterations, output.funcCount, elapsed);
end

%% 5. Best trial
costs = [trials.cost];
[~, best_idx] = min(costs);
best = trials(best_idx);
pos = best.pos;
C_i = best.C_i;
R_a_i = cnst.N_c ./ (cnst.mu_0 * abs(C_i));

fprintf('\n====== Results ======\n');
fprintf('Best: #%d (%s)\n', best_idx, best.name);

% Trial table
fprintf('\n--- Trial convergence ---\n');
for t = 1:n_trials
    m = ''; if t == best_idx, m = ' *'; end
    fprintf('%-16s  cost=%.6e  iters=%d  %.0fs%s\n', trials(t).name, trials(t).cost, trials(t).iters, trials(t).elapsed, m);
end

% Position spread
fprintf('\n--- Position spread ---\n');
all_pos = zeros(3, 6, n_trials);
for t = 1:n_trials, all_pos(:,:,t) = trials(t).pos; end
for i = 1:6
    pts = squeeze(all_pos(:,i,:));
    centroid = mean(pts, 2);
    fprintf('  P%d: max spread = %.2f um\n', i, max(vecnorm(pts - centroid))*1e6);
end

% Per-pole results
fprintf('\n--- Per-pole charge positions and C_i ---\n');
fprintf('%-4s %-6s  %8s  %12s  %12s\n', 'Pole','Type','|c|[um]','C_i','R_a');
for i = 1:6
    if cnst.pole_is_lower(i), ts='Lower'; else, ts='Upper'; end
    fprintf('P%d   %-6s  %8.1f  %12.4e  %12.3e\n', i, ts, norm(pos(:,i))*1e6, C_i(i), R_a_i(i));
end
lower = logical(cnst.pole_is_lower);
fprintf('\nLower avg |c|: %.1f um, R_a: %.3e\n', mean(vecnorm(pos(:,lower)))*1e6, mean(R_a_i(lower)));
fprintf('Upper avg |c|: %.1f um, R_a: %.3e\n', mean(vecnorm(pos(:,~lower)))*1e6, mean(R_a_i(~lower)));

% Fitting error
B_mat = build_basis(pos, w, 1e-7, px, py, pz, N);
b_model = B_mat * C_i;
err_vec = b_model - b_fem;
err_mag = sqrt(err_vec(1:N).^2 + err_vec(N+1:2*N).^2 + err_vec(2*N+1:3*N).^2);
fem_mag = sqrt(b_fem(1:N).^2 + b_fem(N+1:2*N).^2 + b_fem(2*N+1:3*N).^2);
rel = err_mag ./ fem_mag;
fprintf('\n--- Fitting error (%d nodes) ---\n', N);
fprintf('  Mean:   %.4f%%\n', 100*mean(rel));
fprintf('  Median: %.4f%%\n', 100*median(rel));
fprintf('  95th:   %.4f%%\n', 100*prctile(rel, 95));
fprintf('  Max:    %.4f%%\n', 100*max(rel));

%% 6. Save
sv.method = 'B-6x-perC';
sv.pos = pos; sv.C_i = C_i; sv.R_a_i = R_a_i;
sv.w = w; sv.I_apdl = I_apdl; sv.d_hat = d_hat;
sv.cost = best.cost; sv.cost_sphere = cost0;
sv.cube_half = cube_half; sv.N_nodes = N;
sv.err_mean = mean(rel); sv.err_max = max(rel);
sv.trials = trials; sv.best_idx = best_idx;
save(fullfile('..','..','data','B6x_hung_6C.mat'), '-struct', 'sv');
fprintf('\nSaved to data/B6x_hung_6C.mat\n');

%% Local functions

function [cost, C_i] = cost_perC(x, w, px, py, pz, b_fem, N)
    pos = reshape(x, 3, 6);
    B_mat = build_basis(pos, w, 1e-7, px, py, pz, N);
    C_i = B_mat \ b_fem;  % 6x1 linear LS
    cost = sum((B_mat * C_i - b_fem).^2);
end

function B = build_basis(pos, w, k_m, px, py, pz, N)
    % B: 3N x 6 matrix, column i = field from charge i with unit C
    B = zeros(3*N, 6);
    for i = 1:6
        dx = px - pos(1,i);
        dy = py - pos(2,i);
        dz = pz - pos(3,i);
        r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
        B(1:N, i)       = k_m * (-w(i)) * dx ./ r3;
        B(N+1:2*N, i)   = k_m * (-w(i)) * dy ./ r3;
        B(2*N+1:3*N, i) = k_m * (-w(i)) * dz ./ r3;
    end
end
