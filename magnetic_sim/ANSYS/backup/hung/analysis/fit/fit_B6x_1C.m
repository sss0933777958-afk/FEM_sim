%% fit_B6x_1C.m — Method [B-6x] for Hung: alternating superposition
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
%
%  Hung: all tips are SOURCE. Use alternating I_apdl = [+1,-1,+1,-1,-1,+1]
%  so I_diss = [-1,+1,-1,+1,+1,-1], sum=0 → w = K_I × I_diss = I_diss
%  All poles |w|=1, K_I drops out.
%
%  Model: c_i = ell * d_hat_i + delta_i  (18 delta params)
%  C: 1 shared amplitude (VarPro)
%  ell: fixed from Method [A]

%% 1. Constants
cnst = mt_constants();
K_I = eye(6) - ones(6)/6;
cube_half = 50e-6;   % ±50 um

fprintf('=== Method [B-6x] Hung: alternating superposition ===\n\n');

% Alternating APDL current: [+1,-1,+1,-1,-1,+1]
I_apdl = [1; -1; 1; -1; -1; 1];
% All SOURCE: +1A → I_diss=-1, -1A → I_diss=+1
I_diss = -I_apdl;
w = K_I * I_diss;  % should equal I_diss since sum=0

fprintf('I_apdl = [%+d, %+d, %+d, %+d, %+d, %+d]\n', I_apdl);
fprintf('I_diss = [%+d, %+d, %+d, %+d, %+d, %+d]\n', I_diss);
fprintf('w      = [%+.4f, %+.4f, %+.4f, %+.4f, %+.4f, %+.4f]\n', w);
fprintf('sum(I_diss) = %d, max|w - I_diss| = %.2e\n\n', sum(I_diss), max(abs(w-I_diss)));

%% 2. Build alternating superposition from 6 single-coil FEM data
fprintf('Loading 6 coils and computing alternating superposition...\n');

% Load first coil to get coordinates
d1 = import_ansys_data(fullfile('..','..','results','coil1','filleted'), 'wp', 'coil1');
mask = abs(d1.x) < cube_half & abs(d1.y) < cube_half & abs(d1.z) < cube_half;
N = sum(mask);
p_wp = [d1.x(mask), d1.y(mask), d1.z(mask)];

bx_alt = I_apdl(1) * d1.bx(mask);
by_alt = I_apdl(1) * d1.by(mask);
bz_alt = I_apdl(1) * d1.bz(mask);
fprintf('  Coil1 (x%+d): %d nodes\n', I_apdl(1), N);

for k = 2:6
    dk = import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'), ...
        'wp', sprintf('coil%d',k));
    bx_alt = bx_alt + I_apdl(k) * dk.bx(mask);
    by_alt = by_alt + I_apdl(k) * dk.by(mask);
    bz_alt = bz_alt + I_apdl(k) * dk.bz(mask);
    fprintf('  Coil%d (x%+d)\n', k, I_apdl(k));
end
b_fem = [bx_alt; by_alt; bz_alt];

b_mag = sqrt(bx_alt.^2 + by_alt.^2 + bz_alt.^2);
fprintf('\nAlternating |B|: min=%.3f, mean=%.3f, max=%.3f mT\n', ...
    min(b_mag)*1e3, mean(b_mag)*1e3, max(b_mag)*1e3);
fprintf('Data: %d nodes, %d equations\n\n', N, 3*N);

%% 3. Build d_hat
alpha = cnst.alpha;
d_hat = zeros(3, 6);
for i = 1:6
    theta = cnst.pole_angles(i) * pi/180;
    z_sign = sign(cnst.pole_tip_z_wp(i));
    d_hat(:,i) = [cos(theta)*sin(alpha); sin(theta)*sin(alpha); z_sign*cos(alpha)];
end

%% 4. Fix ell from per-layer averages + load [J] positions for init
KI_data = load('../../data/KI_fit.mat');
ell_lower = KI_data.ell_lower;
ell_upper = KI_data.ell_upper;
J_data = load('../../data/J_ideal_fit.mat');
fprintf('ell_lower (fixed): %.1f um, ell_upper: %.1f um\n', ell_lower*1e6, ell_upper*1e6);

%% 5. Sphere baseline (delta=0, per-layer ell)
% Build per-layer ell sphere positions
pos_sphere = zeros(3,6);
for i = 1:6
    if cnst.pole_is_lower(i)
        pos_sphere(:,i) = ell_lower * d_hat(:,i);
    else
        pos_sphere(:,i) = ell_upper * d_hat(:,i);
    end
end

[cost0, C0] = cost_with_C_pos(pos_sphere, w, p_wp, b_fem);
R_a_0 = cnst.N_c / (cnst.mu_0 * abs(C0));
fprintf('Sphere baseline: cost=%.6e, C=%.4e, R_a=%.3e\n\n', cost0, C0, R_a_0);

%% 6. Optimize 18 position parameters (direct positions, not delta)
cost_fn = @(x) cost_pos(x, w, p_wp, b_fem);

opts = optimset('TolX', 1e-10, 'TolFun', 1e-22, ...
    'MaxIter', 50000, 'MaxFunEvals', 50000, 'Display', 'off');

% Initial conditions: sphere + [J] positions + noise
x0_list = {};
x0_names = {};
x0_list{end+1} = pos_sphere(:);          x0_names{end+1} = 'sphere';
x0_list{end+1} = J_data.pos(:);          x0_names{end+1} = '[J] positions';
rng(42); x0_list{end+1} = J_data.pos(:) + 20e-6*randn(18,1); x0_names{end+1} = '[J]+20um noise';

n_trials = length(x0_list);
trials = struct([]);
fprintf('====== Starting %d trials ======\n', n_trials);

for t = 1:n_trials
    fprintf('\n----- Trial %d: %s -----\n', t, x0_names{t});
    tic;
    [x_opt, fval, ~, output] = fminsearch(cost_fn, x0_list{t}, opts);
    elapsed = toc;
    pos = reshape(x_opt, 3, 6);
    trials(t).x_opt = x_opt;
    trials(t).cost = fval;
    trials(t).pos = pos;
    trials(t).iters = output.iterations;
    trials(t).fevals = output.funcCount;
    trials(t).elapsed = elapsed;
    trials(t).name = x0_names{t};
    fprintf('  Cost=%.6e (iters %d, fevals %d, %.0fs)\n', fval, output.iterations, output.funcCount, elapsed);
end

%% 7. Select best
costs = [trials.cost];
[~, best_idx] = min(costs);
best = trials(best_idx);
[~, C_best] = cost_with_C_pos(best.pos, w, p_wp, b_fem);
R_a_best = cnst.N_c / (cnst.mu_0 * abs(C_best));

fprintf('\n====== Results ======\n');
fprintf('Best: #%d (%s)\n', best_idx, best.name);

% Trial table
fprintf('\n--- Trial convergence ---\n');
fprintf('%-16s  %12s  %8s  %6s\n', 'Init', 'Cost', 'Iters', 'Time');
for t = 1:n_trials
    m = ''; if t == best_idx, m = ' *'; end
    fprintf('%-16s  %12.6e  %8d  %5.0fs%s\n', trials(t).name, trials(t).cost, trials(t).iters, trials(t).elapsed, m);
end

% Position consistency
fprintf('\n--- Position spread ---\n');
all_pos = zeros(3, 6, n_trials);
for t = 1:n_trials, all_pos(:,:,t) = trials(t).pos; end
for i = 1:6
    pts = squeeze(all_pos(:,i,:));
    centroid = mean(pts, 2);
    fprintf('  P%d: max spread = %.2f um\n', i, max(vecnorm(pts - centroid))*1e6);
end

%% 8. Detailed results
pos = best.pos;

fprintf('\n--- Optimal parameters ---\n');
fprintf('  C = %.4e, R_a = %.3e A/Wb\n', C_best, R_a_best);

fprintf('\n--- Per-pole charge positions ---\n');
fprintf('%-4s %-6s  %8s %8s %8s  %8s\n', 'Pole','Type','x[um]','y[um]','z[um]','|c|[um]');
for i = 1:6
    ci = pos(:,i);
    if cnst.pole_is_lower(i), ts='Lower'; else, ts='Upper'; end
    fprintf('P%d   %-6s  %8.1f %8.1f %8.1f  %8.1f\n', ...
        i, ts, ci(1)*1e6, ci(2)*1e6, ci(3)*1e6, norm(ci)*1e6);
end
lower_mask = logical(cnst.pole_is_lower);
fprintf('\n  Lower avg |c|: %.1f um\n', mean(vecnorm(pos(:,lower_mask)))*1e6);
fprintf('  Upper avg |c|: %.1f um\n', mean(vecnorm(pos(:,~lower_mask)))*1e6);

% Fitting error
b_unit = eval_field(pos, w, cnst.k_m, p_wp);
b_model = C_best * b_unit;
err_vec = b_model - b_fem;
err_mag = sqrt(err_vec(1:N).^2 + err_vec(N+1:2*N).^2 + err_vec(2*N+1:3*N).^2);
fem_mag = sqrt(b_fem(1:N).^2 + b_fem(N+1:2*N).^2 + b_fem(2*N+1:3*N).^2);
valid = fem_mag > 1e-6 * max(fem_mag);
rel = err_mag ./ fem_mag;

fprintf('\n--- Fitting error (%d nodes) ---\n', N);
fprintf('  Mean:   %.2f%%\n', 100*mean(rel(valid)));
fprintf('  Median: %.2f%%\n', 100*median(rel(valid)));
fprintf('  95th:   %.2f%%\n', 100*prctile(rel(valid), 95));
fprintf('  Max:    %.2f%%\n', 100*max(rel(valid)));

fprintf('\n--- Comparison ---\n');
fprintf('  Sphere: cost=%.6e, R_a=%.3e\n', cost0, R_a_0);
fprintf('  [B-6x]: cost=%.6e, R_a=%.3e\n', best.cost, R_a_best);
fprintf('  Reduction: %.2f%%\n', 100*(1 - best.cost/cost0));

%% 9. Save
sv.method = 'B-6x-hung';
sv.description = 'Alternating superposition [+1,-1,+1,-1,-1,+1] for all-SOURCE Hung';
sv.ell_lower = ell_lower;
sv.ell_upper = ell_upper;
sv.d_hat = d_hat;
sv.pos = best.pos;
sv.C = C_best;
sv.R_a = R_a_best;
sv.cost = best.cost;
sv.cost_sphere = cost0;
sv.w = w;
sv.I_apdl = I_apdl;
sv.I_diss = I_diss;
sv.N_nodes = N;
sv.cube_half = cube_half;
sv.err_mean = mean(rel(valid));
sv.err_max = max(rel(valid));
sv.trials = trials;
save(fullfile('..','..','data','B6x_hung_1C.mat'), '-struct', 'sv');
fprintf('\nSaved to data/B6x_hung_1C.mat\n');

%% ===== Local functions =====

function cost = cost_pos(x, w, p_wp, b_fem)
    pos = reshape(x, 3, 6);
    b_unit = eval_field(pos, w, 1e-7, p_wp);
    C = (b_unit' * b_fem) / (b_unit' * b_unit);
    cost = sum((C * b_unit - b_fem).^2);
end

function [cost, C] = cost_with_C_pos(pos, w, p_wp, b_fem)
    b_unit = eval_field(pos, w, 1e-7, p_wp);
    C = (b_unit' * b_fem) / (b_unit' * b_unit);
    cost = sum((C * b_unit - b_fem).^2);
end

function b = eval_field(pos, w, k_m, p_wp)
    N = size(p_wp, 1);
    bx = zeros(N,1); by = bx; bz = bx;
    for i = 1:6
        dx = p_wp(:,1) - pos(1,i);
        dy = p_wp(:,2) - pos(2,i);
        dz = p_wp(:,3) - pos(3,i);
        r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
        bx = bx + (-w(i)) * dx ./ r3;
        by = by + (-w(i)) * dy ./ r3;
        bz = bz + (-w(i)) * dz ./ r3;
    end
    b = k_m * [bx; by; bz];
end
