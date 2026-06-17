%% fit_all6_with_bias.m — Method [B-6x]: all-excitation b!=0 fitting
%
%  Model: c_i = ell * d_hat_i + delta_i   (i = 1..6)
%    - ell:     fixed from Method [A] (~835 um)
%    - delta_i: 3D offset per pole, 6x3 = 18 free params
%    - C:       1 shared amplitude (solved analytically via VarPro)
%
%  Data: all6 simultaneous excitation (CURR_ARRAY=[1,1,1,1,1,1])
%    Linear material => B(all6) = B(coil1) + ... + B(coil6)
%    100 um cube at WP center
%
%  Weight vector derivation:
%    APDL all +1A => model I_diss = [+1,-1,+1,-1,-1,+1] (coil_sign applied)
%    w = K_I * I_diss = I_diss  (since sum(I_diss) = 0)
%    All 6 poles have |w| = 1  (vs single-coil: active 5/6, passive 1/6)
%
%  Key advantage over [B-sc]: all poles equally constrained (|w|=1)
%  Key advantage over [D]:    single shared C (closer to dissertation model)

%% 1. Load FEM data (all6 or superposition)
cnst = mt_constants();
K_I = eye(6) - ones(6)/6;
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

fprintf('=== Method [B-6x]: all-excitation b!=0 fitting ===\n');

% Build I_diss (model current vector with coil_sign correction)
I_diss = zeros(6,1);
for k = 1:6
    paper_idx = apdl_to_paper_idx(k);
    coil_sign = 2*cnst.pole_is_lower(paper_idx) - 1;
    I_diss(paper_idx) = coil_sign;
end
w = K_I * I_diss;  % = I_diss since sum(I_diss) = 0
fprintf('I_diss = [%+d, %+d, %+d, %+d, %+d, %+d]\n', I_diss);
fprintf('w      = [%+.4f, %+.4f, %+.4f, %+.4f, %+.4f, %+.4f]\n', w);
fprintf('Verify: max|w - I_diss| = %.2e (should be ~0)\n', max(abs(w - I_diss)));

% Try loading all6 directly; fall back to superposition of singles
cube_half = 50e-6;
use_superposition = false;

try
    d = import_ansys_data(fullfile('..', 'results', 'all6'), 'wp', 'all6');
    [air_mask, ~] = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize', false));
    z_wp = d.z - cnst.SPH_OFST;
    mask = air_mask & abs(d.x) < cube_half & abs(d.y) < cube_half & abs(z_wp) < cube_half;

    p_wp = [d.x(mask), d.y(mask), z_wp(mask)];
    b_fem = [d.bx(mask); d.by(mask); d.bz(mask)];
    N = sum(mask);
    fprintf('Loaded all6 directly: %d nodes in 100 um cube\n', N);
catch
    fprintf('all6 data not found, computing from superposition of 6 singles...\n');
    use_superposition = true;

    % Load first coil to get coordinates
    d1 = import_ansys_data(fullfile('..', 'results', 'coil1'), 'wp', 'coil1');
    [air_mask, ~] = filter_iron_nodes(d1.x, d1.y, d1.z, cnst, struct('visualize', false));
    z_wp = d1.z - cnst.SPH_OFST;
    mask = air_mask & abs(d1.x) < cube_half & abs(d1.y) < cube_half & abs(z_wp) < cube_half;

    p_wp = [d1.x(mask), d1.y(mask), z_wp(mask)];
    N = sum(mask);

    % Sum B-fields from all 6 coils
    bx_sum = zeros(N, 1);
    by_sum = zeros(N, 1);
    bz_sum = zeros(N, 1);
    for k = 1:6
        coil_name = sprintf('coil%d', k);
        dk = import_ansys_data(fullfile('..', 'results', coil_name), 'wp', coil_name);
        bx_sum = bx_sum + dk.bx(mask);
        by_sum = by_sum + dk.by(mask);
        bz_sum = bz_sum + dk.bz(mask);
    end
    b_fem = [bx_sum; by_sum; bz_sum];
    fprintf('Superposition computed: %d nodes in 100 um cube\n', N);
end

% Quick field statistics
b_mag = sqrt(b_fem(1:N).^2 + b_fem(N+1:2*N).^2 + b_fem(2*N+1:3*N).^2);
fprintf('|B_fem|: min=%.3f mT, mean=%.3f mT, max=%.3f mT\n', ...
    min(b_mag)*1e3, mean(b_mag)*1e3, max(b_mag)*1e3);

%% 2. Build d_hat (6 direction unit vectors)
alpha = cnst.alpha;
d_hat = zeros(3, 6);
for i = 1:6
    theta = cnst.pole_angles(i) * pi/180;
    z_sign = sign(cnst.pole_tip_z_wp(i));
    d_hat(:,i) = [cos(theta)*sin(alpha); sin(theta)*sin(alpha); z_sign*cos(alpha)];
end

%% 3. Fix ell from Method [A]
mat_A = fullfile('..', 'data', 'charge_model_fit.mat');
if isfile(mat_A)
    prev = load(mat_A);
    ell_fixed = prev.fitA.ell;
    fprintf('\nMethod [A] ell = %.1f um (from charge_model_fit.mat)\n', ell_fixed*1e6);
else
    error('charge_model_fit.mat not found — run fit_charge_model.m first');
end

%% 4. Verify: delta=0 (sphere model) baseline
cost_fn = @(x) cost_bias(x, ell_fixed, d_hat, w, p_wp, b_fem, cnst);
cost0 = cost_fn(zeros(18,1));
[~, C0] = cost_bias_with_C(zeros(18,1), ell_fixed, d_hat, w, p_wp, b_fem, cnst);
R_a_0 = cnst.N_c / (cnst.mu_0 * abs(C0));
fprintf('Sphere baseline: cost = %.6e  (C=%.4e, R_a=%.3e)\n', cost0, C0, R_a_0);

%% 5. Optimize 18 delta parameters
opts = optimset('TolX', 1e-10, 'TolFun', 1e-22, ...
    'MaxIter', 50000, 'MaxFunEvals', 500000, 'Display', 'iter');

% --- Build initial conditions ---
x0_list = {};
x0_names = {};

% (a) zeros = sphere model
x0_list{end+1} = zeros(18,1);
x0_names{end+1} = 'zeros (sphere)';

% (b) +50 um random noise
rng(42);
x0_list{end+1} = 50e-6 * randn(18,1);
x0_names{end+1} = '+50um noise';

% (c) +100 um random noise
x0_list{end+1} = 100e-6 * randn(18,1);
x0_names{end+1} = '+100um noise';

% --- Run all trials ---
n_trials = length(x0_list);
trials = struct([]);
fprintf('\n====== Starting %d optimization trials ======\n', n_trials);

for t = 1:n_trials
    fprintf('\n----- Trial %d/%d: %s -----\n', t, n_trials, x0_names{t});
    tic;
    [x_opt, fval, exitflag, output] = fminsearch(cost_fn, x0_list{t}, opts);
    elapsed = toc;

    delta = reshape(x_opt, 3, 6);
    pos = ell_fixed * d_hat + delta;

    trials(t).x_opt    = x_opt;
    trials(t).cost     = fval;
    trials(t).exitflag = exitflag;
    trials(t).iters    = output.iterations;
    trials(t).fevals   = output.funcCount;
    trials(t).name     = x0_names{t};
    trials(t).delta    = delta;
    trials(t).pos      = pos;
    trials(t).elapsed  = elapsed;

    fprintf('  Cost = %.6e  (iters %d, fevals %d, %.0fs)\n', ...
        fval, output.iterations, output.funcCount, elapsed);
end

%% 6. Select best trial
costs = [trials.cost];
[best_cost, best_idx] = min(costs);
best = trials(best_idx);

fprintf('\n====== Results Summary ======\n');
fprintf('Best trial: #%d (%s)\n', best_idx, best.name);
fprintf('Cost: %.6e -> %.6e  (%.2f%% reduction from sphere)\n', ...
    cost0, best_cost, 100*(1 - best_cost/cost0));

% Trial convergence table
fprintf('\n--- Trial convergence ---\n');
fprintf('%-20s  %12s  %8s  %8s  %6s\n', 'Init', 'Cost', 'Iters', 'FEvals', 'Time');
for t = 1:n_trials
    marker = '';
    if t == best_idx, marker = ' *'; end
    fprintf('%-20s  %12.6e  %8d  %8d  %5.0fs%s\n', ...
        trials(t).name, trials(t).cost, trials(t).iters, trials(t).fevals, ...
        trials(t).elapsed, marker);
end

% Position consistency across trials
if n_trials > 1
    fprintf('\n--- Position consistency (max spread per pole) ---\n');
    all_pos = zeros(3, 6, n_trials);
    for t = 1:n_trials
        all_pos(:,:,t) = trials(t).pos;
    end
    for i = 1:6
        pts = squeeze(all_pos(:,i,:));
        centroid = mean(pts, 2);
        dists = vecnorm(pts - centroid);
        fprintf('  %s: max spread = %.2f um\n', cnst.pole_labels{i}, max(dists)*1e6);
    end
end

%% 7. Detailed analysis of best result
delta = best.delta;
pos   = best.pos;
[~, C_best] = cost_bias_with_C(best.x_opt, ell_fixed, d_hat, w, p_wp, b_fem, cnst);
R_a_best = cnst.N_c / (cnst.mu_0 * abs(C_best));

fprintf('\n--- Optimal parameters ---\n');
fprintf('  ell (fixed): %.1f um\n', ell_fixed*1e6);
fprintf('  C = %.4e,  R_a = %.3e A/Wb\n', C_best, R_a_best);

fprintf('\n--- Per-pole charge positions ---\n');
fprintf('%-4s  %8s %8s %8s  %8s  %8s  %8s %8s\n', ...
    'Pole', 'x[um]', 'y[um]', 'z[um]', '|c|[um]', '|d|[um]', 'd_par', 'd_perp');
for i = 1:6
    ci = pos(:,i);
    di = delta(:,i);
    ell_i = norm(ci);
    d_par = dot(di, d_hat(:,i));
    d_perp = norm(di - d_par * d_hat(:,i));
    fprintf('%-4s  %8.1f %8.1f %8.1f  %8.1f  %8.1f  %8.1f %8.1f\n', ...
        cnst.pole_labels{i}, ci(1)*1e6, ci(2)*1e6, ci(3)*1e6, ...
        ell_i*1e6, norm(di)*1e6, d_par*1e6, d_perp*1e6);
end

lower_mask = logical(cnst.pole_is_lower);
fprintf('\n  Lower avg |c|: %.1f um\n', mean(vecnorm(pos(:, lower_mask)))*1e6);
fprintf('  Upper avg |c|: %.1f um\n', mean(vecnorm(pos(:, ~lower_mask)))*1e6);

% Compute fitting error
b_unit = eval_charge_field(pos, w, cnst.k_m, p_wp);
b_model = C_best * b_unit;
err_vec = b_model - b_fem;
err_mag = sqrt(err_vec(1:N).^2 + err_vec(N+1:2*N).^2 + err_vec(2*N+1:3*N).^2);
fem_mag = sqrt(b_fem(1:N).^2 + b_fem(N+1:2*N).^2 + b_fem(2*N+1:3*N).^2);
valid = fem_mag > 1e-6 * max(fem_mag);
rel = err_mag ./ fem_mag;
rel(~valid) = 0;

fprintf('\n--- Fitting error (100 um cube, %d nodes) ---\n', N);
fprintf('  Mean:   %.2f%%\n', 100*mean(rel(valid)));
fprintf('  Median: %.2f%%\n', 100*median(rel(valid)));
fprintf('  95th:   %.2f%%\n', 100*prctile(rel(valid), 95));
fprintf('  Max:    %.2f%%\n', 100*max(rel(valid)));

%% 8. Comparison with Method [A] (sphere, b=0)
fprintf('\n--- Comparison with Method [A] (sphere, b=0) ---\n');
fprintf('  [A]    cost: %.6e,  R_a = %.3e\n', cost0, R_a_0);
fprintf('  [B-6x] cost: %.6e,  R_a = %.3e\n', best_cost, R_a_best);
fprintf('  Cost reduction: %.2f%%\n', 100*(1 - best_cost/cost0));

%% 9. Save results
save_data.method      = 'B-6x';
save_data.description = 'All-excitation b!=0 fitting (CURR_ARRAY=[1,1,1,1,1,1])';
save_data.ell_fixed   = ell_fixed;
save_data.d_hat       = d_hat;
save_data.delta       = best.delta;
save_data.pos         = best.pos;
save_data.C           = C_best;
save_data.R_a         = R_a_best;
save_data.cost        = best_cost;
save_data.cost_sphere = cost0;
save_data.w           = w;
save_data.I_diss      = I_diss;
save_data.N_nodes     = N;
save_data.err_mean    = mean(rel(valid));
save_data.err_median  = median(rel(valid));
save_data.err_max     = max(rel(valid));
save_data.trials      = trials;
save_data.K_I         = K_I;
save_data.use_superposition = use_superposition;

save(fullfile('..', 'data', 'all6_bias_fit.mat'), '-struct', 'save_data');
fprintf('\nSaved to data/all6_bias_fit.mat\n');

%% ===== Local functions =====

function cost = cost_bias(x, ell, d_hat, w, p_wp, b_fem, cnst)
% Cost function: 18 delta params -> scalar cost (VarPro for C)
    delta = reshape(x, 3, 6);
    pos = ell * d_hat + delta;

    b_unit = eval_charge_field(pos, w, cnst.k_m, p_wp);
    C = (b_unit' * b_fem) / (b_unit' * b_unit);

    r = C * b_unit - b_fem;
    cost = sum(r.^2);
end

function [cost, C] = cost_bias_with_C(x, ell, d_hat, w, p_wp, b_fem, cnst)
% Same as cost_bias but also returns optimal C
    delta = reshape(x, 3, 6);
    pos = ell * d_hat + delta;

    b_unit = eval_charge_field(pos, w, cnst.k_m, p_wp);
    C = (b_unit' * b_fem) / (b_unit' * b_unit);

    r = C * b_unit - b_fem;
    cost = sum(r.^2);
end

function b = eval_charge_field(pos, w, k_m, p_wp)
% B-field from 6 point charges at given positions (unit amplitude C=1)
%   pos:  3x6 charge positions [m]
%   w:    6x1 weight vector (from K_I * I_vec)
%   k_m:  mu_0/(4*pi) = 1e-7
%   p_wp: Nx3 field points [m]
%   Returns: b = [bx; by; bz] (3N x 1)
    N = size(p_wp, 1);
    bx = zeros(N, 1);
    by = zeros(N, 1);
    bz = zeros(N, 1);

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
