%% fit_J_fittedKI.m — Method [J]: Joint 6-coil fitting (Hung)
% VARIANT VERSION: see variants/README or analysis/README for context
%
%  Fit 18 shared charge positions + 6 per-coil C_k (VarPro)
%  simultaneously across all 6 coils.
%
%  Hung-specific: all 6 tips are SOURCES (N1/N2 swap in SOURC36).
%  I_vec = +1 for each coil, C_k will be negative (VarPro handles sign).
%  apdl_to_paper = identity mapping [1,2,3,4,5,6].
%
%  Run from magnetic_sim/hung/analysis/ directory.

%% 1. Load constants and fitted K_I
cnst = mt_constants();

% Use FITTED K_I (from fit_KI_v2.m), not ideal
KI_data = load('../../data/KI_fit.mat');
K_I = KI_data.K_I_fit;
fprintf('=== Method [J]: Joint 6-coil fitting (Hung) ===\n');
fprintf('Using FITTED K_I (from KI_fit.mat)\n');

cube_half = 100e-6;  % ±100 um fitting cube
fprintf('Fitting cube: +/-%.0f um, subsample 5000 nodes\n\n', cube_half*1e6);

coil_data = struct();
for k = 1:6
    coilname = sprintf('coil%d', k);
    data = import_ansys_data(fullfile('..','..','results',coilname,'filleted'), 'wp', coilname);

    % Select nodes in ±cube_half cube, subsample if too many
    mask = abs(data.x) < cube_half & abs(data.y) < cube_half & abs(data.z) < cube_half;
    N_max = 5000;  % subsample for speed
    if sum(mask) > N_max
        idx_all = find(mask);
        rng(0);
        idx_sub = idx_all(randperm(length(idx_all), N_max));
        mask = false(size(mask));
        mask(idx_sub) = true;
    end

    coil_data(k).px = data.x(mask);
    coil_data(k).py = data.y(mask);
    coil_data(k).pz = data.z(mask);
    coil_data(k).bx = data.bx(mask);
    coil_data(k).by = data.by(mask);
    coil_data(k).bz = data.bz(mask);
    coil_data(k).bmag = sqrt(data.bx(mask).^2 + data.by(mask).^2 + data.bz(mask).^2);
    coil_data(k).b_fem = [data.bx(mask); data.by(mask); data.bz(mask)];
    coil_data(k).N = sum(mask);

    % Weight vector: K_I * I_vec_k (identity mapping for Hung)
    I_vec_k = zeros(6, 1);
    I_vec_k(k) = 1;
    coil_data(k).KI_w = K_I * I_vec_k;

    fprintf('  Coil%d (P%d): %d FEM nodes, |B| mean=%.3f mT\n', ...
        k, k, coil_data(k).N, mean(coil_data(k).bmag)*1e3);
end
fprintf('\n');

%% 2. Initial positions from per-coil ell (KI_fit.mat)
alpha = cnst.alpha;
d_hat = zeros(3, 6);
for i = 1:6
    theta = cnst.pole_angles(i) * pi/180;
    z_sign = sign(cnst.pole_tip_z_wp(i));
    d_hat(:,i) = [cos(theta)*sin(alpha); sin(theta)*sin(alpha); z_sign*cos(alpha)];
end

% Use per-layer ell from KI_fit
pos_init = KI_data.pos;  % per-layer ell positions
fprintf('--- Initial positions (per-layer ell) ---\n');
fprintf('  Lower |c|: %.1f um, Upper |c|: %.1f um\n\n', ...
    KI_data.ell_lower*1e6, KI_data.ell_upper*1e6);

%% 3. Joint 6-coil optimization (18 positions, per-coil C_k via VarPro)
fprintf('===== Joint 6-coil optimization =====\n');

cost_fn = @(x) cost_joint(x, coil_data, cnst.k_m);

fms_opts = optimset('TolX', 1e-10, 'TolFun', 1e-22, ...
    'MaxIter', 50000, 'MaxFunEvals', 50000, 'Display', 'off');

% Initial conditions (3 trials to keep runtime reasonable)
x0_sphere = pos_init(:);
rng(42);  x0_n1 = pos_init(:) + 50e-6*randn(18,1);
rng(123); x0_n2 = pos_init(:) + 100e-6*randn(18,1);

inits = {x0_sphere, x0_n1, x0_n2};
init_names = {'Sphere', '+50um noise', '+100um noise'};

n_trials = length(inits);
trials = struct([]);

for t = 1:n_trials
    fprintf('\n----- Trial %d/%d: %s -----\n', t, n_trials, init_names{t});
    tic;
    [x_opt, fval, ~, output] = fminsearch(cost_fn, inits{t}, fms_opts);
    elapsed = toc;

    [~, C_k] = cost_joint(x_opt, coil_data, cnst.k_m);

    trials(t).name    = init_names{t};
    trials(t).pos     = reshape(x_opt, 3, 6);
    trials(t).cost    = fval;
    trials(t).C_k     = C_k;
    trials(t).iters   = output.iterations;
    trials(t).fevals  = output.funcCount;
    trials(t).elapsed = elapsed;

    fprintf('  Cost = %.6e  (iters %d, fevals %d, %.0fs)\n', ...
        fval, output.iterations, output.funcCount, elapsed);
end

%% 4. Select best trial
costs = [trials.cost];
[~, best_idx] = min(costs);
best = trials(best_idx);
pos_best = best.pos;
C_k_best = best.C_k;
R_a_k = cnst.N_c ./ (cnst.mu_0 * abs(C_k_best));

fprintf('\n====== Results Summary ======\n');
fprintf('Best trial: #%d (%s)\n', best_idx, best.name);

%% 5. Trial convergence table
fprintf('\n--- Trial convergence ---\n');
fprintf('%-16s  %12s  %8s  %8s  %6s\n', 'Init', 'Cost', 'Iters', 'FEvals', 'Time');
for t = 1:n_trials
    marker = '';
    if t == best_idx, marker = ' *'; end
    fprintf('%-16s  %12.6e  %8d  %8d  %5.0fs%s\n', ...
        trials(t).name, trials(t).cost, trials(t).iters, trials(t).fevals, ...
        trials(t).elapsed, marker);
end

%% 6. Position consistency (identifiability)
fprintf('\n--- Position consistency (max spread per pole) ---\n');
all_pos = zeros(3, 6, n_trials);
for t = 1:n_trials
    all_pos(:,:,t) = trials(t).pos;
end
for i = 1:6
    pts = squeeze(all_pos(:,i,:));
    centroid = mean(pts, 2);
    dists = vecnorm(pts - centroid);
    fprintf('  P%d: max spread = %.2f um\n', i, max(dists)*1e6);
end

%% 7. Per-pole charge positions
fprintf('\n--- Per-pole charge positions ---\n');
fprintf('%-4s %-6s  %8s %8s %8s  %8s\n', 'Pole', 'Type', 'x[um]', 'y[um]', 'z[um]', '|c|[um]');
for i = 1:6
    ci = pos_best(:,i);
    if cnst.pole_is_lower(i), type_str = 'Lower'; else, type_str = 'Upper'; end
    fprintf('P%d   %-6s  %8.1f %8.1f %8.1f  %8.1f\n', ...
        i, type_str, ci(1)*1e6, ci(2)*1e6, ci(3)*1e6, norm(ci)*1e6);
end
lower_mask = logical(cnst.pole_is_lower);
fprintf('\n  Lower avg |c|: %.1f um\n', mean(vecnorm(pos_best(:, lower_mask)))*1e6);
fprintf('  Upper avg |c|: %.1f um\n', mean(vecnorm(pos_best(:, ~lower_mask)))*1e6);

%% 8. Per-coil C_k and R_a
fprintf('\n--- Per-coil C_k and R_a ---\n');
fprintf('%-8s %-6s  %12s  %12s  %s\n', 'Coil', 'Type', 'C_k', 'R_a [A/Wb]', 'Polarity');
for k = 1:6
    if cnst.pole_is_lower(k), type_str = 'Lower'; else, type_str = 'Upper'; end
    if C_k_best(k) < 0, pol = 'SOURCE'; else, pol = 'SINK'; end
    fprintf('Coil%d   %-6s  %12.4e  %12.3e  %s\n', ...
        k, type_str, C_k_best(k), R_a_k(k), pol);
end
fprintf('\n  Lower avg R_a: %.3e\n', mean(R_a_k(lower_mask)));
fprintf('  Upper avg R_a: %.3e\n', mean(R_a_k(~lower_mask)));

%% 9. Per-coil fitting error
fprintf('\n--- Per-coil fitting error ---\n');
fprintf('%-8s  %10s\n', 'Coil', 'Mean Err[%%]');
err_per_coil = zeros(1, 6);
for k = 1:6
    N_k = coil_data(k).N;
    b_unit_k = eval_charge_field(pos_best, coil_data(k).KI_w, cnst.k_m, ...
        [coil_data(k).px, coil_data(k).py, coil_data(k).pz]);
    b_model_k = C_k_best(k) * b_unit_k;
    err_vec = b_model_k - coil_data(k).b_fem;
    err_mag = sqrt(err_vec(1:N_k).^2 + err_vec(N_k+1:2*N_k).^2 + err_vec(2*N_k+1:3*N_k).^2);
    rel = err_mag ./ coil_data(k).bmag;
    err_per_coil(k) = mean(rel);
    fprintf('Coil%d     %10.2f\n', k, err_per_coil(k)*100);
end
fprintf('Mean:     %10.2f\n', mean(err_per_coil)*100);

%% 10. Save results
save_data.method      = 'J';
save_data.description = 'Joint 6-coil fitting, Hung design (all tips SOURCE)';
save_data.pos         = pos_best;
save_data.C_k         = C_k_best;
save_data.R_a_k       = R_a_k;
save_data.ell_lower   = KI_data.ell_lower;
save_data.ell_upper   = KI_data.ell_upper;
save_data.d_hat       = d_hat;
save_data.K_I         = K_I;
save_data.cube_half   = cube_half;
save_data.N_per_coil  = [coil_data.N];
save_data.err_per_coil = err_per_coil;
save_data.err_mean    = mean(err_per_coil);
save_data.trials      = trials;
save_data.best_idx    = best_idx;

save(fullfile('..','..','data','variants','joint_6coil_fit.mat'), '-struct', 'save_data');
fprintf('\nSaved to data/variants/joint_6coil_fit.mat\n');

%% ===== Local functions =====

function cost = ell_cost(ell, d_hat, w, k_m, p_wp, b_fem)
    pos = ell * d_hat;
    b_unit = eval_charge_field(pos, w, k_m, p_wp);
    C = (b_unit' * b_fem) / (b_unit' * b_unit);
    cost = sum((C * b_unit - b_fem).^2);
end

function [cost, C_k] = cost_joint(x, coil_data, k_m)
    pos = reshape(x, 3, 6);
    cost = 0;
    C_k = zeros(6, 1);
    for k = 1:6
        p_wp = [coil_data(k).px, coil_data(k).py, coil_data(k).pz];
        b_unit = eval_charge_field(pos, coil_data(k).KI_w, k_m, p_wp);
        C_k(k) = (b_unit' * coil_data(k).b_fem) / (b_unit' * b_unit);
        r = C_k(k) * b_unit - coil_data(k).b_fem;
        cost = cost + sum(r.^2);
    end
end

function b = eval_charge_field(pos, w, k_m, p_wp)
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

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
