%% test_joint_6coil_fit.m — Joint 6-coil 19-parameter fitting [J]
%
%  Fit 18 shared charge positions + 6 per-coil scaling factors C_k
%  simultaneously across all 6 coils (19 free params: 18 optimized by
%  fminsearch, 6 C_k solved analytically via Variable Projection).
%
%  The 19th "free" param count comes from the fact that each C_k is a
%  function of R_a_k, so 18 positions + 6 analytically-determined C_k
%  effectively means 18 optimized + 1 R_a per coil (6 total), but the
%  optimizer only sees 18 variables.
%
%  Run from magnetic_sim/hexapole-long2016/analysis/ directory.

clear; clc; close all;

%% 1. Load constants and FEM data from all 6 coils
cnst = mt_constants();
K_I = eye(6) - ones(6)/6;
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];
cube_half = 50e-6;  % 100 um cube fitting region

fprintf('Loading FEM data for all 6 coils...\n');

coil_data = struct();
for k = 1:6
    coilname = sprintf('coil%d', k);
    data = import_ansys_data(fullfile('..', 'results', coilname), 'wp', coilname);
    [air_mask, ~] = filter_iron_nodes(data.x, data.y, data.z, cnst, ...
        struct('visualize', false));
    z_wp = data.z - cnst.SPH_OFST;

    mask = air_mask & ...
        abs(data.x) < cube_half & abs(data.y) < cube_half & abs(z_wp) < cube_half;

    coil_data(k).px = data.x(mask);
    coil_data(k).py = data.y(mask);
    coil_data(k).pz = z_wp(mask);
    coil_data(k).bx = data.bx(mask);
    coil_data(k).by = data.by(mask);
    coil_data(k).bz = data.bz(mask);
    coil_data(k).bmag = sqrt(data.bx(mask).^2 + data.by(mask).^2 + data.bz(mask).^2);
    coil_data(k).b_fem = [data.bx(mask); data.by(mask); data.bz(mask)];
    coil_data(k).N = sum(mask);

    % Build KI_w in paper pole space using correct APDL->paper mapping
    I_vec_k = zeros(6, 1);
    I_vec_k(apdl_to_paper_idx(k)) = 1;
    coil_data(k).KI_w = K_I * I_vec_k;

    fprintf('  Coil %d (%s): %d nodes in fitting region\n', ...
        k, cnst.apdl_to_paper{k}, coil_data(k).N);
end

fprintf('\n');

%% 2. Single-coil baseline: 2-param (ell, R_a) sphere model for initial guess
fprintf('===== Single-coil baseline: sphere (ell, R_a) =====\n');

% Use coil 1 to get initial ell estimate
I_vec_1 = zeros(6,1); I_vec_1(1) = 1;
p_wp_1 = [coil_data(1).px, coil_data(1).py, coil_data(1).pz];

cost_ell_fn = @(ell) ell_cost_sphere(ell, p_wp_1, coil_data(1).b_fem, I_vec_1, K_I, cnst);
ell_scan = linspace(400e-6, 2000e-6, 300);
cost_scan = arrayfun(cost_ell_fn, ell_scan);
[~, imin] = min(cost_scan);
ell_init = fminbnd(cost_ell_fn, max(100e-6, ell_scan(imin)-200e-6), ...
    ell_scan(imin)+200e-6, optimset('TolX', 1e-8, 'Display', 'off'));
[~, C_init] = ell_cost_sphere(ell_init, p_wp_1, coil_data(1).b_fem, I_vec_1, K_I, cnst);
R_a_init = cnst.N_c / (cnst.mu_0 * C_init);

fprintf('  Coil 1 baseline: ell = %.1f um, R_a = %.3e\n\n', ell_init*1e6, R_a_init);

%% 3. Build initial charge positions on sphere
alpha = cnst.alpha;
pos_sphere = zeros(3, 6);
for i = 1:6
    theta = cnst.pole_angles(i) * pi/180;
    z_sign = sign(cnst.pole_tip_z_wp(i));
    pos_sphere(:,i) = ell_init * [cos(theta)*sin(alpha); ...
                                   sin(theta)*sin(alpha); ...
                                   z_sign*cos(alpha)];
end

%% 4. Joint 6-coil fit: 18 shared positions, per-coil C_k (Variable Projection)
fprintf('===== Joint 6-coil fit: 19 params (18 positions + 6 C_k) =====\n');

fms_opts = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-22, ...
    'MaxIter', 50000, 'MaxFunEvals', 500000);

cost_fn = @(x) cost_joint_6coil(x, coil_data, cnst.k_m);

% 5 initial conditions
x0_sphere = pos_sphere(:);
rng(42);  x0_noise1 = pos_sphere(:) + 50e-6*randn(18,1);
rng(123); x0_noise2 = pos_sphere(:) + 100e-6*randn(18,1);
x0_far = 1.5 * pos_sphere(:);
x0_near = 0.7 * pos_sphere(:);

inits = {x0_sphere, x0_noise1, x0_noise2, x0_far, x0_near};
init_names = {'Sphere(ell=835)', '+50um noise', '+100um noise', '1.5x sphere', '0.7x sphere'};

results_J = struct();

fprintf('\n%-22s  %12s  %10s\n', 'Initial', 'Cost', 'Mean Err[%]');
fprintf('%s\n', repmat('-', 1, 50));

for trial = 1:length(inits)
    [x_opt, cost_opt] = fminsearch(cost_fn, inits{trial}, fms_opts);

    % Recover per-coil C_k and R_a
    [~, C_k_opt] = cost_joint_6coil(x_opt, coil_data, cnst.k_m);
    R_a_k = cnst.N_c ./ (cnst.mu_0 * C_k_opt);

    pos_opt = reshape(x_opt, 3, 6);

    % Evaluate per-coil errors
    err_per_coil = zeros(1, 6);
    for k = 1:6
        [bx_m, by_m, bz_m] = eval_joint_model(pos_opt, C_k_opt(k), ...
            coil_data(k).KI_w, cnst.k_m, ...
            [coil_data(k).px, coil_data(k).py, coil_data(k).pz]);
        delta = sqrt((bx_m - coil_data(k).bx).^2 + ...
                      (by_m - coil_data(k).by).^2 + ...
                      (bz_m - coil_data(k).bz).^2);
        err_per_coil(k) = mean(delta ./ coil_data(k).bmag);
    end

    results_J(trial).name = init_names{trial};
    results_J(trial).pos = pos_opt;
    results_J(trial).cost = cost_opt;
    results_J(trial).C_k = C_k_opt;
    results_J(trial).R_a_k = R_a_k;
    results_J(trial).err_per_coil = err_per_coil;
    results_J(trial).mean_err = mean(err_per_coil);

    fprintf('%-22s  %12.6e  %10.2f\n', ...
        init_names{trial}, cost_opt, mean(err_per_coil)*100);
end

%% 5. Select best trial
all_costs = [results_J.cost];
[~, best_idx] = min(all_costs);
best = results_J(best_idx);

fprintf('\n===== Best result: %s =====\n\n', best.name);

%% 6. Summary table: per-pole ell and direction
fprintf('--- Per-pole charge distance (ell) and direction deviation ---\n');
fprintf('%-5s %-6s  %8s  %8s  %6s  %28s\n', ...
    'Pole', 'Layer', 'ell[um]', 'Sph[um]', 'Dev[°]', 'Charge pos [um]');
fprintf('%s\n', repmat('-', 1, 80));

for i = 1:6
    ell_i = norm(best.pos(:,i));
    ell_sph = norm(pos_sphere(:,i));

    % Direction deviation from pole tip direction
    tip = [cnst.pole_tip_x(i); cnst.pole_tip_y(i); cnst.pole_tip_z_wp(i)];
    d_tip = tip / norm(tip);
    d_fit = best.pos(:,i) / ell_i;
    angle_dev = acosd(abs(dot(d_fit, d_tip)));

    layer = 'Upper';
    if cnst.pole_is_lower(i), layer = 'Lower'; end

    fprintf('P%d   %-6s  %8.1f  %8.1f  %6.2f  (%+7.1f, %+7.1f, %+7.1f)\n', ...
        i, layer, ell_i*1e6, ell_sph*1e6, angle_dev, ...
        best.pos(1,i)*1e6, best.pos(2,i)*1e6, best.pos(3,i)*1e6);
end

%% 7. Per-coil R_a table
fprintf('\n--- Per-coil R_a ---\n');
fprintf('%-5s %-6s  %12s  %10s\n', 'Coil', 'Pole', 'R_a [A/Wb]', 'C_k');
fprintf('%s\n', repmat('-', 1, 40));

for k = 1:6
    paper_idx = apdl_to_paper_idx(k);
    fprintf('  %d   %-4s   %12.3e  %10.4e\n', ...
        k, cnst.apdl_to_paper{k}, best.R_a_k(k), best.C_k(k));
end

% Lower/Upper R_a averages
lower_coils = find(cellfun(@(p) cnst.pole_is_lower(apdl_to_paper_idx(p)), num2cell(1:6)));
upper_coils = setdiff(1:6, lower_coils);
fprintf('\nLower avg R_a = %.3e, Upper avg R_a = %.3e\n', ...
    mean(best.R_a_k(lower_coils)), mean(best.R_a_k(upper_coils)));

%% 8. Per-coil error table
fprintf('\n--- Per-coil fitting error ---\n');
fprintf('%-5s %-6s  %10s\n', 'Coil', 'Pole', 'Vec Err[%]');
fprintf('%s\n', repmat('-', 1, 30));

for k = 1:6
    fprintf('  %d   %-4s   %10.2f\n', k, cnst.apdl_to_paper{k}, best.err_per_coil(k)*100);
end
fprintf('  Mean: %.2f%%\n', best.mean_err*100);

%% 9. Identifiability analysis
fprintf('\n===== Identifiability: position spread across 5 trials =====\n');

all_pos = cat(3, results_J.pos);  % 3 x 6 x ntrials

fprintf('\n--- Per-pole distance from WP center [um] ---\n');
fprintf('%-22s', 'Trial');
for i = 1:6, fprintf('  P%d      ', i); end
fprintf('\n%s\n', repmat('-', 1, 82));

for trial = 1:length(results_J)
    fprintf('%-22s', results_J(trial).name);
    for i = 1:6
        fprintf('  %7.1f ', norm(results_J(trial).pos(:,i))*1e6);
    end
    fprintf('\n');
end

fprintf('\n--- Spread (std) of charge position across trials [um] ---\n');
for i = 1:6
    pos_i = squeeze(all_pos(:,i,:)) * 1e6;  % 3 x ntrials
    fprintf('P%d: std(x)=%6.1f  std(y)=%6.1f  std(z)=%6.1f  total=%6.1f\n', ...
        i, std(pos_i(1,:)), std(pos_i(2,:)), std(pos_i(3,:)), ...
        sqrt(sum(std(pos_i, 0, 2).^2)));
end

%% 10. Cone check: are charges inside iron cones?
fprintf('\n--- Cone check: is each charge inside its pole cone? ---\n');
fprintf('%-5s  %8s  %8s  %8s  %s\n', 'Pole', 'ell[um]', 'perp[um]', 'Rcone', 'Inside?');
fprintf('%s\n', repmat('-', 1, 55));

for i = 1:6
    tip = [cnst.pole_tip_x(i); cnst.pole_tip_y(i); cnst.pole_tip_z_wp(i)];
    charge = best.pos(:,i);

    % Vector from tip to charge
    v = charge - tip;
    ax = cnst.pole_axis(:,i);
    s = dot(v, ax);  % projection along cone axis (positive = toward base)
    r_perp = norm(v - s*ax);  % perpendicular distance
    r_cone = cnst.POLE_TIP_R + s * (cnst.POLE_R - cnst.POLE_TIP_R) / cnst.POLE_CONE_LEN;

    inside = (s > 0) && (r_perp < r_cone) && (s < cnst.POLE_CONE_LEN);
    inside_str = 'NO';
    if inside, inside_str = 'YES'; end

    fprintf('P%d    %8.1f  %8.1f  %8.1f  %s\n', ...
        i, norm(charge)*1e6, r_perp*1e6, r_cone*1e6, inside_str);
end

%% 11. Direction deviation summary
fprintf('\n--- Direction deviation from pole tip [degrees] ---\n');
for i = 1:6
    tip = [cnst.pole_tip_x(i); cnst.pole_tip_y(i); cnst.pole_tip_z_wp(i)];
    d_tip = tip / norm(tip);
    d_fit = best.pos(:,i) / norm(best.pos(:,i));
    angle_dev = acosd(abs(dot(d_fit, d_tip)));
    fprintf('  P%d: %.2f deg\n', i, angle_dev);
end

%% 12. Save results
save_path = fullfile('..', 'data', 'joint_6coil_19param_fit.mat');
save(save_path, 'results_J', 'best', 'pos_sphere', 'ell_init', ...
    'coil_data', 'cnst', 'K_I', 'apdl_to_paper_idx');
fprintf('\nResults saved to %s\n', save_path);


%% ===== Local functions =====

function [cost, C_opt] = ell_cost_sphere(ell, p_wp, b_fem, I_vec, K_I, c)
% Sphere model cost for initial ell estimate (single coil)
    R_a_unit = c.N_c / c.mu_0;
    [bx, by, bz] = point_charge_model(p_wp, ell, R_a_unit, I_vec, K_I, c);
    b_unit = [bx; by; bz];
    C_opt = (b_unit' * b_fem) / (b_unit' * b_unit);
    cost = sum((C_opt * b_unit - b_fem).^2);
end

function [cost, C_k] = cost_joint_6coil(x, coil_data, k_m)
% Joint cost function: 18 shared positions, per-coil C_k solved analytically
%   x: 18x1 vector [x1,y1,z1,...,x6,y6,z6] - charge positions in paper pole order
%   coil_data: 1x6 struct array with px, py, pz, b_fem, KI_w per coil
%   k_m: mu_0/(4*pi)
%   Returns total SSE cost and 6x1 vector of per-coil C_k
    pos = reshape(x, 3, 6);
    cost = 0;
    C_k = zeros(6, 1);

    for k = 1:6
        px = coil_data(k).px;
        py = coil_data(k).py;
        pz = coil_data(k).pz;
        b_fem = coil_data(k).b_fem;
        KI_w = coil_data(k).KI_w;
        N = coil_data(k).N;

        % Compute unit B-field (C=1)
        bx_u = zeros(N, 1);
        by_u = zeros(N, 1);
        bz_u = zeros(N, 1);
        for i = 1:6
            dx = px - pos(1,i);
            dy = py - pos(2,i);
            dz = pz - pos(3,i);
            r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
            w = -KI_w(i);  % negative sign convention
            bx_u = bx_u + w * dx ./ r3;
            by_u = by_u + w * dy ./ r3;
            bz_u = bz_u + w * dz ./ r3;
        end
        b_unit = k_m * [bx_u; by_u; bz_u];

        % Analytical C_k (Variable Projection)
        C_k(k) = (b_unit' * b_fem) / (b_unit' * b_unit);
        cost = cost + sum((C_k(k) * b_unit - b_fem).^2);
    end
end

function [bx, by, bz] = eval_joint_model(pos, C, KI_w, k_m, p_wp)
% Evaluate model B-field for a single coil given shared positions and per-coil C
    N = size(p_wp, 1);
    bx = zeros(N, 1);
    by = zeros(N, 1);
    bz = zeros(N, 1);
    for i = 1:6
        dx = p_wp(:,1) - pos(1,i);
        dy = p_wp(:,2) - pos(2,i);
        dz = p_wp(:,3) - pos(3,i);
        r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
        w = -C * KI_w(i);  % negative sign convention
        bx = bx + w * dx ./ r3;
        by = by + w * dy ./ r3;
        bz = bz + w * dz ./ r3;
    end
    bx = k_m * bx;
    by = k_m * by;
    bz = k_m * bz;
end
