%% fit_charge_model.m  -  Fit point-charge model parameters (ell, R_a) to ANSYS FEM data
%  Reproduces Long 2016 dissertation Section 2.2.3, Eq. 2.1-2.4
%  Uses Coil 1 (P1) data only; nominal K_I (Eq. 2.8, symmetric)
%
%  Two fitting ranges compared:
%    A) R < 100 um cube at center (dissertation approach, ~Fig 2.3b region)
%    B) R < 500 um sphere (full physical workspace)
%  Validation: 80 um cube at center (matching dissertation Fig 2.6)

%% 1. Load data
c = mt_constants();

data = import_ansys_data(fullfile('..', 'results', 'coil1'), 'wp', 'coil1');
fprintf('Loaded %d WP nodes from Coil 1\n', length(data.node_id));

%% 2. Convert to WP frame + iron exclusion
opts_filt = struct('visualize', false);
[air_mask, dbg] = filter_iron_nodes(data.x, data.y, data.z, c, opts_filt);

x_wp = data.x;
y_wp = data.y;
z_wp = data.z - c.SPH_OFST;
r_wp = sqrt(x_wp.^2 + y_wp.^2 + z_wp.^2);

bx_all = data.bx;
by_all = data.by;
bz_all = data.bz;

%% 3. Setup model (shared)
K_I = eye(6) - ones(6)/6;   % Eq. 2.8
I_vec = [1; 0; 0; 0; 0; 0]; % Coil 1 = P1

%% 4. Define fitting regions
% Region A: 100 um cube at center (dissertation approach)
cube_half_A = 50e-6;  % ±50 um
mask_A = air_mask & ...
    abs(x_wp) < cube_half_A & abs(y_wp) < cube_half_A & abs(z_wp) < cube_half_A;

% Region B: R < 500 um sphere (full workspace)
R_fit_B = 500e-6;
mask_B = air_mask & (r_wp < R_fit_B);

% Validation region: 80 um cube at center (Fig 2.6)
cube_half_V = 40e-6;  % ±40 um
mask_V = air_mask & ...
    abs(x_wp) < cube_half_V & abs(y_wp) < cube_half_V & abs(z_wp) < cube_half_V;

fprintf('\n--- Fitting regions ---\n');
fprintf('Region A (100 um cube): %d nodes\n', sum(mask_A));
fprintf('Region B (R < 500 um):  %d nodes\n', sum(mask_B));
fprintf('Validation (80 um cube): %d nodes\n', sum(mask_V));

%% 5. Fit both regions
fprintf('\n========== Fit A: 100 um cube ==========\n');
fitA = do_fit(x_wp, y_wp, z_wp, bx_all, by_all, bz_all, mask_A, I_vec, K_I, c);

fprintf('\n========== Fit B: R < 500 um sphere ==========\n');
fitB = do_fit(x_wp, y_wp, z_wp, bx_all, by_all, bz_all, mask_B, I_vec, K_I, c);

%% 6. Evaluate both fits on the 80 um validation cube
fprintf('\n========== Validation: 80 um cube ==========\n');
p_val = [x_wp(mask_V), y_wp(mask_V), z_wp(mask_V)];
bx_val = bx_all(mask_V);
by_val = by_all(mask_V);
bz_val = bz_all(mask_V);

% Fit A predictions in validation region
[bx_mA, by_mA, bz_mA] = point_charge_model(p_val, fitA.ell, fitA.R_a, I_vec, K_I, c);
errA = compute_error(bx_mA, by_mA, bz_mA, bx_val, by_val, bz_val);

% Fit B predictions in validation region
[bx_mB, by_mB, bz_mB] = point_charge_model(p_val, fitB.ell, fitB.R_a, I_vec, K_I, c);
errB = compute_error(bx_mB, by_mB, bz_mB, bx_val, by_val, bz_val);

fprintf('\n--- Comparison in 80 um cube (N = %d nodes) ---\n', sum(mask_V));
fprintf('%-25s  %12s  %12s\n', '', 'Fit A (cube)', 'Fit B (R<500)');
fprintf('%-25s  %12.1f  %12.1f\n', 'ell [um]', fitA.ell*1e6, fitB.ell*1e6);
fprintf('%-25s  %12.2e  %12.2e\n', 'R_a [A/Wb]', fitA.R_a, fitB.R_a);
fprintf('%-25s  %12.2f%%  %11.2f%%\n', 'Mean error', 100*errA.mean_rel, 100*errB.mean_rel);
fprintf('%-25s  %12.2f%%  %11.2f%%\n', 'Median error', 100*errA.median_rel, 100*errB.median_rel);
fprintf('%-25s  %12.2f%%  %11.2f%%\n', '95th percentile', 100*errA.p95_rel, 100*errB.p95_rel);
fprintf('%-25s  %12.2f%%  %11.2f%%\n', 'Max error', 100*errA.max_rel, 100*errB.max_rel);
fprintf('%-25s  %12.3e  %12.3e\n', 'RMS [T]', errA.rms_T, errB.rms_T);

%% 7. Dissertation comparison
fprintf('\n--- Dissertation Fig 2.6 reference ---\n');
fprintf('Dissertation: ell = 900 um, R_a = 6.3e8 A/Wb\n');
fprintf('Dissertation: error in 80 um cube < 2%%, mostly 0.5-1.5%%\n');

%% 8. Figures

% --- Fig 1: Cost landscape comparison ---
figure('Name', 'Cost Landscape', 'Position', [50 50 900 400]);
ell_scan = linspace(400e-6, 2000e-6, 300);

subplot(1,2,1);
plot_cost_landscape(ell_scan, x_wp, y_wp, z_wp, bx_all, by_all, bz_all, ...
    mask_A, I_vec, K_I, c, fitA.ell, 'Fit A: 100 um cube');

subplot(1,2,2);
plot_cost_landscape(ell_scan, x_wp, y_wp, z_wp, bx_all, by_all, bz_all, ...
    mask_B, I_vec, K_I, c, fitB.ell, 'Fit B: R < 500 um');

% --- Fig 2: Validation error comparison (bar chart) ---
figure('Name', 'Validation Error Comparison', 'Position', [100 100 700 450]);
metrics = [errA.mean_rel, errB.mean_rel; ...
           errA.median_rel, errB.median_rel; ...
           errA.p95_rel, errB.p95_rel; ...
           errA.max_rel, errB.max_rel] * 100;
bar(metrics);
set(gca, 'XTickLabel', {'Mean', 'Median', '95th', 'Max'});
ylabel('Relative Error [%]');
legend('Fit A (100 um cube)', 'Fit B (R < 500 um)', 'Location', 'northwest');
title('Validation Error in 80 um Cube');
grid on;

% --- Fig 3: Per-node error in validation cube (like Fig 2.6b) ---
figure('Name', 'Per-Node Error (Fig 2.6b style)', 'Position', [150 150 900 500]);

subplot(1,2,1);
plot(100*errA.rel_per_node, 'b.', 'MarkerSize', 4);
hold on;
yline(100*errA.mean_rel, 'r--', sprintf('mean=%.2f%%', 100*errA.mean_rel), 'LineWidth', 1.2);
hold off;
xlabel('Index of Points'); ylabel('Percent Fitting Error [%]');
title(sprintf('Fit A (100 um cube): \\ell=%.0f um', fitA.ell*1e6));
grid on; ylim([0 max(5, 1.2*100*errA.max_rel)]);

subplot(1,2,2);
plot(100*errB.rel_per_node, 'b.', 'MarkerSize', 4);
hold on;
yline(100*errB.mean_rel, 'r--', sprintf('mean=%.2f%%', 100*errB.mean_rel), 'LineWidth', 1.2);
hold off;
xlabel('Index of Points'); ylabel('Percent Fitting Error [%]');
title(sprintf('Fit B (R<500 um): \\ell=%.0f um', fitB.ell*1e6));
grid on; ylim([0 max(5, 1.2*100*errB.max_rel)]);

% --- Fig 4: FEM vs Model scatter for both fits ---
figure('Name', 'FEM vs Model Scatter', 'Position', [200 200 900 450]);

b_fem_val = sqrt(bx_val.^2 + by_val.^2 + bz_val.^2);
b_modA_val = sqrt(bx_mA.^2 + by_mA.^2 + bz_mA.^2);
b_modB_val = sqrt(bx_mB.^2 + by_mB.^2 + bz_mB.^2);
bmax = max(b_fem_val)*1e3;

subplot(1,2,1);
scatter(b_fem_val*1e3, b_modA_val*1e3, 3, 'b', '.'); hold on;
plot([0 bmax], [0 bmax], 'r-', 'LineWidth', 1.5); hold off;
xlabel('|B_{FEM}| [mT]'); ylabel('|B_{model}| [mT]');
title(sprintf('Fit A: \\ell=%.0f um, R_a=%.2e', fitA.ell*1e6, fitA.R_a));
axis equal; grid on; xlim([0 bmax]); ylim([0 bmax]);

subplot(1,2,2);
scatter(b_fem_val*1e3, b_modB_val*1e3, 3, 'b', '.'); hold on;
plot([0 bmax], [0 bmax], 'r-', 'LineWidth', 1.5); hold off;
xlabel('|B_{FEM}| [mT]'); ylabel('|B_{model}| [mT]');
title(sprintf('Fit B: \\ell=%.0f um, R_a=%.2e', fitB.ell*1e6, fitB.R_a));
axis equal; grid on; xlim([0 bmax]); ylim([0 bmax]);

% --- Fig 5: Spatial error map in validation cube ---
figure('Name', 'Spatial Error Map', 'Position', [250 250 900 450]);
x_v = x_wp(mask_V)*1e6; y_v = y_wp(mask_V)*1e6; z_v = z_wp(mask_V)*1e6;

subplot(1,2,1);
scatter3(x_v, y_v, z_v, 8, 100*errA.rel_per_node, 'filled');
colorbar; colormap(hot); caxis([0 min(5, 100*errA.p95_rel*1.5)]);
xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('z [\mum]');
title('Fit A: Spatial Error [%]'); axis equal; grid on; view(30,25);

subplot(1,2,2);
scatter3(x_v, y_v, z_v, 8, 100*errB.rel_per_node, 'filled');
colorbar; colormap(hot); caxis([0 min(5, 100*errB.p95_rel*1.5)]);
xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('z [\mum]');
title('Fit B: Spatial Error [%]'); axis equal; grid on; view(30,25);

%% 9. Save results
results.fitA = fitA;
results.fitB = fitB;
results.validation.N_nodes = sum(mask_V);
results.validation.cube_size = 80e-6;
results.validation.errA = errA;
results.validation.errB = errB;
results.K_I = K_I;

save(fullfile('..', 'data', 'charge_model_fit.mat'), '-struct', 'results');
fprintf('\nSaved to data/charge_model_fit.mat\n');

%% ===== Local functions =====

function fit = do_fit(x, y, z, bx, by, bz, mask, I_vec, K_I, c)
    p_wp = [x(mask), y(mask), z(mask)];
    b_fem = [bx(mask); by(mask); bz(mask)];
    N = sum(mask);

    cost_fn = @(ell) fit_ell_cost(ell, p_wp, b_fem, I_vec, K_I, c);

    % 1D scan
    ell_scan = linspace(400e-6, 2000e-6, 300);
    cost_scan = zeros(size(ell_scan));
    for k = 1:length(ell_scan)
        cost_scan(k) = cost_fn(ell_scan(k));
    end
    [~, imin] = min(cost_scan);
    ell_0 = ell_scan(imin);
    fprintf('  Scan minimum at ell = %.0f um\n', ell_0*1e6);

    % Refine
    ell_fit = fminbnd(cost_fn, max(100e-6, ell_0 - 200e-6), ell_0 + 200e-6, ...
        optimset('TolX', 1e-8, 'Display', 'off'));

    % Recover R_a
    [~, C_fit] = fit_ell_cost(ell_fit, p_wp, b_fem, I_vec, K_I, c);
    R_a_fit = c.N_c / (c.mu_0 * C_fit);

    % Model at fitted params
    [bx_m, by_m, bz_m] = point_charge_model(p_wp, ell_fit, R_a_fit, I_vec, K_I, c);
    err = compute_error(bx_m, by_m, bz_m, bx(mask), by(mask), bz(mask));

    fprintf('  ell = %.1f um, R_a = %.3e A/Wb\n', ell_fit*1e6, R_a_fit);
    fprintf('  N = %d, mean err = %.2f%%, median = %.2f%%, 95th = %.2f%%\n', ...
        N, 100*err.mean_rel, 100*err.median_rel, 100*err.p95_rel);

    fit.ell = ell_fit;
    fit.R_a = R_a_fit;
    fit.C   = C_fit;
    fit.N   = N;
    fit.err = err;
end

function err = compute_error(bx_m, by_m, bz_m, bx_f, by_f, bz_f)
    dx = bx_m - bx_f;
    dy = by_m - by_f;
    dz = bz_m - bz_f;
    err_mag = sqrt(dx.^2 + dy.^2 + dz.^2);
    fem_mag = sqrt(bx_f.^2 + by_f.^2 + bz_f.^2);

    % Relative error per node: |deltaB| / |B_FEM| (Eq from Fig 2.6c)
    valid = fem_mag > 1e-6 * max(fem_mag);
    rel = err_mag ./ fem_mag;
    rel(~valid) = 0;

    err.rel_per_node = rel;
    err.mean_rel   = mean(rel(valid));
    err.median_rel = median(rel(valid));
    err.p95_rel    = prctile(rel(valid), 95);
    err.max_rel    = max(rel(valid));
    err.rms_T      = sqrt(mean(err_mag.^2));
end

function [cost, C_opt] = fit_ell_cost(ell, p_wp, b_fem, I_vec, K_I, c)
    R_a_unit = c.N_c / c.mu_0;  % C = 1
    [bx, by, bz] = point_charge_model(p_wp, ell, R_a_unit, I_vec, K_I, c);
    b_unit = [bx; by; bz];

    C_opt = (b_unit' * b_fem) / (b_unit' * b_unit);

    residual = C_opt * b_unit - b_fem;
    cost = sum(residual.^2);
end

function plot_cost_landscape(ell_scan, x, y, z, bx, by, bz, mask, I_vec, K_I, c, ell_fit, ttl)
    p_wp = [x(mask), y(mask), z(mask)];
    b_fem = [bx(mask); by(mask); bz(mask)];

    cost_scan = zeros(size(ell_scan));
    R_a_scan  = zeros(size(ell_scan));
    for k = 1:length(ell_scan)
        [cost_scan(k), Ck] = fit_ell_cost(ell_scan(k), p_wp, b_fem, I_vec, K_I, c);
        R_a_scan(k) = c.N_c / (c.mu_0 * Ck);
    end

    yyaxis left;
    plot(ell_scan*1e6, cost_scan, 'b-', 'LineWidth', 1.5); hold on;
    [~, C_fit] = fit_ell_cost(ell_fit, p_wp, b_fem, I_vec, K_I, c);
    cost_fit = sum((C_fit * ...
        point_charge_b_vec(p_wp, ell_fit, I_vec, K_I, c) - b_fem).^2);
    plot(ell_fit*1e6, cost_fit, 'ro', 'MarkerSize', 8, 'LineWidth', 2); hold off;
    ylabel('SSR');

    yyaxis right;
    plot(ell_scan*1e6, R_a_scan*1e-8, 'r--', 'LineWidth', 1);
    ylabel('R_a [10^8 A/Wb]');

    xlabel('\ell [\mum]');
    title(ttl);
    xline(500, 'k:', 'phys. tip'); grid on;
end

function b_unit = point_charge_b_vec(p_wp, ell, I_vec, K_I, c)
    R_a_unit = c.N_c / c.mu_0;
    [bx, by, bz] = point_charge_model(p_wp, ell, R_a_unit, I_vec, K_I, c);
    b_unit = [bx; by; bz];
end
