% VERIFY_SUPERPOSITION  Verify linear superposition of 6-coil ANSYS results
%   Compares B-field from simultaneous 6-coil excitation (CURR_ARRAY=[1,1,1,1,1,1])
%   against sum of 6 individual single-coil results.
%   Linear material (murx=280) => superposition should hold to solver tolerance.

%% Setup
clear; clc;
base_dir = fullfile(fileparts(mfilename('fullpath')), '..');
results_dir = fullfile(base_dir, 'results');
figures_dir = fullfile(base_dir, 'figures');

coil_names = {'coil1','coil2','coil3','coil4','coil5','coil6'};
n_coils = 6;

%% Load all6 simultaneous result
fprintf('=== Loading all6 simultaneous excitation data ===\n');
d_all6     = import_ansys_data(fullfile(results_dir, 'all6'), 'all', 'all6');
d_all6_wp  = import_ansys_data(fullfile(results_dir, 'all6'), 'wp',  'all6');

%% Load 6 single-coil results and sum
fprintf('\n=== Loading 6 single-coil results ===\n');
for i = 1:n_coils
    coil_dir = fullfile(results_dir, coil_names{i});
    d{i}    = import_ansys_data(coil_dir, 'all', coil_names{i});
    d_wp{i} = import_ansys_data(coil_dir, 'wp',  coil_names{i});
end

%% Verify node ID consistency (all dataset)
fprintf('\n=== Node ID consistency check (all nodes) ===\n');
ref_nodes = d_all6.node_id;
for i = 1:n_coils
    assert(isequal(d{i}.node_id, ref_nodes), ...
        'Node ID mismatch: %s vs all6', coil_names{i});
end
fprintf('PASS: All 6 single-coil node IDs match all6 (%d nodes)\n', length(ref_nodes));

%% Verify node ID consistency (WP dataset)
fprintf('\n=== Node ID consistency check (WP nodes) ===\n');
ref_nodes_wp = d_all6_wp.node_id;
for i = 1:n_coils
    assert(isequal(d_wp{i}.node_id, ref_nodes_wp), ...
        'WP Node ID mismatch: %s vs all6', coil_names{i});
end
fprintf('PASS: All 6 single-coil WP node IDs match all6 (%d nodes)\n', length(ref_nodes_wp));

%% Compute superposition sum (all nodes)
fprintf('\n=== Computing superposition (all nodes) ===\n');
bx_super = zeros(size(d_all6.bx));
by_super = zeros(size(d_all6.by));
bz_super = zeros(size(d_all6.bz));
for i = 1:n_coils
    bx_super = bx_super + d{i}.bx;
    by_super = by_super + d{i}.by;
    bz_super = bz_super + d{i}.bz;
end
bsum_super = sqrt(bx_super.^2 + by_super.^2 + bz_super.^2);

%% Compute superposition sum (WP nodes)
fprintf('=== Computing superposition (WP nodes) ===\n');
bx_super_wp = zeros(size(d_all6_wp.bx));
by_super_wp = zeros(size(d_all6_wp.by));
bz_super_wp = zeros(size(d_all6_wp.bz));
for i = 1:n_coils
    bx_super_wp = bx_super_wp + d_wp{i}.bx;
    by_super_wp = by_super_wp + d_wp{i}.by;
    bz_super_wp = bz_super_wp + d_wp{i}.bz;
end
bsum_super_wp = sqrt(bx_super_wp.^2 + by_super_wp.^2 + bz_super_wp.^2);

%% Error metrics (all nodes)
fprintf('\n=== Error metrics (ALL nodes, %d total) ===\n', length(ref_nodes));

% Component errors
delta_bx = d_all6.bx - bx_super;
delta_by = d_all6.by - by_super;
delta_bz = d_all6.bz - bz_super;

fprintf('Max |dBX| = %.4e T\n', max(abs(delta_bx)));
fprintf('Max |dBY| = %.4e T\n', max(abs(delta_by)));
fprintf('Max |dBZ| = %.4e T\n', max(abs(delta_bz)));

% Vector error
delta_vec = sqrt(delta_bx.^2 + delta_by.^2 + delta_bz.^2);
fprintf('Max |dB|  = %.4e T\n', max(delta_vec));
fprintf('Mean |dB| = %.4e T\n', mean(delta_vec));
fprintf('RMS |dB|  = %.4e T\n', sqrt(mean(delta_vec.^2)));

% Relative error (exclude boundary nodes with |B| < 1 uT)
b_ansys = sqrt(d_all6.bx.^2 + d_all6.by.^2 + d_all6.bz.^2);
mask = b_ansys > 1e-6;  % exclude |B| < 1 uT
rel_err = delta_vec(mask) ./ b_ansys(mask);
fprintf('\nRelative error (|B| > 1 uT, %d nodes):\n', sum(mask));
fprintf('Max  rel err = %.4e\n', max(rel_err));
fprintf('Mean rel err = %.4e\n', mean(rel_err));
fprintf('P99  rel err = %.4e\n', prctile(rel_err, 99));

%% Error metrics (WP nodes)
fprintf('\n=== Error metrics (WP nodes, %d total) ===\n', length(ref_nodes_wp));

delta_bx_wp = d_all6_wp.bx - bx_super_wp;
delta_by_wp = d_all6_wp.by - by_super_wp;
delta_bz_wp = d_all6_wp.bz - bz_super_wp;

fprintf('Max |dBX| = %.4e T\n', max(abs(delta_bx_wp)));
fprintf('Max |dBY| = %.4e T\n', max(abs(delta_by_wp)));
fprintf('Max |dBZ| = %.4e T\n', max(abs(delta_bz_wp)));

delta_vec_wp = sqrt(delta_bx_wp.^2 + delta_by_wp.^2 + delta_bz_wp.^2);
fprintf('Max |dB|  = %.4e T\n', max(delta_vec_wp));
fprintf('Mean |dB| = %.4e T\n', mean(delta_vec_wp));
fprintf('RMS |dB|  = %.4e T\n', sqrt(mean(delta_vec_wp.^2)));

b_ansys_wp = sqrt(d_all6_wp.bx.^2 + d_all6_wp.by.^2 + d_all6_wp.bz.^2);
mask_wp = b_ansys_wp > 1e-6;
rel_err_wp = delta_vec_wp(mask_wp) ./ b_ansys_wp(mask_wp);
fprintf('\nRelative error (WP, |B| > 1 uT, %d nodes):\n', sum(mask_wp));
fprintf('Max  rel err = %.4e\n', max(rel_err_wp));
fprintf('Mean rel err = %.4e\n', mean(rel_err_wp));

%% Pass/Fail
fprintf('\n========================================\n');
pass_threshold = 1e-6;
max_rel = max(rel_err);
if max_rel < pass_threshold
    fprintf('PASS: max relative error = %.2e < %.0e threshold\n', max_rel, pass_threshold);
else
    fprintf('FAIL: max relative error = %.2e >= %.0e threshold\n', max_rel, pass_threshold);
end
fprintf('========================================\n');

%% Diagnostic figures
c = mt_constants();

% Compute distance from WP center for each node (all dataset)
r_all = sqrt(d_all6.x.^2 + d_all6.y.^2 + (d_all6.z - c.SPH_OFST).^2);

fig = figure('Position', [100 100 1200 500]);

% Left: Error histogram
subplot(1,2,1);
histogram(log10(delta_vec(delta_vec > 0)), 50, 'FaceColor', [0.2 0.4 0.8]);
xlabel('log_{10}(|B_{ANSYS} - B_{super}|)  [T]');
ylabel('Node count');
title('Superposition error distribution (all nodes)');
grid on;

% Right: Error vs distance scatter
subplot(1,2,2);
scatter(r_all(mask)*1e3, rel_err, 1, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Distance from WP center [mm]');
ylabel('Relative error |dB|/|B|');
title(sprintf('Relative error vs distance (max = %.2e)', max_rel));
set(gca, 'YScale', 'log');
grid on;

sgtitle('Linear Superposition Verification: 6-coil simultaneous vs sum of singles');

% Save figure
saveas(fig, fullfile(figures_dir, 'verify_superposition.png'));
fprintf('\nFigure saved to figures/verify_superposition.png\n');
