%% generate_figures_2_6.m — Fig 2.6(a,b): Point-charge model validation
% Corresponds to Long 2016 dissertation Fig. 2.6
%   (a) FEM vs point-charge model B-field quiver overlay in WP region
%   (b) Per-node fitting error scatter plot (normalized norms)
%
% Prerequisites:
%   - run post_extract_coil1.txt in ANSYS first
%   - run fit_charge_model.m to generate data/charge_model_fit.mat

clear; clc; close all;

%% ---- Unified style parameters ----
FONT_NAME   = 'Helvetica';
FONT_LABEL  = 12;
FONT_TITLE  = 13;
FONT_CB     = 11;
FONT_ANNOT  = 10;
LINE_MAIN   = 1.5;
DPI         = 300;

%% ---- Load constants and data ----
c = mt_constants();
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'coil1');
fig_dir     = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% FEM data
fprintf('Loading WP dataset...\n');
wp = import_ansys_data(results_dir, 'wp');
fprintf('  WP nodes: %d\n', length(wp.node_id));

% Fitted parameters
fit_file = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'charge_model_fit.mat');
fit = load(fit_file);
fprintf('Loaded fit: ell = %.0f um, R_a = %.2e (Fit A)\n', fit.fitA.ell*1e6, fit.fitA.R_a);

% WP-centered coordinates
x_wp = wp.x;
y_wp = wp.y;
z_wp = wp.z - c.SPH_OFST;
r_wp = sqrt(x_wp.^2 + y_wp.^2 + z_wp.^2);

% Iron exclusion
[air_mask, ~] = filter_iron_nodes(wp.x, wp.y, wp.z, c, struct('visualize', false));

% Use Fit A (100 um cube — matches dissertation approach)
ell = fit.fitA.ell;
R_a = fit.fitA.R_a;
K_I   = eye(6) - ones(6)/6;     % Eq. 2.8
I_vec = [1; 0; 0; 0; 0; 0];    % Coil 1 = P1

%% ======== Fig 2.6(a): FEM vs Model quiver overlay (80 um cube) ========
cube_half = 40e-6;   % ±40 um = 80 um cube (dissertation Fig 2.6 region)
mask_disp = air_mask & ...
	abs(x_wp) < cube_half & abs(y_wp) < cube_half & abs(z_wp) < cube_half;
fprintf('Fig 2.6(a): %d nodes in 80 um cube\n', sum(mask_disp));

p_disp = [x_wp(mask_disp), y_wp(mask_disp), z_wp(mask_disp)];
bx_fem = wp.bx(mask_disp);
by_fem = wp.by(mask_disp);
bz_fem = wp.bz(mask_disp);

[bx_mod, by_mod, bz_mod] = point_charge_model(p_disp, ell, R_a, I_vec, K_I, c);

fig_a = figure('Name', 'Fig 2.6(a)', 'Position', [50 50 800 700], 'Color', 'w');
set(fig_a, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

px = p_disp(:,1)*1e6;   % um
py = p_disp(:,2)*1e6;
pz = p_disp(:,3)*1e6;

% Normalize arrows for visibility
arrow_len = 6;  % um
b_fem_norm = sqrt(bx_fem.^2 + by_fem.^2 + bz_fem.^2);
b_fem_norm(b_fem_norm == 0) = 1;
b_mod_norm = sqrt(bx_mod.^2 + by_mod.^2 + bz_mod.^2);
b_mod_norm(b_mod_norm == 0) = 1;

hold on;
% FEM arrows (red)
q1 = quiver3(px, py, pz, ...
	bx_fem./b_fem_norm*arrow_len, by_fem./b_fem_norm*arrow_len, bz_fem./b_fem_norm*arrow_len, ...
	0, 'Color', [0.85 0.2 0.2], 'LineWidth', 1.3, 'MaxHeadSize', 0.3);

% Model arrows (blue, slight offset for visual separation)
offset = 1;  % 1 um offset
q2 = quiver3(px + offset, py + offset, pz + offset, ...
	bx_mod./b_mod_norm*arrow_len, by_mod./b_mod_norm*arrow_len, bz_mod./b_mod_norm*arrow_len, ...
	0, 'Color', [0.2 0.3 0.85], 'LineWidth', 1.0, 'MaxHeadSize', 0.3);
hold off;

legend([q1, q2], {'FEM', 'Point-charge model'}, ...
	'Location', 'northeast', 'FontSize', FONT_ANNOT, 'FontName', FONT_NAME);
axis equal; grid on; box on;
xlabel(['x [' char(956) 'm]'], 'FontSize', FONT_LABEL);
ylabel(['y [' char(956) 'm]'], 'FontSize', FONT_LABEL);
zlabel(['z [' char(956) 'm]'], 'FontSize', FONT_LABEL);
title(sprintf('(a) FEM vs point-charge model (%c = %.0f %cm)', ...
	char(8467), ell*1e6, char(956)), 'FontSize', FONT_TITLE);
view([-37.5, 30]);
lim = max(abs([px; py; pz])) * 1.3;
if isempty(lim) || lim == 0, lim = 50; end
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);
set(gca, 'Color', 'w');

exportgraphics(fig_a, fullfile(fig_dir, 'fig2_6a.png'), 'Resolution', DPI);
fprintf('Saved fig2_6a.png\n');

%% ======== Fig 2.6(b): Per-node fitting error scatter ========
% Evaluate over R < 400 um sphere (matching dissertation validation range)
R_err = 400e-6;
mask_err = air_mask & (r_wp < R_err);
fprintf('Fig 2.6(b): %d nodes in R < %.0f um sphere\n', sum(mask_err), R_err*1e6);

p_err = [x_wp(mask_err), y_wp(mask_err), z_wp(mask_err)];
bx_f  = wp.bx(mask_err);
by_f  = wp.by(mask_err);
bz_f  = wp.bz(mask_err);

[bx_m, by_m, bz_m] = point_charge_model(p_err, ell, R_a, I_vec, K_I, c);

% Per-node error: |B_model - B_FEM| / |B_FEM| * 100%
err_vec = sqrt((bx_m - bx_f).^2 + (by_m - by_f).^2 + (bz_m - bz_f).^2);
fem_mag = sqrt(bx_f.^2 + by_f.^2 + bz_f.^2);
valid   = fem_mag > 1e-6 * max(fem_mag);
pct_err = 100 * err_vec ./ fem_mag;
pct_err(~valid) = 0;

mean_err   = mean(pct_err(valid));
median_err = median(pct_err(valid));

fig_b = figure('Name', 'Fig 2.6(b)', 'Position', [100 100 900 500], 'Color', 'w');
set(fig_b, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

scatter(1:sum(mask_err), pct_err, 4, [0.3 0.3 0.8], 'filled', ...
	'MarkerFaceAlpha', 0.5);
hold on;
yline(mean_err, 'r-', sprintf('Mean = %.1f%%', mean_err), ...
	'LineWidth', LINE_MAIN, 'FontSize', FONT_ANNOT, ...
	'LabelVerticalAlignment', 'bottom');
yline(median_err, '--', sprintf('Median = %.1f%%', median_err), ...
	'Color', [0 0.6 0], 'LineWidth', 1.0, 'FontSize', FONT_ANNOT, ...
	'LabelVerticalAlignment', 'top');
hold off;

xlabel('Index of Points', 'FontSize', FONT_LABEL);
ylabel('Fitting Error [%]', 'FontSize', FONT_LABEL);
title(sprintf('(b) Per-node fitting error (%c = %.0f %cm, R < %.0f %cm)', ...
	char(8467), ell*1e6, char(956), R_err*1e6, char(956)), 'FontSize', FONT_TITLE);
grid on;
ylim([0 max(5, 1.1*max(pct_err(valid)))]);
set(gca, 'Color', 'w');

exportgraphics(fig_b, fullfile(fig_dir, 'fig2_6b.png'), 'Resolution', DPI);
fprintf('Saved fig2_6b.png\n');

fprintf('\nAll Fig 2.6 panels saved to %s\n', fig_dir);
