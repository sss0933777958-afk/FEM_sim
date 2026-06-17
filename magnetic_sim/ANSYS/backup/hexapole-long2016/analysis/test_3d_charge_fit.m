%% test_3d_charge_fit.m  -  3D vector charge fitting for Coil 1
%
%  Model: B(p) = Q * (p - c) / |p - c|^3
%  For fixed c = (cx, cy, cz), Q solved analytically (normal equation).
%  fminsearch optimizes c using vector SSE as cost.
%
%  Outputs:
%    Fig 1: Dissertation-style (a) quiver overlay, (b) error scatter
%    Fig 2: xz-plane geometry — pole shape, tip, charge position
%    Verification: separated vs joint 4-param optimization equivalence

clear; clc; close all;

%% ---- Unified style parameters (match generate_figures_2_6.m) ----
FONT_NAME   = 'Helvetica';
FONT_LABEL  = 12;
FONT_TITLE  = 13;
FONT_CB     = 11;
FONT_ANNOT  = 10;
LINE_MAIN   = 1.5;
DPI         = 300;
fig_dir     = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');

%% 1. Setup and load data
cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

k = 1;  % APDL coil index
coilname = sprintf('coil%d', k);
data = import_ansys_data(fullfile('..', 'results', coilname), 'wp', coilname);
[air_mask, ~] = filter_iron_nodes(data.x, data.y, data.z, cnst, struct('visualize', false));
z_wp = data.z - cnst.SPH_OFST;
r_wp = sqrt(data.x.^2 + data.y.^2 + z_wp.^2);

paper_idx = apdl_to_paper_idx(k);
tip = [cnst.pole_tip_x(paper_idx); cnst.pole_tip_y(paper_idx); cnst.pole_tip_z_wp(paper_idx)];
d_pole = tip / norm(tip);
R_norm = norm(tip);

%% 2. Fit: 3D vector charge model
cube_half = 50e-6;
mask_fit = air_mask & ...
    abs(data.x) < cube_half & abs(data.y) < cube_half & abs(z_wp) < cube_half;

px = data.x(mask_fit);  py = data.y(mask_fit);  pz = z_wp(mask_fit);
bx_f = data.bx(mask_fit); by_f = data.by(mask_fit); bz_f = data.bz(mask_fit);
bmag_f = sqrt(bx_f.^2 + by_f.^2 + bz_f.^2);
bvec_fem = [bx_f; by_f; bz_f];
N_fit = sum(mask_fit);

% 1D baseline for initial guess
b_scan = linspace(0, 800e-6, 2000);
cost_scan = zeros(size(b_scan));
for j = 1:length(b_scan)
    cp = (R_norm + b_scan(j)) * d_pole;
    r2 = (px - cp(1)).^2 + (py - cp(2)).^2 + (pz - cp(3)).^2;
    m = 1 ./ r2;
    a2 = (m' * bmag_f) / (m' * m);
    cost_scan(j) = sum((a2 * m - bmag_f).^2);
end
[~, imin] = min(cost_scan);
charge_1d = (R_norm + b_scan(imin)) * d_pole;

% 3D vector fit via fminsearch (variable projection)
cost_fn = @(pos) vector_charge_cost(pos, px, py, pz, bvec_fem);
opts = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-22, 'MaxIter', 10000);
[c1, f1] = fminsearch(cost_fn, charge_1d(:)', opts);
[c2, f2] = fminsearch(cost_fn, -charge_1d(:)', opts);
if f2 < f1, charge_opt = c2; else, charge_opt = c1; end

% Recover Q
[dx, dy, dz] = deal(px - charge_opt(1), py - charge_opt(2), pz - charge_opt(3));
r2 = dx.^2 + dy.^2 + dz.^2;  r = sqrt(r2);  r3 = r2 .* r;
mvec = [dx ./ r3; dy ./ r3; dz ./ r3];
Q_opt = (mvec' * bvec_fem) / (mvec' * mvec);

fprintf('============================================================\n');
fprintf('  3D Vector Charge Fit — Coil 1 (P1)\n');
fprintf('  Fitting: 100 um cube, N = %d\n', N_fit);
fprintf('============================================================\n');
fprintf('Tip:    [%+.1f, %+.1f, %+.1f] um   |tip| = %.1f um\n', tip*1e6, R_norm*1e6);
fprintf('Charge: [%+.1f, %+.1f, %+.1f] um   |c|   = %.1f um\n', charge_opt*1e6, norm(charge_opt)*1e6);
fprintf('Q = %.4e\n', Q_opt);
fprintf('Deviation from pole axis: %.2f deg\n', ...
    acosd(abs(dot(charge_opt(:)/norm(charge_opt), d_pole))));

%% 3. Verify: separated vs joint 4-param optimization
% Joint: fminsearch over [Q, cx, cy, cz] simultaneously
cost_joint = @(params) joint_cost(params, px, py, pz, bvec_fem);
x0_joint = [Q_opt * 0.8, charge_1d(:)'];  % perturbed initial guess
opts_j = optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-24, 'MaxIter', 20000);
[params_joint, fval_joint] = fminsearch(cost_joint, x0_joint, opts_j);
Q_joint = params_joint(1);
charge_joint = params_joint(2:4);

fprintf('\n--- Verification: separated vs joint 4-param ---\n');
fprintf('Separated: c = [%.2f, %.2f, %.2f] um, Q = %.6e, cost = %.6e\n', ...
    charge_opt*1e6, Q_opt, cost_fn(charge_opt));
fprintf('Joint:     c = [%.2f, %.2f, %.2f] um, Q = %.6e, cost = %.6e\n', ...
    charge_joint*1e6, Q_joint, cost_joint(params_joint));
fprintf('Difference in cost: %.2e (ratio: %.6f)\n', ...
    abs(cost_fn(charge_opt) - fval_joint), fval_joint / cost_fn(charge_opt));

%% 4. Evaluate on larger region for Fig 2.6(b)
R_err = 400e-6;
mask_err = air_mask & (r_wp < R_err);
N_err = sum(mask_err);

p_err_x = data.x(mask_err);
p_err_y = data.y(mask_err);
p_err_z = z_wp(mask_err);
bx_fe = data.bx(mask_err);
by_fe = data.by(mask_err);
bz_fe = data.bz(mask_err);
bmag_fe = sqrt(bx_fe.^2 + by_fe.^2 + bz_fe.^2);

% Predict
dx_e = p_err_x - charge_opt(1);
dy_e = p_err_y - charge_opt(2);
dz_e = p_err_z - charge_opt(3);
r2_e = dx_e.^2 + dy_e.^2 + dz_e.^2;
r_e  = sqrt(r2_e);
r3_e = r2_e .* r_e;
bx_mod = Q_opt * dx_e ./ r3_e;
by_mod = Q_opt * dy_e ./ r3_e;
bz_mod = Q_opt * dz_e ./ r3_e;

% Per-node vector error
err_vec = sqrt((bx_mod - bx_fe).^2 + (by_mod - by_fe).^2 + (bz_mod - bz_fe).^2);
valid = bmag_fe > 1e-6 * max(bmag_fe);
pct_err = 100 * err_vec ./ bmag_fe;
pct_err(~valid) = 0;

mean_err = mean(pct_err(valid));
median_err = median(pct_err(valid));
fprintf('\nError on R < 400 um sphere (%d nodes):\n', N_err);
fprintf('  Mean:   %.2f%%\n', mean_err);
fprintf('  Median: %.2f%%\n', median_err);
fprintf('  Max:    %.2f%%\n', max(pct_err(valid)));

%% 5. Also evaluate on 100 um cube (fitting region)
bx_fit = Q_opt * dx ./ r3;
by_fit = Q_opt * dy ./ r3;
bz_fit = Q_opt * dz ./ r3;
err_fit = sqrt((bx_fit - bx_f).^2 + (by_fit - by_f).^2 + (bz_fit - bz_f).^2);
pct_fit = 100 * err_fit ./ bmag_f;
fprintf('\nError on 100 um cube (%d nodes):\n', N_fit);
fprintf('  Mean:   %.2f%%\n', mean(pct_fit));
fprintf('  Median: %.2f%%\n', median(pct_fit));
fprintf('  Max:    %.2f%%\n', max(pct_fit));

%% ======== Fig 2.6(a): FEM vs Model quiver overlay ========
cube_disp = 40e-6;  % 80 um cube for display
mask_disp = air_mask & ...
    abs(data.x) < cube_disp & abs(data.y) < cube_disp & abs(z_wp) < cube_disp;

p_dx = data.x(mask_disp)*1e6;
p_dy = data.y(mask_disp)*1e6;
p_dz = z_wp(mask_disp)*1e6;

bx_fd = data.bx(mask_disp);
by_fd = data.by(mask_disp);
bz_fd = data.bz(mask_disp);

% Model prediction at display points
dd_x = data.x(mask_disp) - charge_opt(1);
dd_y = data.y(mask_disp) - charge_opt(2);
dd_z = z_wp(mask_disp)   - charge_opt(3);
rr3 = (dd_x.^2 + dd_y.^2 + dd_z.^2).^(3/2);
bx_md = Q_opt * dd_x ./ rr3;
by_md = Q_opt * dd_y ./ rr3;
bz_md = Q_opt * dd_z ./ rr3;

% Normalize arrows
arrow_len = 6;
bn_fem = sqrt(bx_fd.^2 + by_fd.^2 + bz_fd.^2);  bn_fem(bn_fem==0) = 1;
bn_mod = sqrt(bx_md.^2 + by_md.^2 + bz_md.^2);  bn_mod(bn_mod==0) = 1;

fig_a = figure('Name', 'Fig 2.6(a)', 'Position', [50 50 800 700], 'Color', 'w');
set(fig_a, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);
hold on;
q1 = quiver3(p_dx, p_dy, p_dz, ...
    bx_fd./bn_fem*arrow_len, by_fd./bn_fem*arrow_len, bz_fd./bn_fem*arrow_len, ...
    0, 'Color', [0.85 0.2 0.2], 'LineWidth', 1.3, 'MaxHeadSize', 0.3);
ofs = 1;
q2 = quiver3(p_dx+ofs, p_dy+ofs, p_dz+ofs, ...
    bx_md./bn_mod*arrow_len, by_md./bn_mod*arrow_len, bz_md./bn_mod*arrow_len, ...
    0, 'Color', [0.2 0.3 0.85], 'LineWidth', 1.0, 'MaxHeadSize', 0.3);
hold off;
legend([q1, q2], {'Flux Density(FEM)', 'Flux density(fit)'}, ...
    'Location', 'northeast', 'FontSize', FONT_ANNOT, 'FontName', FONT_NAME);
axis equal; grid on; box on;
xlabel(['x(' char(956) 'm)'], 'FontSize', FONT_LABEL);
ylabel(['y(' char(956) 'm)'], 'FontSize', FONT_LABEL);
zlabel(['z(' char(956) 'm)'], 'FontSize', FONT_LABEL);
title('(a)', 'FontSize', FONT_TITLE);
view([-37.5, 30]);
lm = max(abs([p_dx; p_dy; p_dz])) * 1.3;
xlim([-lm lm]); ylim([-lm lm]); zlim([-lm lm]);
set(gca, 'Color', 'w');
exportgraphics(fig_a, fullfile(fig_dir, 'fig2_6a_vecfit.png'), 'Resolution', DPI);
fprintf('\nSaved fig2_6a_vecfit.png\n');

%% ======== Fig 2.6(b): Error scatter (100 um cube — fitting region) ========
mean_fit = mean(pct_fit);
median_fit = median(pct_fit);

fig_b = figure('Name', 'Fig 2.6(b)', 'Position', [100 100 900 500], 'Color', 'w');
set(fig_b, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

scatter(1:N_fit, pct_fit, 4, [0.3 0.3 0.8], 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
yline(mean_fit, 'r-', sprintf('Mean = %.2f%%', mean_fit), ...
    'LineWidth', LINE_MAIN, 'FontSize', FONT_ANNOT, ...
    'LabelVerticalAlignment', 'bottom');
yline(median_fit, '--', sprintf('Median = %.2f%%', median_fit), ...
    'Color', [0 0.6 0], 'LineWidth', 1.0, 'FontSize', FONT_ANNOT, ...
    'LabelVerticalAlignment', 'top');
hold off;
xlabel('Index of Points', 'FontSize', FONT_LABEL);
ylabel('Percent fitting error', 'FontSize', FONT_LABEL);
title('(b)', 'FontSize', FONT_TITLE);
grid on;
ylim([0 max(5, 1.1*max(pct_fit))]);
set(gca, 'Color', 'w');
exportgraphics(fig_b, fullfile(fig_dir, 'fig2_6b_vecfit.png'), 'Resolution', DPI);
fprintf('Saved fig2_6b_vecfit.png\n');

%% ======== Fig 3: xz-plane geometry diagram ========
% Shows: pole cone shape (P1 lower, P2 upper), WP region, tip, charge
fig_c = figure('Name', 'xz Geometry', 'Position', [150 50 1000 800], 'Color', 'w');
set(fig_c, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);
hold on;

% -- Draw P1 (lower, angle=0, horizontal) --
cone_len = cnst.POLE_CONE_LEN;
tip_r = cnst.POLE_TIP_R;
base_r = cnst.POLE_R;

% P1 lower: horizontal, milled flat (upper edge flat at z=tip_z)
p1_tip_x = cnst.R_norm_xy;
p1_tip_z = -cnst.R_norm_z;
x_cone = [p1_tip_x, p1_tip_x + cone_len] * 1e3;
z_upper = [p1_tip_z, p1_tip_z] * 1e3;
z_lower = [p1_tip_z - tip_r, p1_tip_z - base_r] * 1e3;
fill([x_cone, fliplr(x_cone)], [z_upper, fliplr(z_lower)], ...
    [0.75 0.75 0.75], 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.2);

% -- Draw P2 (upper, angle=180, inclined upward) --
p2_tip_x = cnst.pole_tip_x(2);
p2_tip_z = cnst.pole_tip_z_wp(2);
inc = cnst.upper_incline;
t_pts = linspace(0, cone_len, 50);
cx_center = p2_tip_x - cos(inc) * t_pts;
cz_center = p2_tip_z + sin(inc) * t_pts;
r_at_t = tip_r + (base_r - tip_r) * t_pts / cone_len;
perp_x = sin(inc);  perp_z = cos(inc);
x_top = (cx_center + r_at_t * perp_x) * 1e3;
z_top = (cz_center + r_at_t * perp_z) * 1e3;
x_bot = (cx_center - r_at_t * perp_x) * 1e3;
z_bot = (cz_center - r_at_t * perp_z) * 1e3;
fill([x_top, fliplr(x_bot)], [z_top, fliplr(z_bot)], ...
    [0.75 0.75 0.75], 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.2);

% -- Workspace sphere (R = 500 um) --
theta_ws = linspace(0, 2*pi, 200);
h_ws = plot(cnst.R_norm * cos(theta_ws) * 1e3, cnst.R_norm * sin(theta_ws) * 1e3, ...
    'b--', 'LineWidth', 1.0);

% -- WP center --
h_wp = plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2.5);

% -- 100 um cube (fitting region) --
cb = cube_half * 1e3;
rectangle('Position', [-cb, -cb, 2*cb, 2*cb], ...
    'EdgeColor', [0 0.6 0], 'LineStyle', '--', 'LineWidth', 1.0);

% -- Tip positions --
h_tip1 = plot(tip(1)*1e3, tip(3)*1e3, 'r^', 'MarkerSize', 13, 'MarkerFaceColor', 'r');
plot(cnst.pole_tip_x(2)*1e3, cnst.pole_tip_z_wp(2)*1e3, 'r^', ...
    'MarkerSize', 13, 'MarkerFaceColor', 'r');

% -- Line from WP center through charge (ell vector) --
plot([0, charge_opt(1)*1e3], [0, charge_opt(3)*1e3], ...
    'b-', 'LineWidth', 1.2);
% Dashed extension to show pole axis
ext = 1.2;
plot([0, d_pole(1)*ext], [0, d_pole(3)*ext], ...
    'k:', 'LineWidth', 0.8);

% -- Charge position (3D vector fit) --
h_ch3d = plot(charge_opt(1)*1e3, charge_opt(3)*1e3, 'bp', ...
    'MarkerSize', 15, 'MarkerFaceColor', [0.2 0.4 0.9], 'LineWidth', 1.5);

% -- 1D charge position --
h_ch1d = plot(charge_1d(1)*1e3, charge_1d(3)*1e3, 'gs', ...
    'MarkerSize', 11, 'MarkerFaceColor', [0.3 0.8 0.3], 'LineWidth', 1.0);

% -- Annotations with better placement --
% P1 tip label (above)
text(tip(1)*1e3 - 0.10, tip(3)*1e3 + 0.07, ...
    sprintf('P1 tip (%.0f, %.0f) %cm', tip(1)*1e6, tip(3)*1e6, char(956)), ...
    'FontSize', FONT_ANNOT, 'FontName', FONT_NAME, 'Color', [0.7 0 0], 'FontWeight', 'bold');
% P2 tip label
text(cnst.pole_tip_x(2)*1e3 + 0.03, cnst.pole_tip_z_wp(2)*1e3 + 0.06, ...
    sprintf('P2 tip (%.0f, %.0f) %cm', cnst.pole_tip_x(2)*1e6, cnst.pole_tip_z_wp(2)*1e6, char(956)), ...
    'FontSize', FONT_ANNOT, 'FontName', FONT_NAME, 'Color', [0.7 0 0], 'FontWeight', 'bold');

% 3D charge label (below-left)
text(charge_opt(1)*1e3 - 0.45, charge_opt(3)*1e3 - 0.07, ...
    sprintf('3D charge (%.0f, %.0f, %.0f) %cm', charge_opt*1e6, char(956)), ...
    'FontSize', FONT_ANNOT, 'FontName', FONT_NAME, 'Color', [0.1 0.2 0.8], 'FontWeight', 'bold');

% 1D charge label (above-right)
text(charge_1d(1)*1e3 + 0.04, charge_1d(3)*1e3 + 0.05, ...
    sprintf('1D charge (%.0f, %.0f, %.0f) %cm', charge_1d*1e6, char(956)), ...
    'FontSize', FONT_ANNOT, 'FontName', FONT_NAME, 'Color', [0.1 0.5 0.1], 'FontWeight', 'bold');

% WP center label
text(0.03, 0.06, 'WP center', 'FontSize', FONT_ANNOT, 'FontName', FONT_NAME, 'FontWeight', 'bold');

% ell annotation along the line
ell_mid = charge_opt(:) * 0.45;
text(ell_mid(1)*1e3 - 0.15, ell_mid(3)*1e3 + 0.04, ...
    sprintf('%c = %.0f %cm', char(8467), norm(charge_opt)*1e6, char(956)), ...
    'FontSize', FONT_ANNOT+1, 'FontName', FONT_NAME, 'Color', 'b', 'FontWeight', 'bold');

% Legend
legend([h_tip1, h_ch3d, h_ch1d, h_ws], ...
    {'Pole tip', '3D vector fit charge', '1D scalar fit charge', ...
    sprintf('Workspace sphere (R=%.0f %cm)', cnst.R_norm*1e6, char(956))}, ...
    'Location', 'southeast', 'FontSize', FONT_ANNOT, 'FontName', FONT_NAME);

hold off;
xlabel('x (mm)', 'FontSize', FONT_LABEL);
ylabel('z (mm)', 'FontSize', FONT_LABEL);
title('(c) xz-plane: pole geometry, tip & charge positions (Coil 1 / P1)', ...
    'FontSize', FONT_TITLE);
axis equal; grid on; box on;
xlim([-1.2 1.5]);
ylim([-0.9 0.9]);
set(gca, 'Color', 'w');

exportgraphics(fig_c, fullfile(fig_dir, 'fig_xz_geometry_coil1.png'), 'Resolution', DPI);
fprintf('Saved fig_xz_geometry_coil1.png\n');

%% Local functions
function cost = vector_charge_cost(pos, px, py, pz, bvec_fem)
    dx = px - pos(1);  dy = py - pos(2);  dz = pz - pos(3);
    r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
    mvec = [dx ./ r3; dy ./ r3; dz ./ r3];
    Q = (mvec' * bvec_fem) / (mvec' * mvec);
    residual = Q * mvec - bvec_fem;
    cost = sum(residual.^2);
end

function cost = joint_cost(params, px, py, pz, bvec_fem)
    Q = params(1);
    pos = params(2:4);
    dx = px - pos(1);  dy = py - pos(2);  dz = pz - pos(3);
    r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
    bpred = Q * [dx ./ r3; dy ./ r3; dz ./ r3];
    cost = sum((bpred - bvec_fem).^2);
end
