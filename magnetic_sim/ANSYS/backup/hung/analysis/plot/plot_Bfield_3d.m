%% plot_Bfield_3d.m — 3D B-field vector plot for Hung hexapole
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
% Shows magnetic flux flow INSIDE iron structure only.

clear; clc; close all;

%% ---- Style ----
FONT_NAME  = 'Helvetica';
FONT_LABEL = 12;
FONT_TITLE = 13;
FONT_ANNOT = 10;
LINE_MAIN  = 1.5;
LINE_THIN  = 0.8;
DPI        = 400;

%% ---- Load data ----
c = mt_constants();
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'coil1');
fig_dir     = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

fprintf('Loading full dataset...\n');
d = import_ansys_data(results_dir, 'all');
fprintf('  Nodes: %d\n', length(d.node_id));

%% ---- Filter: only iron nodes (high B threshold) ----
% Iron regions: B > 5 mT (capture more of yoke/GP/pole flux)
iron_mask = d.bsum > 0.005;
fprintf('  Iron nodes (B > 20 mT): %d\n', sum(iron_mask));

x_fe = d.x(iron_mask) * 1e3;
y_fe = d.y(iron_mask) * 1e3;
z_fe = d.z(iron_mask) * 1e3;
bx_fe = d.bx(iron_mask);
by_fe = d.by(iron_mask);
bz_fe = d.bz(iron_mask);
bsum_fe = d.bsum(iron_mask);

%% ---- 3D spatial binning (iron only) ----
grid_cells = 50;  % finer binning to resolve GP/yoke/pole details
z_lo = -15;  z_hi = 50;
x_edges = linspace(-80, 80, grid_cells+1);
y_edges = linspace(-80, 80, grid_cells+1);
z_edges = linspace(z_lo, z_hi, grid_cells+1);

ix = discretize(x_fe, x_edges);
iy = discretize(y_fe, y_edges);
iz = discretize(z_fe, z_edges);
valid = ~isnan(ix) & ~isnan(iy) & ~isnan(iz);

cell_id = zeros(length(x_fe), 1);
cell_id(valid) = (ix(valid)-1)*grid_cells^2 + (iy(valid)-1)*grid_cells + iz(valid);

idx = zeros(grid_cells^3, 1);
n_found = 0;
unique_cells = unique(cell_id(cell_id > 0));
for k = 1:length(unique_cells)
    cell_nodes = find(cell_id == unique_cells(k));
    [~, best] = max(bsum_fe(cell_nodes));
    n_found = n_found + 1;
    idx(n_found) = cell_nodes(best);
end
idx = idx(1:n_found)';
fprintf('3D plot: %d arrows (iron only)\n', length(idx));

%% ---- Power-law scaling ----
bmag = bsum_fe(idx);
bvec = sqrt(bx_fe(idx).^2 + by_fe(idx).^2 + bz_fe(idx).^2);
bvec(bvec == 0) = 1e-10;

arrow_max = 5.0;
scale = arrow_max * bmag.^0.25 ./ bvec;
u = bx_fe(idx) .* scale;
v = by_fe(idx) .* scale;
w = bz_fe(idx) .* scale;

xp = x_fe(idx);
yp = y_fe(idx);
zp = z_fe(idx);

%% ---- Plot ----
fig = figure('Name', '3D B-field (iron)', 'Position', [50 50 1600 1200], 'Color', 'w');
set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

n_bins = 20;
bmax = max(bmag);
edges_b = linspace(0, bmax, n_bins+1);
cmap_q = turbo(n_bins);
lw_range = [0.4, 1.8];

hold on;
for k = 1:n_bins
    in_bin = bmag >= edges_b(k) & bmag < edges_b(k+1);
    if k == n_bins
        in_bin = in_bin | (bmag >= edges_b(end));
    end
    if any(in_bin)
        lw = lw_range(1) + (k-1)/(n_bins-1)*(lw_range(2)-lw_range(1));
        quiver3(xp(in_bin), yp(in_bin), zp(in_bin), ...
                u(in_bin), v(in_bin), w(in_bin), 0, ...
                'Color', cmap_q(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.3);
    end
end

% WP marker only
plot3(0, 0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);

% Yoke ring at correct Z height
% fc_yoke_zbot ~ 44mm for TILT_UP=35°, COIL_H=14mm
yoke_z_mm = 44.0;
theta_c = linspace(0, 2*pi, 100);
plot3(c.YOKE_RI*1e3*cos(theta_c), c.YOKE_RI*1e3*sin(theta_c), ...
    ones(size(theta_c))*yoke_z_mm, 'r--', 'LineWidth', LINE_MAIN);
plot3(c.YOKE_RO*1e3*cos(theta_c), c.YOKE_RO*1e3*sin(theta_c), ...
    ones(size(theta_c))*yoke_z_mm, 'r--', 'LineWidth', LINE_MAIN);

hold off;

% Colorbar
colormap(turbo);
clim([0 1.0]);
cb = colorbar;
ylabel(cb, 'Tesla', 'FontSize', 11);

% Axes
axis equal;
xlim([-80 80]); ylim([-80 80]); zlim([z_lo z_hi]);
xlabel('x [mm]', 'FontSize', FONT_LABEL);
ylabel('y [mm]', 'FontSize', FONT_LABEL);
zlabel('z [mm]', 'FontSize', FONT_LABEL);
title('3D B-field in iron structure (Coil1 excited)', 'FontSize', FONT_TITLE);

view(135, 25);
grid on;
set(gca, 'BoxStyle', 'full', 'Box', 'on');

exportgraphics(fig, fullfile(fig_dir, 'fig_3d_bfield.png'), 'Resolution', DPI);
fprintf('Saved fig_3d_bfield.png\n');
