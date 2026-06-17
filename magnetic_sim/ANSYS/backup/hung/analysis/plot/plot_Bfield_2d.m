%% plot_Bfield_2d.m — Hung hexapole B-field figures
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
% Generates 3 figures matching Long 2016 dissertation style:
%   (a) Fig 2.3a: Top view XY B-field vector distribution
%   (b) Fig 2.4a: |B| contour on XY plane (z=0), +/-300 um
%   (c) Fig 2.4b: |B| contour on XZ plane (y=0), +/-300 um
%
% Prerequisites: run post_export_data.txt in ANSYS first.

clear; clc; close all;

%% ---- Unified style parameters (matching Long2016) ----
FONT_NAME   = 'Helvetica';
FONT_LABEL  = 12;
FONT_TITLE  = 13;
FONT_CB     = 11;
FONT_ANNOT  = 10;
LINE_MAIN   = 1.5;
LINE_THIN   = 0.8;
DPI         = 300;
CMAP        = turbo(256);

%% ---- Load constants and data ----
c = mt_constants();
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'coil1');
fig_dir     = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

fprintf('Loading full dataset...\n');
d = import_ansys_data(results_dir, 'all');
fprintf('  Nodes: %d\n', length(d.node_id));

fprintf('Loading WP dataset...\n');
wp = import_ansys_data(results_dir, 'wp');
fprintf('  WP nodes: %d\n', length(wp.node_id));

%% ======== Fig 2.3(a): XY top view — B-field vector ========
% Style: power-law arrow scaling, binned linewidth, threshold filter
fig_a = figure('Name', 'Fig 2.3(a)', 'Position', [50 50 850 750], 'Color', 'w');
set(fig_a, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

% Spatial binning: 80x80 grid, pick highest |B| node per cell
grid_cells = 80;
x_mm = d.x * 1e3;
y_mm = d.y * 1e3;
x_edges = linspace(-80, 80, grid_cells+1);
y_edges = linspace(-80, 80, grid_cells+1);

ix_bin = discretize(x_mm, x_edges);
iy_bin = discretize(y_mm, y_edges);
valid_bin = ~isnan(ix_bin) & ~isnan(iy_bin);
cell_id = zeros(length(x_mm), 1);
cell_id(valid_bin) = (ix_bin(valid_bin)-1)*grid_cells + iy_bin(valid_bin);

idx = zeros(grid_cells*grid_cells, 1);
n_found = 0;
unique_cells = unique(cell_id(cell_id > 0));
for k = 1:length(unique_cells)
    cell_nodes = find(cell_id == unique_cells(k));
    [~, best] = max(d.bsum(cell_nodes));
    n_found = n_found + 1;
    idx(n_found) = cell_nodes(best);
end
idx = idx(1:n_found)';

% Filter: threshold > 0.1 mT
valid = d.bsum(idx) > 1e-4;
idx = idx(valid);
x_q = d.x(idx)*1e3;
y_q = d.y(idx)*1e3;
bmag_q = d.bsum(idx);

% Power-law arrow scaling: B^0.25 compresses dynamic range
arrow_max = 4.0;
bxy = sqrt(d.bx(idx).^2 + d.by(idx).^2);
bxy(bxy == 0) = 1e-10;
scale_factor = arrow_max * bmag_q.^0.25 ./ bxy;
bx_q = d.bx(idx) .* scale_factor;
by_q = d.by(idx) .* scale_factor;
fprintf('Fig 2.3(a): %d arrows after binning + threshold\n', length(idx));

% Color bins with varying linewidth
n_bins = 28;
bmax = max(bmag_q);
edges = linspace(0, bmax, n_bins + 1);
cmap_a = turbo(n_bins);
lw_range = [0.40, 1.75];

hold on;
for k = 1:n_bins
    in_bin = bmag_q >= edges(k) & bmag_q < edges(k+1);
    if k == n_bins
        in_bin = in_bin | (bmag_q >= edges(end));
    end
    if any(in_bin)
        lw = lw_range(1) + (k-1)/(n_bins-1) * (lw_range(2) - lw_range(1));
        quiver(x_q(in_bin), y_q(in_bin), bx_q(in_bin), by_q(in_bin), 0, ...
            'Color', cmap_a(k,:), 'LineWidth', lw, 'MaxHeadSize', 0.3);
    end
end

% Yoke outlines (red dashed)
theta_c = linspace(0, 2*pi, 200);
plot(c.YOKE_RI*1e3*cos(theta_c), c.YOKE_RI*1e3*sin(theta_c), ...
    'r--', 'LineWidth', LINE_THIN);
plot(c.YOKE_RO*1e3*cos(theta_c), c.YOKE_RO*1e3*sin(theta_c), ...
    'r--', 'LineWidth', LINE_THIN);

% Pole markers at GP/block positions (R~43mm, not at tip R~0.5mm)
% P1-P6 azimuth angles
pole_azim = [0, 180, 120, 300, 60, 240];
pole_labels = {'P1','P2','P3','P4','P5','P6'};
pole_R_mm = c.POLE_TOTAL_LEN * 1e3;  % ~43mm from center (approximate block position)
GP_R_mm = 4.0;  % guide post radius for marker circle

for k = 1:6
    ang = pole_azim(k) * pi/180;
    cx = pole_R_mm * cos(ang);
    cy = pole_R_mm * sin(ang);
    plot(cx + GP_R_mm*cos(theta_c), cy + GP_R_mm*sin(theta_c), ...
        'k-', 'LineWidth', LINE_MAIN);
    text(cx*1.15, cy*1.15, pole_labels{k}, ...
        'HorizontalAlignment', 'center', 'FontSize', FONT_ANNOT, ...
        'FontWeight', 'bold', 'FontName', FONT_NAME);
end

% WP center crosshair
plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
hold off;

% Colorbar (linear Tesla, 0 to 1)
colormap(turbo);
clim([0, 1.0]);
cb = colorbar;
ylabel(cb, 'Tesla', 'FontSize', FONT_CB);
set(cb, 'FontSize', FONT_CB);

axis equal; grid on;
xlim([-80 80]); ylim([-80 80]);
xlabel('x [mm]', 'FontSize', FONT_LABEL);
ylabel('y [mm]', 'FontSize', FONT_LABEL);
title('(a) Top view: B-field vector distribution (unit: Tesla)', 'FontSize', FONT_TITLE);

exportgraphics(fig_a, fullfile(fig_dir, 'fig2_3a.png'), 'Resolution', DPI);
fprintf('Saved fig2_3a.png\n');

%% ======== Fig 2.4(a): |B| contour on XY plane (z=0) ========
x_wp = wp.x;
y_wp = wp.y;
z_wp = wp.z;

slice_tol = 50e-6;
grid_n    = 200;
range_um  = 300;
n_levels  = 20;

mask_xy = abs(z_wp) < slice_tol;
fprintf('Fig 2.4(a): %d nodes in XY slice (|z| < %.0f um)\n', sum(mask_xy), slice_tol*1e6);

if sum(mask_xy) < 50
    slice_tol_xy = 100e-6;
    mask_xy = abs(z_wp) < slice_tol_xy;
    fprintf('  Expanded to %.0f um: %d nodes\n', slice_tol_xy*1e6, sum(mask_xy));
end

xq = linspace(-range_um, range_um, grid_n) * 1e-6;
yq = linspace(-range_um, range_um, grid_n) * 1e-6;
[Xg, Yg] = meshgrid(xq, yq);

F_xy  = scatteredInterpolant(x_wp(mask_xy), y_wp(mask_xy), wp.bsum(mask_xy), ...
    'natural', 'none');
Bg_xy = F_xy(Xg, Yg) * 1e3;   % mT

fig_b = figure('Name', 'Fig 2.4(a)', 'Position', [50 50 700 600], 'Color', 'w');
set(fig_b, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

max_xy_mT = max(Bg_xy(:), [], 'omitnan');
min_xy_mT = min(Bg_xy(:), [], 'omitnan');
contourf(Xg*1e6, Yg*1e6, Bg_xy, n_levels, 'LineColor', 'none');
colormap(CMAP);
clim([min_xy_mT, max_xy_mT]);
cb = colorbar;
ylabel(cb, '|B| [mT]', 'FontSize', FONT_CB);
% Force max value into colorbar ticks
cb_ticks = get(cb, 'Ticks');
cb_ticks = unique([cb_ticks, round(max_xy_mT, 1)]);
set(cb, 'Ticks', cb_ticks);
fprintf('  fig2_4a range = %.1f ~ %.1f mT\n', min_xy_mT, max_xy_mT);

axis equal; grid on;
xlim([-range_um range_um]); ylim([-range_um range_um]);
xlabel(['x [' char(956) 'm]'], 'FontSize', FONT_LABEL);
ylabel(['y [' char(956) 'm]'], 'FontSize', FONT_LABEL);
title('(a) magnetic flux density on x-y plane (z=0)', 'FontSize', FONT_TITLE);

exportgraphics(fig_b, fullfile(fig_dir, 'fig2_4a.png'), 'Resolution', DPI);
fprintf('Saved fig2_4a.png\n');

%% ======== Fig 2.4(b): |B| contour on XZ plane (y=0) ========
mask_xz = abs(y_wp) < slice_tol;
fprintf('Fig 2.4(b): %d nodes in XZ slice (|y| < %.0f um)\n', sum(mask_xz), slice_tol*1e6);

if sum(mask_xz) < 50
    slice_tol_xz = 100e-6;
    mask_xz = abs(y_wp) < slice_tol_xz;
    fprintf('  Expanded to %.0f um: %d nodes\n', slice_tol_xz*1e6, sum(mask_xz));
end

zq = linspace(-range_um, range_um, grid_n) * 1e-6;
[Xg2, Zg2] = meshgrid(xq, zq);

F_xz  = scatteredInterpolant(x_wp(mask_xz), z_wp(mask_xz), wp.bsum(mask_xz), ...
    'natural', 'none');
Bg_xz = F_xz(Xg2, Zg2) * 1e3;   % mT

% Clamp to reasonable range for comparison with Long2016
Bg_xz_clamped = min(Bg_xz, 35);  % cap at 35 mT (Long2016 max ~30 mT)

fig_c = figure('Name', 'Fig 2.4(b)', 'Position', [100 100 700 600], 'Color', 'w');
set(fig_c, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

max_xz_mT = max(Bg_xz_clamped(:), [], 'omitnan');
min_xz_mT = min(Bg_xz_clamped(:), [], 'omitnan');
contourf(Xg2*1e6, Zg2*1e6, Bg_xz_clamped, n_levels, 'LineColor', 'none');
colormap(CMAP);
clim([min_xz_mT, max_xz_mT]);
cb2 = colorbar;
ylabel(cb2, '|B| [mT]', 'FontSize', FONT_CB);
cb2_ticks = get(cb2, 'Ticks');
cb2_ticks = unique([cb2_ticks, round(max_xz_mT, 1)]);
set(cb2, 'Ticks', cb2_ticks);
fprintf('  fig2_4b range = %.1f ~ %.1f mT (clamped)\n', min_xz_mT, max_xz_mT);

axis equal; grid on;
xlim([-range_um range_um]); ylim([-range_um range_um]);
xlabel(['x [' char(956) 'm]'], 'FontSize', FONT_LABEL);
ylabel(['z [' char(956) 'm]'], 'FontSize', FONT_LABEL);
title('(b) magnetic flux density on x-z plane (y=0)', 'FontSize', FONT_TITLE);

exportgraphics(fig_c, fullfile(fig_dir, 'fig2_4b.png'), 'Resolution', DPI);
fprintf('Saved fig2_4b.png\n');

fprintf('\nAll figures saved to %s\n', fig_dir);
