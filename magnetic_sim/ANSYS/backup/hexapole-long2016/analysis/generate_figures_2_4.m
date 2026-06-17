%% generate_figures_2_4.m — Fig 2.4(a,b): |B| contour plots
% Corresponds to Long 2016 dissertation Fig. 2.4
%   (a) Horizontal contour: |B| in XY plane at z_wp = 0
%   (b) Vertical contour:   |B| in XZ plane at y_wp = 0
%
% Dissertation uses +/-300 um range with Gauss units.
% We use +/-300 um range with mT units (1 mT = 10 Gauss).
% Each panel has its own color scale (not clamped).
%
% Prerequisites: run post_extract_coil1.txt in ANSYS first.

clear; clc; close all;

%% ---- Unified style parameters ----
FONT_NAME   = 'Helvetica';
FONT_LABEL  = 12;
FONT_TITLE  = 13;
FONT_CB     = 11;
FONT_ANNOT  = 10;
LINE_MAIN   = 1.5;
DPI         = 300;
CMAP        = turbo(256);

%% ---- Load constants and data ----
c = mt_constants();
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'coil1');
fig_dir     = fullfile(fileparts(mfilename('fullpath')), '..', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

fprintf('Loading WP dataset...\n');
wp = import_ansys_data(results_dir, 'wp');
fprintf('  WP nodes: %d\n', length(wp.node_id));

% WP-centered coordinates
x_wp = wp.x;
y_wp = wp.y;
z_wp = wp.z - c.SPH_OFST;

%% ---- Common interpolation parameters ----
slice_tol = 50e-6;     % half-thickness for slice selection [m]
grid_n    = 200;        % interpolation grid resolution
range_um  = 300;        % plot range: +/- 300 um (matches dissertation)
n_levels  = 20;         % number of contour levels

%% ======== Fig 2.4(a): Horizontal contour — XY at z_wp = 0 ========
mask_xy = abs(z_wp) < slice_tol;
fprintf('Fig 2.4(a): %d nodes in XY slice (|z_wp| < %.0f um)\n', ...
	sum(mask_xy), slice_tol*1e6);

if sum(mask_xy) < 50
	slice_tol_xy = 100e-6;
	mask_xy = abs(z_wp) < slice_tol_xy;
	fprintf('  Expanded to %.0f um: %d nodes\n', slice_tol_xy*1e6, sum(mask_xy));
end

xq = linspace(-range_um, range_um, grid_n) * 1e-6;   % [m]
yq = linspace(-range_um, range_um, grid_n) * 1e-6;
[Xg, Yg] = meshgrid(xq, yq);

F_xy  = scatteredInterpolant(x_wp(mask_xy), y_wp(mask_xy), wp.bsum(mask_xy), ...
	'natural', 'none');
Bg_xy = F_xy(Xg, Yg) * 1e3;   % mT

fig_a = figure('Name', 'Fig 2.4(a)', 'Position', [50 50 700 600], 'Color', 'w');
set(fig_a, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

contourf(Xg*1e6, Yg*1e6, Bg_xy, n_levels, 'LineColor', 'none');
colormap(CMAP);
cb = colorbar;
ylabel(cb, '|B| [mT]', 'FontSize', FONT_CB);

axis equal; grid on;
xlim([-range_um range_um]); ylim([-range_um range_um]);
xlabel(['x [' char(956) 'm]'], 'FontSize', FONT_LABEL);
ylabel(['y [' char(956) 'm]'], 'FontSize', FONT_LABEL);
title('(a) magnetic flux density on x-y plane (z=0)', 'FontSize', FONT_TITLE);

exportgraphics(fig_a, fullfile(fig_dir, 'fig2_4a.png'), 'Resolution', DPI);
fprintf('Saved fig2_4a.png\n');

%% ======== Fig 2.4(b): Vertical contour — XZ at y_wp = 0 ========
mask_xz = abs(y_wp) < slice_tol;
fprintf('Fig 2.4(b): %d nodes in XZ slice (|y_wp| < %.0f um)\n', ...
	sum(mask_xz), slice_tol*1e6);

if sum(mask_xz) < 50
	slice_tol_xz = 100e-6;
	mask_xz = abs(y_wp) < slice_tol_xz;
	fprintf('  Expanded to %.0f um: %d nodes\n', slice_tol_xz*1e6, sum(mask_xz));
end

zq = linspace(-range_um, range_um, grid_n) * 1e-6;
[Xg2, Zg2] = meshgrid(xq, zq);

F_xz  = scatteredInterpolant(x_wp(mask_xz), z_wp(mask_xz), wp.bsum(mask_xz), ...
	'natural', 'none');
Bg_xz = F_xz(Xg2, Zg2) * 1e3;  % mT

fig_b = figure('Name', 'Fig 2.4(b)', 'Position', [100 100 700 600], 'Color', 'w');
set(fig_b, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

contourf(Xg2*1e6, Zg2*1e6, Bg_xz, n_levels, 'LineColor', 'none');
colormap(CMAP);
cb2 = colorbar;
ylabel(cb2, '|B| [mT]', 'FontSize', FONT_CB);

axis equal; grid on;
xlim([-range_um range_um]); ylim([-range_um range_um]);
xlabel(['x [' char(956) 'm]'], 'FontSize', FONT_LABEL);
ylabel(['z [' char(956) 'm]'], 'FontSize', FONT_LABEL);
title('(b) magnetic flux density on x-z plane (y=0)', 'FontSize', FONT_TITLE);

exportgraphics(fig_b, fullfile(fig_dir, 'fig2_4b.png'), 'Resolution', DPI);
fprintf('Saved fig2_4b.png\n');

fprintf('\nAll Fig 2.4 panels saved to %s\n', fig_dir);
