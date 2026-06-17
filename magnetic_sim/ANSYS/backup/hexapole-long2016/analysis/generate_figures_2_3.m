%% generate_figures_2_3.m — Fig 2.3(a,b,c): B-field vector plots
% Corresponds to Long 2016 dissertation Fig. 2.3
%   (a) Top view (XY): dense colored quiver (ANSYS PLVECT style)
%   (b) WP center: 100 um cube 3D quiver (actuation coordinate)
%   (c) P1 pole tip: flux convergence with cone outline
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
LINE_THIN   = 0.8;
DPI         = 300;

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

% P1 pole tip in WP frame
tip_x_wp = c.R_norm_xy;    % ~408 um
tip_z_wp = -c.R_norm_z;    % ~-289 um

%% ======== Fig 2.3(a): XY top view — reference .fig style ========
% Style: AutoScale='off', U/V = B in Tesla, linewidth varies by |B|
fig_a = figure('Name', 'Fig 2.3(a)', 'Position', [50 50 850 750], 'Color', 'w');
set(fig_a, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

% Spatial binning: divide XY domain into grid cells, pick one node per cell
% This distributes arrows evenly regardless of mesh density
grid_cells = 80;   % 80x80 grid over ±80 mm => 2 mm/cell
x_mm = d.x * 1e3;
y_mm = d.y * 1e3;
x_edges = linspace(-80, 80, grid_cells+1);
y_edges = linspace(-80, 80, grid_cells+1);

% Assign each node to a grid cell using discretize (vectorized, fast)
ix_bin = discretize(x_mm, x_edges);
iy_bin = discretize(y_mm, y_edges);
valid_bin = ~isnan(ix_bin) & ~isnan(iy_bin);
cell_id = zeros(length(x_mm), 1);
cell_id(valid_bin) = (ix_bin(valid_bin)-1)*grid_cells + iy_bin(valid_bin);

% For each occupied cell, pick node with highest |B|
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

% Filter: keep only nodes with meaningful |B|
valid = d.bsum(idx) > 1e-4;   % > 0.1 mT threshold (lower for outer arrows)
idx = idx(valid);
x_q = d.x(idx)*1e3;   % mm
y_q = d.y(idx)*1e3;
bmag_q = d.bsum(idx);

% Power-law arrow scaling: B^0.25 compresses 1000:1 range to ~6:1
% Makes weak-field arrows visible while strong ones still stand out
arrow_max = 4.0;   % max arrow length in mm (at |B|=1T)
bxy = sqrt(d.bx(idx).^2 + d.by(idx).^2);   % in-plane magnitude
bxy(bxy == 0) = 1e-10;
scale_factor = arrow_max * bmag_q.^0.25 ./ bxy;  % direction from bx,by; length from |B|^0.25
bx_q = d.bx(idx) .* scale_factor;
by_q = d.by(idx) .* scale_factor;
fprintf('Fig 2.3(a): %d arrows after spatial binning + threshold\n', length(idx));

% Color bins with varying linewidth (matching reference style)
n_bins = 28;
bmax = max(bmag_q);
edges = linspace(0, bmax, n_bins + 1);
cmap_a = turbo(n_bins);
lw_range = [0.40, 1.75];   % thin for weak, thick for strong

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

% Yoke outlines (red dashed, matching reference)
theta_c = linspace(0, 2*pi, 200);
plot(c.YOKE_IN_R*1e3*cos(theta_c), c.YOKE_IN_R*1e3*sin(theta_c), ...
	'r--', 'LineWidth', LINE_THIN);
plot(c.YOKE_OUT_R*1e3*cos(theta_c), c.YOKE_OUT_R*1e3*sin(theta_c), ...
	'r--', 'LineWidth', LINE_THIN);

% Coil protrusion circles + pole labels
for k = 1:6
	ang = c.pole_angles(k) * pi/180;
	cx = c.YOKE_MID_R * cos(ang) * 1e3;
	cy = c.YOKE_MID_R * sin(ang) * 1e3;
	plot(cx + c.PROT_R*1e3*cos(theta_c), cy + c.PROT_R*1e3*sin(theta_c), ...
		'k-', 'LineWidth', LINE_MAIN);
	text(cx*1.15, cy*1.15, c.pole_labels{k}, ...
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

axis equal; grid on;
set(gca, 'Color', 'w');
xlim([-80 80]); ylim([-80 80]);
xlabel('x [mm]', 'FontSize', FONT_LABEL);
ylabel('y [mm]', 'FontSize', FONT_LABEL);
title('(a) Top view: B-field vector distribution (unit: Tesla)', ...
	'FontSize', FONT_TITLE);

exportgraphics(fig_a, fullfile(fig_dir, 'fig2_3a.png'), 'Resolution', DPI);
fprintf('Saved fig2_3a.png\n');

%% ======== Fig 2.3(b): 3D quiver at WP center (100 um cube) ========
fig_b = figure('Name', 'Fig 2.3(b)', 'Position', [100 100 750 700], 'Color', 'w');
set(fig_b, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

half_cube = 50e-6;   % +/-50 um
mask_b = abs(wp.x) < half_cube & ...
	abs(wp.y) < half_cube & ...
	abs(wp.z - c.SPH_OFST) < half_cube;

fprintf('Fig 2.3(b): %d nodes in 100 um cube\n', sum(mask_b));
if sum(mask_b) < 5
	half_cube = 100e-6;
	mask_b = abs(wp.x) < half_cube & ...
		abs(wp.y) < half_cube & ...
		abs(wp.z - c.SPH_OFST) < half_cube;
	fprintf('  Expanded to 200 um cube: %d nodes\n', sum(mask_b));
end

% Subsample to ~200 arrows for clarity
n_b = sum(mask_b);
sub_b = max(1, floor(n_b / 200));
idx_b = find(mask_b);
idx_b = idx_b(1:sub_b:end);

qx = wp.x(idx_b)*1e6;          % um
qy = wp.y(idx_b)*1e6;
qz = (wp.z(idx_b) - c.SPH_OFST)*1e6;
bx_b = wp.bx(idx_b);
by_b = wp.by(idx_b);
bz_b = wp.bz(idx_b);

% Small arrows: normalize to unit direction, scale to ~4 um length
bmag_raw = sqrt(bx_b.^2 + by_b.^2 + bz_b.^2);
bmag_raw(bmag_raw == 0) = 1;
arrow_len = 4;  % um
ux_b = bx_b ./ bmag_raw * arrow_len;
uy_b = by_b ./ bmag_raw * arrow_len;
uz_b = bz_b ./ bmag_raw * arrow_len;

hold on;
quiver3(qx, qy, qz, ux_b, uy_b, uz_b, 0, ...
	'b', 'LineWidth', 0.5, 'MaxHeadSize', 0.2);
hold off;

axis equal; grid on; box on;
set(gca, 'Color', 'w');
xlabel(['x(' char(956) 'm)'], 'FontSize', FONT_LABEL);
ylabel(['y(' char(956) 'm)'], 'FontSize', FONT_LABEL);
zlabel(['z(' char(956) 'm)'], 'FontSize', FONT_LABEL);
title('(b) Magnetic flux density vectors near workspace center', ...
	'FontSize', FONT_TITLE);
view([-37.5, 30]);

lim_b = 60;
xlim([-lim_b lim_b]); ylim([-lim_b lim_b]); zlim([-lim_b lim_b]);

exportgraphics(fig_b, fullfile(fig_dir, 'fig2_3b.png'), 'Resolution', DPI);
fprintf('Saved fig2_3b.png\n');

%% ======== Fig 2.3(c): Flux convergence near P1 tip (XY plane, WP frame) ========
fig_c = figure('Name', 'Fig 2.3(c)', 'Position', [150 150 800 700], 'Color', 'w');
set(fig_c, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

% XY slice at z_wp ~ tip_z (P1 tip level), WP-centered coordinates
% P1 tip in WP frame: (R_norm_xy, 0, -R_norm_z) = (~408, 0, ~-289) um
z_slice_wp = tip_z_wp;           % z_wp of P1 tip
z_tol  = 0.3e-3;                 % half-thickness [m]
xrange_wp = [-100e-6, 1200e-6];  % WP frame [m], focus on P1 region
yrange_wp = [-500e-6, 500e-6];

mask_c = abs((d.z - c.SPH_OFST) - z_slice_wp) < z_tol & ...
	d.x > xrange_wp(1) & d.x < xrange_wp(2) & ...
	d.y > yrange_wp(1) & d.y < yrange_wp(2);

fprintf('Fig 2.3(c): %d nodes in XY slice at z_wp = %.0f um\n', ...
	sum(mask_c), z_slice_wp*1e6);
if sum(mask_c) < 200
	z_tol  = 0.6e-3;
	xrange_wp = [-200e-6, 1500e-6];
	yrange_wp = [-700e-6, 700e-6];
	mask_c = abs((d.z - c.SPH_OFST) - z_slice_wp) < z_tol & ...
		d.x > xrange_wp(1) & d.x < xrange_wp(2) & ...
		d.y > yrange_wp(1) & d.y < yrange_wp(2);
	fprintf('  Expanded: %d nodes\n', sum(mask_c));
end

% WP-centered coordinates
x_c_all = d.x(mask_c)*1e6;          % um
y_c_all = d.y(mask_c)*1e6;          % um
bx_c_all  = d.bx(mask_c);
by_c_all  = d.by(mask_c);
bsum_c_all = d.bsum(mask_c);

% Spatial binning for even arrow distribution
gc_x = 50;  gc_y = 40;
xedge_c = linspace(min(x_c_all), max(x_c_all), gc_x+1);
yedge_c = linspace(min(y_c_all), max(y_c_all), gc_y+1);
ix_c = discretize(x_c_all, xedge_c);
iy_c = discretize(y_c_all, yedge_c);
valid_c = ~isnan(ix_c) & ~isnan(iy_c);
cell_c = zeros(length(x_c_all), 1);
cell_c(valid_c) = (ix_c(valid_c)-1)*gc_y + iy_c(valid_c);
qi = [];
ucells = unique(cell_c(cell_c > 0));
for kk = 1:length(ucells)
	cn = find(cell_c == ucells(kk));
	[~, best] = max(bsum_c_all(cn));
	qi(end+1) = cn(best); %#ok<SAGROW>
end
fprintf('  Spatial binning: %d arrows\n', length(qi));

% Two-tone: red = iron (>50 mT), blue = air (>0.5 mT)
bsum_qi = bsum_c_all(qi);
iron = bsum_qi > 0.05;
air  = bsum_qi > 0.0005 & ~iron;

% Arrow scale: normalize to unit direction, fixed length
bxy_qi = sqrt(bx_c_all(qi).^2 + by_c_all(qi).^2);
bxy_qi(bxy_qi == 0) = 1;
arrow_um = 20;   % arrow length in um
ux_c = bx_c_all(qi) ./ bxy_qi * arrow_um;
uy_c = by_c_all(qi) ./ bxy_qi * arrow_um;

hold on;
% Draw P1 cone outline in XY view at z = tip_z
% Cone axis: +x from tip. Cross-section at this z is y = ±r(s)
s_cone = linspace(0, 0.8e-3, 100);
r_cone = c.POLE_TIP_R + s_cone * (c.POLE_R - c.POLE_TIP_R) / c.POLE_CONE_LEN;
x_cone = (tip_x_wp + s_cone) * 1e6;   % um
y_upper = r_cone * 1e6;
y_lower = -r_cone * 1e6;
fill([x_cone, fliplr(x_cone)], [y_upper, fliplr(y_lower)], ...
	[0.85 0.85 0.85], 'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 1.0, ...
	'FaceAlpha', 0.4);

% Air arrows (blue)
if any(air)
	quiver(x_c_all(qi(air)), y_c_all(qi(air)), ux_c(air), uy_c(air), 0, ...
		'Color', [0.2 0.3 0.85], 'LineWidth', 0.4, 'MaxHeadSize', 0.2);
end
% Iron arrows (red, on top)
if any(iron)
	quiver(x_c_all(qi(iron)), y_c_all(qi(iron)), ux_c(iron), uy_c(iron), 0, ...
		'Color', [0.85 0.15 0.15], 'LineWidth', 0.6, 'MaxHeadSize', 0.2);
end

% Mark tip position
plot(tip_x_wp*1e6, 0, 'ko', 'MarkerSize', 5, ...
	'MarkerFaceColor', 'k', 'LineWidth', 1);
hold off;

axis equal; grid on;
set(gca, 'Color', 'w');
xlabel(['x [' char(956) 'm]'], 'FontSize', FONT_LABEL);
ylabel(['y [' char(956) 'm]'], 'FontSize', FONT_LABEL);
title(sprintf('(c) B-field vectors near P1 tip (z_{wp} = %.0f %cm)', ...
	z_slice_wp*1e6, char(956)), 'FontSize', FONT_TITLE);

exportgraphics(fig_c, fullfile(fig_dir, 'fig2_3c.png'), 'Resolution', DPI);
fprintf('Saved fig2_3c.png\n');

fprintf('\nAll Fig 2.3 panels saved to %s\n', fig_dir);
