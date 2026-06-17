% plot_upper_pole_sensor_placement_4p472.m
% ----------------------------------------------------------------------
% Variant: sensor moved 100 um closer to tip (4.572 -> 4.472 mm along cone slant).
% All other geometry/style identical to plot_upper_pole_sensor_placement.m.
% Output: upper_pole_sensor_placement_4p472.png
% ----------------------------------------------------------------------

clear; close all;

% ---- Geometry (CAD intent, mm) ---------------------------------------
POLE_R         = 3.175;
POLE_CONE_LEN  = 15.876;
cyl_len        = 14.124;
half_angle_deg = atan2d(POLE_R, POLE_CONE_LEN);   % 11.31 deg
beta           = deg2rad(half_angle_deg);

tilt_deg = 90 - 54.7356;                          % 35.2644 deg
tilt     = deg2rad(tilt_deg);

ux =  cos(tilt);
uy =  sin(tilt);
px = -sin(tilt);
py =  cos(tilt);

sx_slant = cos(tilt + beta);
sy_slant = sin(tilt + beta);
nx = -sin(tilt + beta);
ny =  cos(tilt + beta);

% ---- Pole outline vertices (full natural cone + cylinder) ------------
A   = [0; 0];
J   = POLE_CONE_LEN     * [ux; uy];
T   = (POLE_CONE_LEN + cyl_len) * [ux; uy];
Ju  = J + POLE_R * [px; py];
Jl  = J - POLE_R * [px; py];
Tu  = T + POLE_R * [px; py];
Tl  = T - POLE_R * [px; py];
pole_x = [A(1), Ju(1), Tu(1), Tl(1), Jl(1), A(1)];
pole_y = [A(2), Ju(2), Tu(2), Tl(2), Jl(2), A(2)];

% ---- Sensor position (CHANGED: 4.572 -> 4.472 mm) -------------------
s_slant  = 4.472;   % mm ALONG the cone slant surface from apex
s_offset = 0.41;    % mm PERPENDICULAR to the cone surface, outward

foot = s_slant * [sx_slant; sy_slant];
sx_sensor = foot(1) + s_offset * nx;
sy_sensor = foot(2) + s_offset * ny;

% ---- Figure ----------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 900 800]);
ax  = axes(fig); hold(ax, 'on');

fill(ax, pole_x, pole_y, [0.78 0.80 0.84], ...
     'EdgeColor', [0.20 0.20 0.25], 'LineWidth', 1.8);

plot(ax, [Ju(1) Jl(1)], [Ju(2) Jl(2)], ...
     ':', 'Color', [0.55 0.55 0.60], 'LineWidth', 0.8);

axis_end = T + 1.5 * [ux; uy];
plot(ax, [A(1) axis_end(1)], [A(2) axis_end(2)], ...
     '-.', 'Color', [0.78 0.78 0.83], 'LineWidth', 0.4);

plot(ax, sx_sensor, sy_sensor, 'o', 'MarkerSize', 9, ...
     'MarkerFaceColor', [0.85 0.10 0.10], 'MarkerEdgeColor', [0.55 0.05 0.05]);

% ---- Dimension 1: 4.472 mm along the cone slant -----------------------
dim_off = 2.6;
D1a = A    + dim_off * [nx; ny];
D1b = foot + dim_off * [nx; ny];
plot(ax, [D1a(1) D1b(1)], [D1a(2) D1b(2)], '-', ...
     'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
tick_len = 0.18;
plot(ax, [D1a(1) - tick_len*sx_slant, D1a(1) + tick_len*sx_slant], ...
        [D1a(2) - tick_len*sy_slant, D1a(2) + tick_len*sy_slant], ...
     '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
plot(ax, [D1b(1) - tick_len*sx_slant, D1b(1) + tick_len*sx_slant], ...
        [D1b(2) - tick_len*sy_slant, D1b(2) + tick_len*sy_slant], ...
     '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
plot(ax, [A(1)    D1a(1)], [A(2)    D1a(2)], '--', ...
     'Color', [0.85 0.55 0.55], 'LineWidth', 0.6);
plot(ax, [foot(1) D1b(1)], [foot(2) D1b(2)], '--', ...
     'Color', [0.85 0.55 0.55], 'LineWidth', 0.6);
mid_D1 = (D1a + D1b)/2 + 0.35 * [nx; ny];
text(ax, mid_D1(1), mid_D1(2), '4.472 mm', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Color', [0.65 0.05 0.05], 'FontSize', 12, 'FontWeight', 'bold', ...
     'Rotation', tilt_deg + half_angle_deg);

% ---- Dimension 2: 0.41 mm perpendicular (cone slant -> sensor) -------
tick2 = 0.18;
plot(ax, [foot(1) - tick2*sx_slant, foot(1) + tick2*sx_slant], ...
        [foot(2) - tick2*sy_slant, foot(2) + tick2*sy_slant], ...
     '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
plot(ax, [sx_sensor - tick2*sx_slant, sx_sensor + tick2*sx_slant], ...
        [sy_sensor - tick2*sy_slant, sy_sensor + tick2*sy_slant], ...
     '-', 'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
plot(ax, [foot(1) sx_sensor], [foot(2) sy_sensor], '-', ...
     'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
lbl2 = [sx_sensor; sy_sensor] + 0.6 * [nx; ny];
text(ax, lbl2(1), lbl2(2), '0.41 mm', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Color', [0.65 0.05 0.05], 'FontSize', 12, 'FontWeight', 'bold', ...
     'Rotation', tilt_deg + half_angle_deg);

% ---- Annotations ------------------------------------------------------
plot(ax, A(1), A(2), 'k.', 'MarkerSize', 14);
text(ax, A(1) - 0.4*ux - 0.7*px, A(2) - 0.4*uy - 0.7*py, 'Apex (0, 0)', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
     'FontSize', 10, 'Color', [0.2 0.2 0.2]);

text(ax, sx_sensor + 0.5*nx + 0.8*sx_slant, ...
        sy_sensor + 0.5*ny + 0.8*sy_slant, ...
     'Hall sensor (in air)', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
     'FontSize', 10, 'Color', [0.65 0.05 0.05], ...
     'Rotation', tilt_deg + half_angle_deg);

mid_cone = 0.55 * (A + Ju) - 0.4 * [nx; ny];
text(ax, mid_cone(1), mid_cone(2), ...
     sprintf('Cone face (natural, %.2f\\circ half-angle)', half_angle_deg), ...
     'HorizontalAlignment', 'center', 'FontSize', 9, ...
     'Color', [0.20 0.20 0.30], 'Rotation', tilt_deg + half_angle_deg);

mid_cyl = 0.5 * (Ju + Tu) - 0.5 * [px; py];
text(ax, mid_cyl(1), mid_cyl(2), 'Cylinder body (R = 3.175 mm)', ...
     'HorizontalAlignment', 'center', 'FontSize', 9, ...
     'Color', [0.20 0.20 0.30], 'Rotation', tilt_deg);

text(ax, axis_end(1) + 0.3*ux, axis_end(2) + 0.3*uy, ...
     sprintf('Pole axis  (tilt %.2f\\circ above horizontal)', tilt_deg), ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
     'FontSize', 9, 'Color', [0.45 0.45 0.55], ...
     'Rotation', tilt_deg);

text(ax, 12, -3.0, ...
     ['Hexapole assembly: upper pole at polar angle ' ...
      '\alpha = 54.74\circ from +z (apex closest to WP, body up to yoke).'], ...
     'HorizontalAlignment', 'center', 'FontSize', 8.5, ...
     'Color', [0.35 0.35 0.45], 'FontAngle', 'italic');
text(ax, 12, -4.0, ...
     ['Pole is NOT milled (full natural cone). ' ...
      '4.472 mm is measured along the cone slant; 0.41 mm is perpendicular to that surface.'], ...
     'HorizontalAlignment', 'center', 'FontSize', 8.5, ...
     'Color', [0.35 0.35 0.45], 'FontAngle', 'italic');

% --- Axes formatting ---------------------------------------------------
axis(ax, 'equal');
xlim(ax, [-4, 30]);
ylim(ax, [-5.5, 24]);
xlabel(ax, 'x [mm]');
ylabel(ax, 'y [mm]');
title(ax, ['Upper Pole Sensor Placement — Side View ', ...
           '(full cone, hexapole tilt \alpha = 54.74\circ, mm) — sensor 4.472 mm']);
set(ax, 'FontSize', 10, 'Box', 'on', 'Layer', 'top');
grid(ax, 'off');

% --- Save --------------------------------------------------------------
out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
out_path = fullfile(out_dir, 'upper_pole_sensor_placement_4p472.png');
exportgraphics(fig, out_path, 'Resolution', 200);
fprintf('Saved: %s\n', out_path);
fprintf('Sensor center (mm): (%.4f, %.4f)\n', sx_sensor, sy_sensor);
fprintf('Foot on cone slant (mm): (%.4f, %.4f)\n', foot(1), foot(2));
