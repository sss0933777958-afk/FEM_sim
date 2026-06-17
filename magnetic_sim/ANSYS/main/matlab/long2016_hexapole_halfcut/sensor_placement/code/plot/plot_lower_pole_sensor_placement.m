% plot_lower_pole_sensor_placement.m
% ----------------------------------------------------------------------
% Side view of Long Fei lower magnetic pole as MILLED HALF-CONE +
% HALF-CYLINDER (the upper half along the pole axis has been ground
% away, leaving a D-shaped cross-section).  The flat top is the
% sensing surface; the Hall sensor sits 0.41 mm above it in air.
%
% Pole geometry (CAD intent, mm):
%   - apex at (0, 0); pole axis along +x; y = 0 = milled sensing surface
%   - upper boundary: horizontal line y = 0 from apex to cylinder tail
%   - lower boundary: cone slant from (0,0) to (15.876, -3.175),
%     then horizontal at y = -3.175 to (cyl_end, -3.175)
%   - POLE_R = 3.175 mm (1/8 inch), POLE_CONE_LEN = 15.876 mm,
%     half-angle = atan(POLE_R / POLE_CONE_LEN) = 11.31 deg
%
% Sensor position (design intent):
%   - along pole axis from apex: 4.572 mm
%   - perpendicular to sensing surface (vertical here): 0.41 mm
%   - sensor center = (4.572, 0.41)
%
% Output:  magnetic_sim/ANSYS/main/figures/long2016_hexapole_halfcut/lower_pole_sensor_placement.png
% ----------------------------------------------------------------------

clear; close all;

% ---- Geometry (CAD intent, mm) ---------------------------------------
POLE_R         = 3.175;
POLE_CONE_LEN  = 15.876;
cyl_end        = 30.0;
half_angle_deg = atan2d(POLE_R, POLE_CONE_LEN);   % 11.31 deg

% Sensor design intent
s_along  = 4.572;        % mm from apex along axis
s_offset = 0.41;         % mm perpendicular (upward) from sensing surface

% Sensor center (very simple now: just above y=0 by 0.41 mm)
sx = s_along;
sy = s_offset;

% ---- Pole outline polygon (milled half-cone + half-cylinder) --------
% Trace counter-clockwise starting at apex:
%   apex(0,0) -> bottom of cone (15.876, -R) -> bottom of cyl (cyl_end, -R)
%   -> tail top (cyl_end, 0) -> back to apex along sensing surface (y=0)
pole_x = [0, POLE_CONE_LEN, cyl_end, cyl_end, 0];
pole_y = [0, -POLE_R,       -POLE_R, 0,       0];

% ---- Figure ----------------------------------------------------------
fig = figure('Color', 'w', 'Position', [100 100 1100 420]);
ax  = axes(fig); hold(ax, 'on');

% Pole material (filled, light grey)
fill(ax, pole_x, pole_y, [0.78 0.80 0.84], ...
     'EdgeColor', [0.20 0.20 0.25], 'LineWidth', 1.8);

% Cone -> cylinder junction line (light dotted)
plot(ax, [POLE_CONE_LEN POLE_CONE_LEN], [-POLE_R 0], ...
     ':', 'Color', [0.55 0.55 0.60], 'LineWidth', 0.8);

% Pole axis reference (dash-dot, very faint, where the axis WOULD be
% in the un-milled cone — i.e. y = 0, same as sensing surface here)
plot(ax, [0, cyl_end], [0 0], '-.', 'Color', [0.75 0.75 0.80], 'LineWidth', 0.4);

% Sensor (red filled circle, marker visible at this zoom)
plot(ax, sx, sy, 'o', 'MarkerSize', 9, ...
     'MarkerFaceColor', [0.85 0.10 0.10], 'MarkerEdgeColor', [0.55 0.05 0.05]);

% --- Dimension 1: 4.572 mm horizontal (apex along sensing surface)
dim_y_top = 1.35;                                  % above the pole, in clear air
plot(ax, [0 s_along], [dim_y_top dim_y_top], '-', ...
     'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
plot(ax, [0 0],             dim_y_top + [-0.10 0.10], '-', ...
     'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
plot(ax, [s_along s_along], dim_y_top + [-0.10 0.10], '-', ...
     'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
text(ax, s_along/2, dim_y_top + 0.15, '4.572 mm', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'Color', [0.65 0.05 0.05], 'FontSize', 12, 'FontWeight', 'bold');

% --- Dimension 2: 0.41 mm vertical (sensing surface -> sensor center)
dim_x_right = s_along + 1.8;                       % to the right of sensor
plot(ax, [dim_x_right dim_x_right], [0 s_offset], '-', ...
     'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
plot(ax, dim_x_right + [-0.12 0.12], [0 0], '-', ...
     'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
plot(ax, dim_x_right + [-0.12 0.12], [s_offset s_offset], '-', ...
     'Color', [0.80 0.10 0.10], 'LineWidth', 1.4);
% Light dashed line from sensor across to the dimension stack
plot(ax, [sx dim_x_right], [sy sy], '--', ...
     'Color', [0.85 0.55 0.55], 'LineWidth', 0.6);
plot(ax, [sx dim_x_right], [0  0 ], '--', ...
     'Color', [0.85 0.55 0.55], 'LineWidth', 0.6);
text(ax, dim_x_right + 0.25, s_offset/2, '0.41 mm', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
     'Color', [0.65 0.05 0.05], 'FontSize', 12, 'FontWeight', 'bold');

% --- Annotations -------------------------------------------------------
% Apex
plot(ax, 0, 0, 'k.', 'MarkerSize', 14);
text(ax, -0.3, -0.35, 'Apex (0, 0)', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 10, 'Color', [0.2 0.2 0.2]);

% Sensor label
text(ax, sx + 0.35, sy + 0.05, 'Hall sensor (in air)', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
     'FontSize', 10, 'Color', [0.65 0.05 0.05]);

% Sensing surface label — placed in the empty cylinder-top stretch
text(ax, 22, 0.55, 'Sensing surface — milled flat (y = 0, through pole axis)', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'FontSize', 9.5, 'Color', [0.20 0.20 0.30]);
% small leader from label down to the surface
plot(ax, [22 22], [0.50 0.05], '-', 'Color', [0.45 0.45 0.50], 'LineWidth', 0.5);

% Lower cone label
text(ax, 7.5, -2.4, ...
     sprintf('Lower cone face (slope -0.2, %.2f\\circ)', half_angle_deg), ...
     'HorizontalAlignment', 'center', 'FontSize', 9, ...
     'Color', [0.25 0.25 0.30], 'Rotation', -half_angle_deg);

% Cylinder body label
text(ax, (POLE_CONE_LEN + cyl_end)/2, -1.55, ...
     'Half-cylinder body (R = 3.175 mm)', ...
     'HorizontalAlignment', 'center', 'FontSize', 9, ...
     'Color', [0.25 0.25 0.30]);

% Cone -> cylinder junction marker
text(ax, POLE_CONE_LEN, -POLE_R - 0.35, ...
     sprintf('x = %.3f mm\n(cone\\rightarrowcyl)', POLE_CONE_LEN), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'FontSize', 8, 'Color', [0.45 0.45 0.50]);

% --- Axes formatting ---------------------------------------------------
axis(ax, 'equal');
xlim(ax, [-2, cyl_end + 2]);
ylim(ax, [-POLE_R - 1.3, 2.0]);
xlabel(ax, 'x [mm]  (along pole axis from apex)');
ylabel(ax, 'y [mm]');
title(ax, 'Lower Pole Sensor Placement — Side View (milled half-cone, mm)');
set(ax, 'FontSize', 10, 'Box', 'on', 'Layer', 'top');
grid(ax, 'off');

% --- Save --------------------------------------------------------------
out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
out_path = fullfile(out_dir, 'lower_pole_sensor_placement.png');
exportgraphics(fig, out_path, 'Resolution', 200);
fprintf('Saved: %s\n', out_path);
