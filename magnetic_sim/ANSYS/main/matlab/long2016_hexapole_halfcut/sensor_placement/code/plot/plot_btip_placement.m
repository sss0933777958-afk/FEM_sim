%% plot_btip_placement.m — visualize B_surface + B_tip sensor placement on P1
%
%  Side view of P1 lower pole (milled half-cone, same as image 26):
%   - Pole apex at (0,0)
%   - Cone face going down-right (-11.31° slope) to x=15.876
%   - Half-cylinder body from x=15.876 to x=30
%   - Milled flat at y=0 (extends from apex along +x)
%
%  Mark BOTH sensors:
%   - B_surface @ pole-local (4.572, 0.41, 0) — 既有 (image 26)
%   - B_tip     @ pole-local (-0.25, 0, 0)   — new

clear; clc; close all;

fig = figure('Position',[60 60 1400 380],'Color','w');
ax = axes; hold on;

%% --- Draw lower pole (milled half-cone) side view ---
% Cone face: apex (0,0) → (15.876, -3.175)  (slope tan(-11.31°) = -0.2)
% Half-cylinder body: x ∈ [15.876, 30], y ∈ [-3.175, 0]
% Milled flat at y = 0 from x=0 to x=30 (through axis)

% Pole outline polygon
pole_x = [0, 15.876, 30, 30, 0];
pole_y = [0, -3.175, -3.175, 0, 0];
fill(pole_x, pole_y, [0.85 0.87 0.92], 'EdgeColor', [0.25 0.25 0.30], 'LineWidth', 1.5);

%% --- Apex marker ---
plot(0, 0, 'k.', 'MarkerSize', 18);
text(-0.5, -0.55, 'Apex (0, 0)', 'FontSize', 11, 'HorizontalAlignment','right');

%% --- B_surface sensor (existing, image 26 position) ---
xs = 4.572; ys = 0.41;
plot(xs, ys, 'o', 'MarkerSize', 14, 'MarkerFaceColor', [0.85 0.20 0.20], ...
                                    'MarkerEdgeColor', 'k', 'LineWidth', 1);
text(xs+0.7, ys+0.9, 'B_{surface}', 'FontSize', 12, 'Color', [0.7 0 0], 'FontWeight','bold');
text(xs+0.7, ys+0.45, '(4.572, 0.41) mm', 'FontSize', 9, 'Color', [0.7 0 0]);
% sensor normal arrow (up, +y)
quiver(xs, ys, 0, 0.7, 0, 'Color', [0.85 0.20 0.20], 'LineWidth', 2, 'MaxHeadSize', 0.7);
text(xs-2.5, ys+1.25, 'n = (0, +1, 0)', 'FontSize', 9, 'Color', [0.85 0.20 0.20]);
% dimension lines for B_surface
plot([0 xs], [ys+0.05 ys+0.05], '--', 'Color', [0.85 0.20 0.20], 'LineWidth', 0.8);
text(xs/2, ys+0.32, '4.572 mm', 'FontSize', 10, 'Color', [0.85 0.20 0.20], ...
     'HorizontalAlignment','center', 'FontWeight','bold');
plot([xs xs], [0 ys], '--', 'Color', [0.85 0.20 0.20], 'LineWidth', 0.8);
text(xs+0.1, ys/2, ' 0.41 mm', 'FontSize', 9, 'Color', [0.85 0.20 0.20]);

%% --- B_tip sensor (NEW, 250 µm in front of apex, normal → +x toward apex) ---
xt = -0.25; yt = 0;
plot(xt, yt, 'o', 'MarkerSize', 14, 'MarkerFaceColor', [0.20 0.50 0.85], ...
                                    'MarkerEdgeColor', 'k', 'LineWidth', 1);
text(xt-0.7, yt+1.4, 'B_{tip}', 'FontSize', 12, 'Color', [0 0.2 0.7], 'FontWeight','bold');
text(xt-0.7, yt+0.95, '(-0.25, 0) mm', 'FontSize', 9, 'Color', [0 0.2 0.7]);
% sensor normal arrow: face → apex, so points in +x (from sensor toward apex)
quiver(xt+0.04, yt, 0.18, 0, 0, 'Color', [0.20 0.50 0.85], 'LineWidth', 2.5, 'MaxHeadSize', 4);
text(xt+0.2, yt-0.6, 'n = (+1, 0, 0)  toward apex', 'FontSize', 9, 'Color', [0.20 0.50 0.85]);
% dimension: 0.25 mm from apex
plot([xt 0], [yt-1.0 yt-1.0], '-', 'Color', [0.20 0.50 0.85], 'LineWidth', 1.2);
plot([xt xt], [yt-1.1 yt-0.9], '-', 'Color', [0.20 0.50 0.85], 'LineWidth', 1.2);
plot([0 0], [yt-1.1 yt-0.9], '-', 'Color', [0.20 0.50 0.85], 'LineWidth', 1.2);
text(xt/2, yt-1.4, '0.25 mm', 'FontSize', 10, 'Color', [0.20 0.50 0.85], ...
     'HorizontalAlignment','center', 'FontWeight','bold');

%% --- Labels and styling ---
text(22, -1.5, 'Half-cylinder body (R = 3.175 mm)', 'FontSize', 10, ...
     'Color', [0.30 0.30 0.35], 'HorizontalAlignment','center');
text(8, -1.0, 'Lower cone face (slope -0.2, 11.31°)', 'FontSize', 10, ...
     'Color', [0.30 0.30 0.35], 'Rotation', -11.31);
text(20, 0.5, 'Sensing surface — milled flat (y = 0, through pole axis)', ...
     'FontSize', 10, 'Color', [0.30 0.30 0.35], 'HorizontalAlignment','center');
text(15.876, -3.6, 'x = 15.876 mm (cone\rightarrowcyl)', 'FontSize',8, ...
     'HorizontalAlignment','center', 'Color', [0.50 0.50 0.55]);
plot([15.876 15.876], [-3.175 0], ':', 'Color', [0.5 0.5 0.55], 'LineWidth', 0.6);

% WP direction indicator (left side, since pole apex faces WP)
text(-1.5, -2.5, '\leftarrow WP (toward -x)', 'FontSize', 11, 'Color', [0 0.5 0], ...
     'FontWeight','bold', 'HorizontalAlignment','left');

%% --- Axes ---
xlim([-2 31]); ylim([-4.5 2]);
xlabel('x [mm]  (along pole axis from apex)', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
title('P1 Lower Pole — B_{surface} (red) and B_{tip} (blue) sensor placement', ...
      'FontSize', 13, 'FontWeight','bold');
grid on; box on;
set(ax, 'GridAlpha', 0.3);
axis equal;
xlim([-2 31]); ylim([-4.5 2]);

%% --- Save ---
out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
out_path = fullfile(out_dir, 'btip_sensor_placement.png');
exportgraphics(fig, out_path, 'Resolution', 300);
fprintf('Saved: %s\n', out_path);
