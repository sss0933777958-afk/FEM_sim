%% sweep_J_validity_hung.m — Sweep [J] fit cube sizes for hung hexapole
%   Mirrors magnetic_sim/ANSYS/main/analysis/fit/sweep_J_validity.m but for 6-pole hexapole.
%
%   Output: magnetic_sim/ANSYS/backup/hung/figures/coil1/fit_J_validity_boundary_hung.png

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
clear; clc;

%% ---- Paths ----
this_dir = fileparts(mfilename('fullpath'));
results_root = fullfile(this_dir, '..', '..', 'results');
fig_dir      = fullfile(this_dir, '..', '..', 'figures', 'coil1');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ---- Geometry constants (hung hexapole, l=250 — only set with all 6 coils) ----
% Use l=250 data because only coil1 has filleted_conv (l=500); all 6 have l250
R_norm     = 250e-6;             % 250 µm
R_norm_xy  = R_norm * sqrt(2/3); % 204 µm
R_norm_z   = R_norm / sqrt(3);   % 144 µm
POLE_TIP_R = 0.040e-3;           % pole tip cylinder radius

% Pole angles & layer (same as mt_constants but scaled to R_norm=250)
pole_angles = [0, 180, 120, 300, 60, 240];   % deg
pole_is_lower = [1, 0, 1, 0, 0, 1];          % P1-P6
tip_x = R_norm_xy * cosd(pole_angles);
tip_y = R_norm_xy * sind(pole_angles);

%% ---- Load 6 coils (l=250) ----
fprintf('Loading 6 hung coil datasets (filleted_l250_conv)...\n');
coil_raw = struct();
for k = 1:6
    dk = import_ansys_data(fullfile(results_root, sprintf('coil%d', k), 'filleted_l250_conv'), ...
                           'wp', sprintf('coil%d', k));
    coil_raw(k).x = dk.x; coil_raw(k).y = dk.y; coil_raw(k).z = dk.z;
    coil_raw(k).bx = dk.bx; coil_raw(k).by = dk.by; coil_raw(k).bz = dk.bz;
    coil_raw(k).bsum = dk.bsum;
    fprintf('  coil%d: %d nodes\n', k, length(dk.x));
end

%% ---- Sweep cube sizes ----
sweep_xy_um = [20, 30, 40, 50, 60, 75, 100, 125, 150, 175, 200, 225, 250];
n_sweep = length(sweep_xy_um);
err_pct  = nan(1, n_sweep);
r_fit_um = nan(1, n_sweep);
n_pts_arr = nan(1, n_sweep);

N_POLES = 6;
K_I = eye(N_POLES) - ones(N_POLES)/N_POLES;

% Initial guess: pole tip 3D positions
pos_init = zeros(3, N_POLES);
for k = 1:N_POLES
    pos_init(1, k) = tip_x(k);
    pos_init(2, k) = tip_y(k);
    if pole_is_lower(k)
        pos_init(3, k) = -R_norm_z;
    else
        pos_init(3, k) = +R_norm_z;
    end
end

fprintf('\n=== Sweeping fit cube sizes (hung hexapole, l=500) ===\n');
fprintf('  half_xy   N_pts   error    |r|_fit  Tip→Charge\n');
fprintf('  ------------------------------------------------\n');

for s = 1:n_sweep
    half_xy = sweep_xy_um(s) * 1e-6;
    half_z  = half_xy;          % cube (no z-slab limit for hexapole; poles span both layers)

    coil_data = struct();
    valid_run = true;
    for k = 1:N_POLES
        mask = abs(coil_raw(k).x) < half_xy & abs(coil_raw(k).y) < half_xy & ...
               abs(coil_raw(k).z) < half_z;
        n_in = sum(mask);
        N_max = 5000;
        if n_in > N_max
            idx_all = find(mask);
            rng(0);
            idx_sub = idx_all(randperm(n_in, N_max));
            mask = false(size(mask));
            mask(idx_sub) = true;
        end
        if sum(mask) < 30
            valid_run = false;
            break;
        end
        coil_data(k).px = coil_raw(k).x(mask);
        coil_data(k).py = coil_raw(k).y(mask);
        coil_data(k).pz = coil_raw(k).z(mask);
        coil_data(k).bmag  = sqrt(coil_raw(k).bx(mask).^2 + coil_raw(k).by(mask).^2 + coil_raw(k).bz(mask).^2);
        coil_data(k).b_fem = [coil_raw(k).bx(mask); coil_raw(k).by(mask); coil_raw(k).bz(mask)];
        coil_data(k).N     = sum(mask);
        Iv = zeros(N_POLES, 1); Iv(k) = 1;
        coil_data(k).KI_w  = K_I * Iv;
    end

    if ~valid_run
        fprintf('  %5.0f µm   too few nodes, skip\n', sweep_xy_um(s));
        continue;
    end

    % Fit
    opts = optimset('TolX', 1e-10, 'TolFun', 1e-22, ...
                    'MaxIter', 50000, 'MaxFunEvals', 50000, 'Display', 'off');
    cost_fn = @(x) cjoint(x, coil_data, 1e-7, N_POLES);
    [xo, ~] = fminsearch(cost_fn, pos_init(:), opts);
    pos = reshape(xo, 3, N_POLES);
    [~, Ck] = cjoint(xo, coil_data, 1e-7, N_POLES);

    epc = zeros(1, N_POLES);
    for k = 1:N_POLES
        Nk = coil_data(k).N;
        pw = [coil_data(k).px, coil_data(k).py, coil_data(k).pz];
        bu = evcf(pos, coil_data(k).KI_w, 1e-7, pw, N_POLES);
        bm = Ck(k) * bu;
        res = bm - coil_data(k).b_fem;
        em = sqrt(res(1:Nk).^2 + res(Nk+1:2*Nk).^2 + res(2*Nk+1:3*Nk).^2);
        epc(k) = mean(em ./ coil_data(k).bmag);
    end
    err_pct(s) = 100 * mean(epc);
    r_fit_um(s) = mean(vecnorm(pos)) * 1e6;
    n_pts_arr(s) = sum([coil_data.N]);

    fprintf('  %5.0f µm   %5d   %5.2f%%  %6.1f   %5.1f\n', ...
            sweep_xy_um(s), coil_data(1).N, err_pct(s), r_fit_um(s), r_fit_um(s)-R_norm*1e6);
end
fprintf('\n');

%% ---- Find validity thresholds ----
valid_mask = ~isnan(err_pct);
xs = sweep_xy_um(valid_mask);
ys = err_pct(valid_mask);
[xs, idx_s] = sort(xs);
ys = ys(idx_s);

threshold_5  = interp_threshold(xs, ys, 5);
threshold_10 = interp_threshold(xs, ys, 10);
threshold_20 = interp_threshold(xs, ys, 20);

fprintf('=== Validity boundaries (hung hexapole l=500 µm) ===\n');
if ~isnan(threshold_5),  fprintf('  5%%  threshold: cube_half ≈ %.0f µm  (%.0f%% of R_norm)\n',  threshold_5,  threshold_5/(R_norm*1e6)*100); end
if ~isnan(threshold_10), fprintf('  10%% threshold: cube_half ≈ %.0f µm  (%.0f%%)\n',  threshold_10, threshold_10/(R_norm*1e6)*100); end
if ~isnan(threshold_20), fprintf('  20%% threshold: cube_half ≈ %.0f µm  (%.0f%%)\n',  threshold_20, threshold_20/(R_norm*1e6)*100); end
fprintf('  Largest cube tested (%.0f µm): error ≈ %.1f%%\n', xs(end), ys(end));

%% ---- Plot 1×2 figure ----
FONT_NAME  = 'Helvetica';
FONT_LABEL = 12;
FONT_TITLE = 13;
LINE_W     = 1.8;
DPI        = 300;

fig = figure('Position', [60 60 1400 600], 'Color', 'w');
set(fig, 'DefaultAxesFontName', FONT_NAME, 'DefaultAxesFontSize', FONT_LABEL);

%% Left panel: error curve
ax1 = subplot(1, 2, 1);
hold on; box on; grid on;

plot(xs, ys, 'o-', 'Color', [0.85 0.10 0.10], 'MarkerFaceColor', [0.85 0.10 0.10], ...
     'MarkerSize', 7, 'LineWidth', LINE_W);

xl = [0, max(xs)*1.1];
yl = [0, max(50, max(ys)*1.1)];

plot(xl, [5 5],   '--', 'Color', [0.10 0.65 0.10], 'LineWidth', 1.2);
plot(xl, [10 10], '--', 'Color', [0.95 0.55 0.10], 'LineWidth', 1.2);
plot(xl, [20 20], '--', 'Color', [0.85 0.10 0.10], 'LineWidth', 1.2);

if ~isnan(threshold_5)
    plot([threshold_5 threshold_5], yl, ':', 'Color', [0.10 0.65 0.10], 'LineWidth', 1.0);
end
if ~isnan(threshold_10)
    plot([threshold_10 threshold_10], yl, ':', 'Color', [0.95 0.55 0.10], 'LineWidth', 1.0);
end

text(xl(2)*0.97, 5,  '5% (valid)',         'FontSize', 9, 'Color', [0.10 0.65 0.10], 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');
text(xl(2)*0.97, 10, '10% (transitional)', 'FontSize', 9, 'Color', [0.95 0.55 0.10], 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');
text(xl(2)*0.97, 20, '20% (invalid)',      'FontSize', 9, 'Color', [0.85 0.10 0.10], 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');

if ~isnan(threshold_5)
    text(threshold_5, yl(2)*0.92, sprintf('valid boundary\n%.0f µm', threshold_5), ...
        'FontSize', 9, 'Color', [0.10 0.65 0.10], 'HorizontalAlignment', 'center', ...
        'BackgroundColor', 'w', 'EdgeColor', [0.10 0.65 0.10], 'Margin', 2);
end

xlim(xl); ylim(yl);
xlabel('Fit cube xy half-edge  [µm]', 'FontSize', FONT_LABEL);
ylabel('Mean fit error  [%]', 'FontSize', FONT_LABEL);
title('(a) [J] 6-charge model error vs fit-region size  (hung hexapole, l=250 µm)', ...
      'FontSize', FONT_TITLE);

%% Right panel: top-down geometry + valid region
ax2 = subplot(1, 2, 2);
hold on; box on; grid on;

% Pole circles (POLE_TIP_R radius), top-down projection
theta_c = linspace(0, 2*pi, 100);
for k = 1:N_POLES
    cx = tip_x(k) * 1e6;
    cy = tip_y(k) * 1e6;
    if pole_is_lower(k)
        face_col = [0.85 0.85 1.00];   % light blue (Lower)
        layer_lab = 'L';
    else
        face_col = [1.00 0.90 0.85];   % light orange (Upper)
        layer_lab = 'U';
    end
    fill(cx + POLE_TIP_R*1e6*cos(theta_c), cy + POLE_TIP_R*1e6*sin(theta_c), face_col, ...
         'EdgeColor', 'k', 'LineWidth', 1.0);
    text(cx, cy, sprintf('P%d (%s)', k, layer_lab), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

% Validity region circles
if ~isnan(threshold_5)
    th = linspace(0, 2*pi, 200);
    fill(threshold_5*cos(th), threshold_5*sin(th), [0.6 0.95 0.6], ...
         'FaceAlpha', 0.4, 'EdgeColor', [0.10 0.65 0.10], 'LineWidth', 1.5);
end
if ~isnan(threshold_10)
    th = linspace(0, 2*pi, 200);
    plot(threshold_10*cos(th), threshold_10*sin(th), '-', 'Color', [0.95 0.55 0.10], 'LineWidth', 1.5);
end

% Physical workspace boundary (R_norm_xy in xy, where poles project to)
plot(R_norm_xy*1e6*cos(theta_c), R_norm_xy*1e6*sin(theta_c), '-.', 'Color', [0.85 0.10 0.10], 'LineWidth', 1.2);
% R_norm 3D for reference
plot(R_norm*1e6*cos(theta_c), R_norm*1e6*sin(theta_c), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);

% WP center
plot(0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);

% Legend annotations (bottom)
text_y = -400;
text(-380, text_y,        sprintf('• Green fill: valid region (err<5%%, r<%.0f µm)', threshold_5), ...
     'FontSize', 9, 'Color', [0.10 0.65 0.10]);
text(-380, text_y - 40,   sprintf('• Orange ring: transitional (err<10%%, r<%.0f µm)', threshold_10), ...
     'FontSize', 9, 'Color', [0.95 0.55 0.10]);
text(-380, text_y - 80,   sprintf('• Red dash-dot: physical xy boundary (R_{norm,xy} = %.0f µm)', R_norm_xy*1e6), ...
     'FontSize', 9, 'Color', [0.85 0.10 0.10]);
text(-380, text_y - 120,  sprintf('• Gray dash: 3D R_{norm} = %.0f µm', R_norm*1e6), ...
     'FontSize', 9, 'Color', [0.6 0.6 0.6]);

axis equal;
xlim([-400, 400]); ylim([-550, 400]);
xlabel('x  [µm]', 'FontSize', FONT_LABEL);
ylabel('y  [µm]', 'FontSize', FONT_LABEL);
title('(b) Spatial validity in xy plane (z = 0)', 'FontSize', FONT_TITLE);

%% Save
out_path = fullfile(fig_dir, 'fit_J_validity_boundary_hung.png');
exportgraphics(fig, out_path, 'Resolution', DPI);
fprintf('\nSaved: %s\n', out_path);


%% ============== Helpers ==============
function thr = interp_threshold(xs, ys, target)
    thr = NaN;
    for i = 1:length(xs)-1
        if (ys(i) <= target && ys(i+1) > target) || (ys(i) > target && ys(i+1) <= target)
            thr = xs(i) + (target - ys(i)) / (ys(i+1) - ys(i)) * (xs(i+1) - xs(i));
            return;
        end
    end
end

function [c, Ck] = cjoint(x, cd, km, N_POLES)
    p = reshape(x, 3, N_POLES);
    c = 0; Ck = zeros(N_POLES, 1);
    for k = 1:N_POLES
        pw = [cd(k).px, cd(k).py, cd(k).pz];
        bu = evcf(p, cd(k).KI_w, km, pw, N_POLES);
        Ck(k) = (bu' * cd(k).b_fem) / (bu' * bu);
        r = Ck(k)*bu - cd(k).b_fem;
        c = c + sum(r.^2);
    end
end

function b = evcf(p, w, km, pw, N_POLES)
    N = size(pw, 1);
    bx = zeros(N, 1); by = bx; bz = bx;
    for i = 1:N_POLES
        dx = pw(:,1) - p(1,i);
        dy = pw(:,2) - p(2,i);
        dz = pw(:,3) - p(3,i);
        r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
        bx = bx + (-w(i)) * dx ./ r3;
        by = by + (-w(i)) * dy ./ r3;
        bz = bz + (-w(i)) * dz ./ r3;
    end
    b = km * [bx; by; bz];
end
