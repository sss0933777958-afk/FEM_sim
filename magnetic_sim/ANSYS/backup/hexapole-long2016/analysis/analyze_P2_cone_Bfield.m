%% analyze_P2_cone_Bfield.m — P2 cone interior B-field analysis
%
%  Method A coordinate system:
%    n_hat  (fixed): cone surface outward normal, perpendicular to surface
%    tau(d) (variable): always points from current position toward tip
%    n_hat perp tau only at d = d_max (surface); angle ~78.8-90 deg
%
%  Part 1 (Sections 0-4): Iron interior analysis using PATH data
%  Part 2 (Sections 5-7): Air-side comparison using full-mesh data
%
%  Run from magnetic_sim/hexapole-long2016/analysis/ directory.

clear; clc; close all;

%% 0. Constants, style, geometry
cnst = mt_constants();

% Figure style (match generate_figures_2_3.m)
FONT_LABEL = 12;
FONT_TITLE = 13;
LINE_MAIN  = 1.5;
LINE_THIN  = 1.2;
DPI        = 300;
fig_dir    = fullfile('..', 'figures');

% P2 geometry
pole_idx = 2;
tip_wp   = [cnst.pole_tip_x(pole_idx); cnst.pole_tip_y(pole_idx); cnst.pole_tip_z_wp(pole_idx)];
tip_ans  = tip_wp + [0; 0; cnst.SPH_OFST];  % tip in ANSYS coords
ax       = cnst.pole_axis(:, pole_idx);       % tip -> base
beta     = atan((cnst.POLE_R - cnst.POLE_TIP_R) / cnst.POLE_CONE_LEN);

% Transverse directions for the two azimuthal path families
t_hat = cross(ax, [0;1;0]); t_hat = t_hat/norm(t_hat);  % phi=0
y_hat = cross(ax, t_hat);   y_hat = y_hat/norm(y_hat);   % phi=90

% n_hat for each azimuthal direction (fixed per path family)
n_hat_t = -sin(beta)*ax + cos(beta)*t_hat;  n_hat_t = n_hat_t/norm(n_hat_t);
n_hat_y = -sin(beta)*ax + cos(beta)*y_hat;  n_hat_y = n_hat_y/norm(n_hat_y);

s_vals  = [2 3 4 5 6 7 8 10];  % [mm]
n_s     = length(s_vals);
results_dir = fullfile('..', 'results', 'coil5');

fprintf('P2 cone B-field analysis (Method A)\n');
fprintf('  beta = %.2f deg\n', rad2deg(beta));
fprintf('  tip (ANSYS) = (%.4e, %.4e, %.4e) m\n', tip_ans);
fprintf('  ax = (%.4f, %.4f, %.4f)\n', ax);

%% 1. PPATH endpoint coordinates (from post_path_P2_cone.txt)
% p1 = axis point, p2 = surface point, in ANSYS global coords [m]
PP = struct();
PP(1).s=2;  PP(1).t1=[-2.0829721382e-03, 0, -1.1179366119e-02]; PP(1).t2=[-2.2731968536e-03, 0, -1.1579499704e-02]; PP(1).y2=[-2.0141017236e-03, -4.3466666667e-04, -1.1230494224e-02];
PP(2).s=3;  PP(2).t1=[-2.9171651779e-03, 0, -1.0560076833e-02]; PP(2).t2=[-3.1937495801e-03, 0, -1.1141866156e-02]; PP(2).y2=[-2.8170284401e-03, -6.3200000000e-04, -1.0634416471e-02];
PP(3).s=4;  PP(3).t1=[-3.7513582177e-03, 0, -9.9407875471e-03]; PP(3).t2=[-4.1143023067e-03, 0, -1.0704232608e-02]; PP(3).y2=[-3.6199551567e-03, -8.2933333333e-04, -1.0038338718e-02];
PP(4).s=5;  PP(4).t1=[-4.5855512574e-03, 0, -9.3214982613e-03]; PP(4).t2=[-5.0348550332e-03, 0, -1.0266599061e-02]; PP(4).y2=[-4.4228818732e-03, -1.0266666667e-03, -9.4422609649e-03];
PP(5).s=6;  PP(5).t1=[-5.4197442972e-03, 0, -8.7022089755e-03]; PP(5).t2=[-5.9554077597e-03, 0, -9.8289655128e-03]; PP(5).y2=[-5.2258085898e-03, -1.2240000000e-03, -8.8461832117e-03];
PP(6).s=7;  PP(6).t1=[-6.2539373369e-03, 0, -8.0829196897e-03]; PP(6).t2=[-6.8759604862e-03, 0, -9.3913319650e-03]; PP(6).y2=[-6.0287353063e-03, -1.4213333333e-03, -8.2501054585e-03];
PP(7).s=8;  PP(7).t1=[-7.0881303767e-03, 0, -7.4636304039e-03]; PP(7).t2=[-7.7965132127e-03, 0, -8.9536984172e-03]; PP(7).y2=[-6.8316620229e-03, -1.6186666667e-03, -7.6540277053e-03];
PP(8).s=10; PP(8).t1=[-8.7565164562e-03, 0, -6.2250518322e-03]; PP(8).t2=[-9.6376186658e-03, 0, -8.0784313216e-03]; PP(8).y2=[-8.4375154560e-03, -2.0133333333e-03, -6.4618721989e-03];

%% 2-3. Load PATH data, compute vectors, decompose B
% Storage for t_hat direction paths
Btau_t = cell(n_s,1); Bn_t = cell(n_s,1); Bsum_t = cell(n_s,1);
dd_t   = cell(n_s,1); ang_t = cell(n_s,1);
% Storage for y_hat direction paths
Btau_y = cell(n_s,1); Bn_y = cell(n_s,1); Bsum_y = cell(n_s,1);
dd_y   = cell(n_s,1); ang_y = cell(n_s,1);

fprintf('\n--- Iron interior PATH decomposition ---\n');
for i = 1:n_s
    s_mm = s_vals(i);

    % --- t_hat direction ---
    fname = fullfile(results_dir, sprintf('path_d_t_s%02d.dat', s_mm));
    [d, BX, BY, BZ, BS] = load_prpath(fname);
    npts = length(d);
    frac = d / max(d);  % 0=axis, 1=surface
    p1 = PP(i).t1(:); p2 = PP(i).t2(:);

    Bt_i = zeros(npts,1); Bn_i = zeros(npts,1); ag_i = zeros(npts,1);
    for j = 1:npts
        pt = p1 + frac(j)*(p2 - p1);
        tau_d = (tip_ans - pt); tau_d = tau_d/norm(tau_d);
        B = [BX(j); BY(j); BZ(j)];
        Bt_i(j) = dot(B, tau_d);
        Bn_i(j) = dot(B, n_hat_t);
        ag_i(j) = acosd(abs(dot(tau_d, n_hat_t)));
    end
    Btau_t{i} = Bt_i; Bn_t{i} = Bn_i; Bsum_t{i} = BS;
    dd_t{i} = frac; ang_t{i} = ag_i;

    % --- y_hat direction ---
    fname = fullfile(results_dir, sprintf('path_d_y_s%02d.dat', s_mm));
    [d, BX, BY, BZ, BS] = load_prpath(fname);
    npts = length(d);
    frac = d / max(d);
    p1y = PP(i).t1(:); p2y = PP(i).y2(:);  % same axis point, different surface

    Bt_i = zeros(npts,1); Bn_i = zeros(npts,1); ag_i = zeros(npts,1);
    for j = 1:npts
        pt = p1y + frac(j)*(p2y - p1y);
        tau_d = (tip_ans - pt); tau_d = tau_d/norm(tau_d);
        B = [BX(j); BY(j); BZ(j)];
        Bt_i(j) = dot(B, tau_d);
        Bn_i(j) = dot(B, n_hat_y);
        ag_i(j) = acosd(abs(dot(tau_d, n_hat_y)));
    end
    Btau_y{i} = Bt_i; Bn_y{i} = Bn_i; Bsum_y{i} = BS;
    dd_y{i} = frac; ang_y{i} = ag_i;

    fprintf('  s=%2d mm: Btau=%.0f~%.0f, Bn=%.0f~%.0f mT, angle=%.1f~%.1f deg\n', ...
        s_mm, min(Bt_i)*1e3, max(Btau_t{i})*1e3, ...
        min(Bn_t{i})*1e3, max(Bn_t{i})*1e3, min(ang_t{i}), max(ang_t{i}));
end

%% 4. Iron interior figures
colors = lines(n_s);

% --- Fig 1: B_tau and B_n vs d ---
fig1 = figure('Position', [30 50 1100 450], 'Color', 'w');

subplot(1,2,1); hold on;
for i = 1:n_s
    plot(dd_t{i}, Btau_t{i}*1e3, '-', 'LineWidth', LINE_MAIN, 'Color', colors(i,:));
end
xlabel('d / d_{max}', 'FontSize', FONT_LABEL);
ylabel('B_\tau [mT]', 'FontSize', FONT_LABEL);
title('B_\tau  (toward tip)', 'FontSize', FONT_TITLE);
grid on; set(gca, 'FontSize', 11);

subplot(1,2,2); hold on;
for i = 1:n_s
    plot(dd_t{i}, Bn_t{i}*1e3, '-', 'LineWidth', LINE_MAIN, 'Color', colors(i,:));
end
xlabel('d / d_{max}', 'FontSize', FONT_LABEL);
ylabel('B_n [mT]', 'FontSize', FONT_LABEL);
title('B_n  (perpendicular to surface)', 'FontSize', FONT_TITLE);
grid on; set(gca, 'FontSize', 11);
legend(arrayfun(@(x) sprintf('s=%dmm',x), s_vals, 'Uni', 0), ...
    'Location', 'best', 'FontSize', 8);

exportgraphics(fig1, fullfile(fig_dir, 'P2_cone_Btau_Bn_vs_d.png'), 'Resolution', DPI);
fprintf('Saved Fig 1: P2_cone_Btau_Bn_vs_d.png\n');

% --- Fig 2: Energy fraction vs d ---
fig2 = figure('Position', [30 530 1000 450], 'Color', 'w');
hold on;
for i = 1:n_s
    pct_tau = (Btau_t{i} ./ Bsum_t{i}).^2 * 100;
    pct_n   = (Bn_t{i}   ./ Bsum_t{i}).^2 * 100;
    plot(dd_t{i}, pct_tau, '-',  'LineWidth', LINE_MAIN, 'Color', colors(i,:));
    plot(dd_t{i}, pct_n,   '--', 'LineWidth', LINE_THIN, 'Color', colors(i,:));
end
yline(sin(beta)^2*100, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(0.02, sin(beta)^2*100 + 2, sprintf('sin^2(\\beta) = %.1f%%', sin(beta)^2*100), ...
    'FontSize', 10, 'Color', [0.5 0.5 0.5]);
xlabel('d / d_{max}   (0 = axis, 1 = surface)', 'FontSize', FONT_LABEL);
ylabel('(B_i / |B|)^2  [%]', 'FontSize', FONT_LABEL);
title('Energy fraction: B_\tau^2 (solid) + B_n^2 (dashed)', 'FontSize', FONT_TITLE);
grid on; set(gca, 'FontSize', 11); ylim([0 105]);
text(0.02, 92, 'B_\tau^2 / |B|^2', 'FontSize', 11, 'Color', [0.3 0.3 0.3]);
text(0.02, 20, 'B_n^2 / |B|^2', 'FontSize', 11, 'Color', [0.3 0.3 0.3]);
legend(arrayfun(@(x) sprintf('s=%dmm',x), s_vals, 'Uni', 0), ...
    'Location', 'east', 'FontSize', 8);

exportgraphics(fig2, fullfile(fig_dir, 'P2_cone_energy_fraction_vs_d.png'), 'Resolution', DPI);
fprintf('Saved Fig 2: P2_cone_energy_fraction_vs_d.png\n');

% --- Fig 3: Angle between tau and n vs d ---
fig3 = figure('Position', [30 50 800 400], 'Color', 'w');
hold on;
for i = 1:n_s
    plot(dd_t{i}, ang_t{i}, '-', 'LineWidth', LINE_MAIN, 'Color', colors(i,:));
end
yline(90, ':', '90 deg', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xlabel('d / d_{max}', 'FontSize', FONT_LABEL);
ylabel('angle(tau, n) [deg]', 'FontSize', FONT_LABEL);
title('Angle between tau(d) and n  (90 deg = orthogonal)', 'FontSize', FONT_TITLE);
grid on; set(gca, 'FontSize', 11); ylim([75 92]);
legend(arrayfun(@(x) sprintf('s=%dmm',x), s_vals, 'Uni', 0), ...
    'Location', 'best', 'FontSize', 8);

exportgraphics(fig3, fullfile(fig_dir, 'P2_cone_angle_tau_n_vs_d.png'), 'Resolution', DPI);
fprintf('Saved Fig 3: P2_cone_angle_tau_n_vs_d.png\n');

%% 5. Load full-mesh sph10 data for air-side analysis
fprintf('\n--- Loading full-mesh sph10 data for air-side analysis ---\n');
data = import_ansys_data(fullfile('..', 'results', 'coil5_sph10'), 'all', 'coil5_sph10');

x_wp = data.x;
y_wp = data.y;
z_wp = data.z - cnst.SPH_OFST;

% Cone coordinates relative to P2
tip = [cnst.pole_tip_x(pole_idx); cnst.pole_tip_y(pole_idx); cnst.pole_tip_z_wp(pole_idx)];
vx = x_wp - tip(1);
vy = y_wp - tip(2);
vz = z_wp - tip(3);

s_all   = vx*ax(1) + vy*ax(2) + vz*ax(3);            % axial distance
dist2   = vx.^2 + vy.^2 + vz.^2;
r_perp  = sqrt(max(0, dist2 - s_all.^2));             % perpendicular distance
r_cone  = cnst.POLE_TIP_R + s_all * (cnst.POLE_R - cnst.POLE_TIP_R) / cnst.POLE_CONE_LEN;
dr      = r_perp - r_cone;                             % positive = air side
bmag    = data.bsum;

% Limit to fine-mesh sphere
r_wp = sqrt(x_wp.^2 + y_wp.^2 + z_wp.^2);
in_sphere = r_wp < 10e-3;  % SPH_FINE_R = 10mm for sph10
fprintf('  Nodes in fine sphere: %d / %d\n', sum(in_sphere), numel(r_wp));

%% 6. Binned iron-shell and air-shell analysis near P2 surface
s_edges  = (1:0.5:12) * 1e-3;  % s bins: 1 to 12 mm, 0.5 mm width
n_sbins  = length(s_edges) - 1;

% Storage
iron_Btau = zeros(n_sbins,1); iron_Bn = zeros(n_sbins,1); iron_Bmag = zeros(n_sbins,1);
air_Btau  = zeros(n_sbins,1); air_Bn  = zeros(n_sbins,1); air_Bmag  = zeros(n_sbins,1);
iron_nn   = zeros(n_sbins,1); air_nn  = zeros(n_sbins,1);
bin_s     = zeros(n_sbins,1);

fprintf('\n--- Surface shell analysis (iron & air) ---\n');
for b = 1:n_sbins
    in_bin = in_sphere & (s_all >= s_edges(b)) & (s_all < s_edges(b+1));
    idx_bin = find(in_bin);
    if isempty(idx_bin), continue; end

    dr_bin   = dr(idx_bin);
    bmag_bin = bmag(idx_bin);
    bin_s(b) = (s_edges(b) + s_edges(b+1)) / 2;

    % --- Iron shell: dr < 0, |B| > 50 mT, closest 20% ---
    iron_mask = (dr_bin < 0) & (bmag_bin > 50e-3) & (s_all(idx_bin) > 0);
    if sum(iron_mask) >= 3
        iron_sub = idx_bin(iron_mask);
        [~, si] = sort(abs(dr(iron_sub)));
        n_take = max(3, round(0.2 * length(si)));
        shell = iron_sub(si(1:min(n_take, length(si))));
        [bt, bn, bm] = compute_shell_Btau_Bn(data, shell, tip, ax, beta, cnst);
        iron_Btau(b) = mean(abs(bt));
        iron_Bn(b)   = mean(abs(bn));
        iron_Bmag(b)  = mean(bm);
        iron_nn(b)   = length(shell);
    end

    % --- Air shell: dr > 0, |B| < 50 mT, closest 20% ---
    air_mask = (dr_bin > 0) & (bmag_bin < 50e-3) & (s_all(idx_bin) > 0);
    if sum(air_mask) >= 3
        air_sub = idx_bin(air_mask);
        [~, si] = sort(dr(air_sub));
        n_take = max(3, round(0.2 * length(si)));
        shell = air_sub(si(1:min(n_take, length(si))));
        [bt, bn, bm] = compute_shell_Btau_Bn(data, shell, tip, ax, beta, cnst);
        air_Btau(b) = mean(abs(bt));
        air_Bn(b)   = mean(abs(bn));
        air_Bmag(b)  = mean(bm);
        air_nn(b)   = length(shell);
    end
end

valid = iron_nn > 0 & air_nn > 0;
fprintf('  Valid bins: %d / %d\n', sum(valid), n_sbins);

%% 7. Iron-to-air transition figure (KEY)
fig4 = figure('Position', [30 50 1200 500], 'Color', 'w');

% Left: |B_tau| and |B_n| on both sides vs s
subplot(1,2,1); hold on;
plot(bin_s(valid)*1e3, iron_Btau(valid)*1e3, '-o', 'LineWidth', LINE_MAIN, ...
    'Color', [0.85 0.15 0.15], 'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerSize', 4);
plot(bin_s(valid)*1e3, iron_Bn(valid)*1e3, '--s', 'LineWidth', LINE_THIN, ...
    'Color', [0.85 0.15 0.15], 'MarkerFaceColor', 'w', 'MarkerSize', 4);
plot(bin_s(valid)*1e3, air_Bn(valid)*1e3, '-o', 'LineWidth', LINE_MAIN, ...
    'Color', [0.15 0.3 0.85], 'MarkerFaceColor', [0.15 0.3 0.85], 'MarkerSize', 4);
plot(bin_s(valid)*1e3, air_Btau(valid)*1e3, '--s', 'LineWidth', LINE_THIN, ...
    'Color', [0.15 0.3 0.85], 'MarkerFaceColor', 'w', 'MarkerSize', 4);
xlabel('s [mm] (distance from tip)', 'FontSize', FONT_LABEL);
ylabel('|B| [mT]', 'FontSize', FONT_LABEL);
title('B components at surface (iron vs air)', 'FontSize', FONT_TITLE);
legend({'Iron |B_\tau|', 'Iron |B_n|', 'Air |B_n|', 'Air |B_\tau|'}, ...
    'FontSize', 9, 'Location', 'northeast');
grid on; set(gca, 'FontSize', 11);
set(gca, 'YScale', 'log');

% Right: B_n^2/|B|^2 on both sides
subplot(1,2,2); hold on;
iron_pct = (iron_Bn(valid) ./ iron_Bmag(valid)).^2 * 100;
air_pct  = (air_Bn(valid)  ./ air_Bmag(valid)).^2  * 100;
plot(bin_s(valid)*1e3, iron_pct, '-o', 'LineWidth', LINE_MAIN, ...
    'Color', [0.85 0.15 0.15], 'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerSize', 5);
plot(bin_s(valid)*1e3, air_pct, '-o', 'LineWidth', LINE_MAIN, ...
    'Color', [0.15 0.3 0.85], 'MarkerFaceColor', [0.15 0.3 0.85], 'MarkerSize', 5);
yline(50, ':', 'Color', [0.5 0.5 0.5]);
xlabel('s [mm]', 'FontSize', FONT_LABEL);
ylabel('B_n^2 / |B|^2  [%]', 'FontSize', FONT_LABEL);
title('Surface-normal energy fraction', 'FontSize', FONT_TITLE);
legend({'Iron side', 'Air side'}, 'FontSize', 10, 'Location', 'best');
grid on; set(gca, 'FontSize', 11); ylim([0 100]);

exportgraphics(fig4, fullfile(fig_dir, 'P2_cone_iron_air_transition_vs_s.png'), 'Resolution', DPI);
fprintf('Saved Fig 4: P2_cone_iron_air_transition_vs_s.png\n');

%% 8. Print summary
fprintf('\n========== SUMMARY ==========\n');
fprintf('Iron interior (PATH data, d=0 to surface):\n');
fprintf('  B_tau dominates: B_tau^2/|B|^2 ~ 90-100%%\n');
fprintf('  B_n small: B_n^2/|B|^2 ~ 3-15%%\n');
fprintf('  tau(d) and n non-orthogonal at interior (angle 78.8-90 deg)\n');
fprintf('\nIron-air transition at surface:\n');
fprintf('  Iron side: mean |B_tau| = %.1f mT, mean |B_n| = %.1f mT\n', ...
    mean(iron_Btau(valid))*1e3, mean(iron_Bn(valid))*1e3);
fprintf('  Air side:  mean |B_tau| = %.1f mT, mean |B_n| = %.1f mT\n', ...
    mean(air_Btau(valid))*1e3, mean(air_Bn(valid))*1e3);
fprintf('  B_tau ratio (iron/air) ~ %.0fx\n', ...
    mean(iron_Btau(valid)) / mean(air_Btau(valid)));
fprintf('  B_n ratio (iron/air) ~ %.1fx (should be ~1 if continuous)\n', ...
    mean(iron_Bn(valid)) / mean(air_Bn(valid)));
fprintf('\nAll figures saved to %s\n', fig_dir);
fprintf('Done.\n');


%% ========== Local function ==========
function [Btau_arr, Bn_arr, Bmag_arr] = compute_shell_Btau_Bn(data, node_idx, tip_wp, ax, beta, cnst)
% Compute B_tau (toward tip) and B_n (surface normal) for selected nodes
    x = data.x(node_idx);
    y = data.y(node_idx);
    z = data.z(node_idx) - cnst.SPH_OFST;  % WP frame
    bx = data.bx(node_idx);
    by = data.by(node_idx);
    bz = data.bz(node_idx);

    n = length(node_idx);
    Btau_arr = zeros(n,1);
    Bn_arr   = zeros(n,1);
    Bmag_arr = sqrt(bx.^2 + by.^2 + bz.^2);

    for k = 1:n
        pt = [x(k); y(k); z(k)];
        B  = [bx(k); by(k); bz(k)];

        % tau: toward tip
        v_to_tip = tip_wp - pt;
        tau_d = v_to_tip / norm(v_to_tip);
        Btau_arr(k) = dot(B, tau_d);

        % n: local surface normal (depends on azimuthal direction of node)
        v_from_tip = pt - tip_wp;
        s_ax = dot(v_from_tip, ax);
        v_perp = v_from_tip - s_ax * ax;
        rp = norm(v_perp);
        if rp > 1e-10
            r_hat_local = v_perp / rp;
        else
            r_hat_local = [0;0;0];  % on axis, no defined r_hat
        end
        n_hat_local = cos(beta) * r_hat_local - sin(beta) * ax;
        n_hat_local = n_hat_local / max(norm(n_hat_local), 1e-10);
        Bn_arr(k) = dot(B, n_hat_local);
    end
end
