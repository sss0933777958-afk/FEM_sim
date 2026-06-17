%% verify_B_perpendicularity.m — Check B ⊥ P1 pole surface (air side)
%
%  Verifies that B is perpendicular to the P1 cone surface on the air
%  side, as expected from the iron-air boundary condition (murx=280).
%
%  Approach:
%    1. Load Coil1 full-mesh FEM data (~494k nodes)
%    2. For each axial bin along the P1 cone (s = 1-6 mm from tip),
%       select the closest shell of confirmed-air nodes on the
%       LOWER half of the cone (phi < -10°, avoiding the VSBV flat cut)
%    3. Compute the analytical cone surface normal at each node
%    4. Measure the angle between B and the surface normal
%
%  Key definitions:
%    s     — axial distance from P1 tip along cone axis [m]
%    dr    — radial distance from analytical cone surface [m]
%            (positive = outside cone = air side)
%    phi   — azimuthal angle around cone axis [deg]
%            (phi < 0 = lower half, phi > 0 = upper half / flat-cut)
%    angle — deviation of B from surface normal [deg]
%            (0° = B perfectly perpendicular to surface)
%
%  Run from magnetic_sim/ANSYS/backup/hexapole-long2016/analysis/ directory.

clear; clc; close all;

%% 1. Load constants and full mesh data
cnst = mt_constants();

fprintf('Loading full mesh data for Coil1...\n');
data = import_ansys_data(fullfile('..', 'results', 'coil1'), 'all', 'coil1');

% Convert to WP frame
x_wp = data.x;
y_wp = data.y;
z_wp = data.z - cnst.SPH_OFST;

%% 2. Select nodes within fine mesh sphere (7 mm)
SPH_FINE_R = 7e-3;  % [m] — matches original APDL SPH_FINE_R (use 10e-3 for sph10 data)
r_wp = sqrt(x_wp.^2 + y_wp.^2 + z_wp.^2);
in_sphere = r_wp < SPH_FINE_R;
idx_sph = find(in_sphere);
fprintf('Nodes in fine sphere: %d / %d\n', sum(in_sphere), numel(r_wp));

%% 3. Compute cone coordinates for all sphere nodes
% P1: pole index 1, lower pole, angle 0°, axis = [1,0,0]
pole_idx = 1;
tip = [cnst.pole_tip_x(pole_idx); cnst.pole_tip_y(pole_idx); cnst.pole_tip_z_wp(pole_idx)];
ax = cnst.pole_axis(:, pole_idx);  % unit vector, tip -> base

% Vector from tip to each node
vx = x_wp(idx_sph) - tip(1);
vy = y_wp(idx_sph) - tip(2);
vz = z_wp(idx_sph) - tip(3);

% Axial projection (s > 0 = toward base)
s = vx*ax(1) + vy*ax(2) + vz*ax(3);

% Perpendicular distance from axis
dist = sqrt(vx.^2 + vy.^2 + vz.^2);
r_perp = sqrt(max(0, dist.^2 - s.^2));

% Analytical cone radius and distance from surface
r_cone = cnst.POLE_TIP_R + s * (cnst.POLE_R - cnst.POLE_TIP_R) / cnst.POLE_CONE_LEN;
dr = r_perp - r_cone;  % positive = air side

% B magnitude and azimuthal angle around cone axis
b_full = sqrt(data.bx(idx_sph).^2 + data.by(idx_sph).^2 + data.bz(idx_sph).^2);
phi = atan2d(vz, vy);  % azimuth in (y,z) plane perpendicular to axis

%% 4. Binned analysis: closest air shell per axial bin
s_edges = (1:0.25:6.25) * 1e-3;  % s = 1 to 6 mm
n_bins = length(s_edges) - 1;
beta = atan((cnst.POLE_R - cnst.POLE_TIP_R) / cnst.POLE_CONE_LEN);  % cone half-angle ~11.2°

% Storage
bin_s         = zeros(n_bins, 1);
bin_angle_mean = zeros(n_bins, 1);
bin_angle_med  = zeros(n_bins, 1);
bin_angle_std  = zeros(n_bins, 1);
bin_bmag_mean  = zeros(n_bins, 1);
bin_n_nodes    = zeros(n_bins, 1);
bin_dr_mean    = zeros(n_bins, 1);
bin_valid      = false(n_bins, 1);

for b = 1:n_bins
    in_bin = (s >= s_edges(b)) & (s < s_edges(b+1));
    idx_bin = find(in_bin);
    if isempty(idx_bin), continue; end

    dr_bin  = dr(idx_bin);
    b_bin   = b_full(idx_bin);
    phi_bin = phi(idx_bin);

    % Confirmed air nodes on lower half of cone (avoid VSBV flat cut)
    air_mask = (b_bin < 50e-3) & (dr_bin > 0) & (phi_bin < -10);
    if sum(air_mask) == 0, continue; end

    air_sub    = idx_bin(air_mask);
    dr_air_sub = dr(air_sub);

    % Adaptive shell: take the closest 20% of air nodes (min 3)
    n_take = max(3, round(0.2 * length(dr_air_sub)));
    [~, si] = sort(dr_air_sub);
    shell_nodes = air_sub(si(1:min(n_take, length(si))));

    % Compute angle for each node in shell
    angles = zeros(length(shell_nodes), 1);
    bmags  = zeros(length(shell_nodes), 1);
    for j = 1:length(shell_nodes)
        k  = shell_nodes(j);
        gi = idx_sph(k);

        v_perp_vec = [vx(k); vy(k); vz(k)] - s(k) * ax;
        rp = norm(v_perp_vec);
        if rp < 1e-10, angles(j) = NaN; continue; end

        rad_dir = v_perp_vec / rp;
        % Outward normal: cos(beta)*radial - sin(beta)*axial
        n_hat = cos(beta) * rad_dir - sin(beta) * ax;
        n_hat = n_hat / norm(n_hat);

        B    = [data.bx(gi); data.by(gi); data.bz(gi)];
        Bmag = norm(B);
        bmags(j)  = Bmag;
        ca = min(1, abs(dot(B, n_hat)) / Bmag);
        angles(j) = acosd(ca);
    end

    valid  = ~isnan(angles);
    angles = angles(valid);
    bmags  = bmags(valid);
    if isempty(angles), continue; end

    bin_s(b)          = (s_edges(b) + s_edges(b+1)) / 2;
    bin_angle_mean(b) = mean(angles);
    bin_angle_med(b)  = median(angles);
    bin_angle_std(b)  = std(angles);
    bin_bmag_mean(b)  = mean(bmags);
    bin_n_nodes(b)    = length(angles);
    bin_dr_mean(b)    = mean(dr(shell_nodes));
    bin_valid(b)      = true;
end

v = bin_valid;

%% 5. Print results
fprintf('\n===== B Perpendicularity on P1 Cone Surface (air side, lower half) =====\n');
fprintf('Cone half-angle beta: %.2f°\n', rad2deg(beta));
fprintf('Selection: phi < -10° (lower half, excludes VSBV flat cut)\n');
fprintf('Shell: closest 20%% of air nodes per bin (min 3)\n\n');

fprintf('%-8s  %6s  %8s  %8s  %8s  %8s  %10s\n', ...
    's[mm]', 'Nodes', 'Mean[°]', 'Med[°]', 'Std[°]', 'dr[um]', '|B|[mT]');
fprintf('%s\n', repmat('-', 1, 64));
for b = find(v')
    fprintf('%6.2f    %5d  %8.2f  %8.2f  %8.2f  %8.1f  %10.1f\n', ...
        bin_s(b)*1e3, bin_n_nodes(b), bin_angle_mean(b), bin_angle_med(b), ...
        bin_angle_std(b), bin_dr_mean(b)*1e6, bin_bmag_mean(b)*1e3);
end

fprintf('\nOverall (s=1-6mm): mean angle = %.2f°, range = %.2f° - %.2f°\n', ...
    mean(bin_angle_mean(v)), min(bin_angle_mean(v)), max(bin_angle_mean(v)));

%% 6. Figures

% --- Fig 1: Angle vs s (mean ± std) ---
figure('Name', 'B Perpendicularity: Angle vs s', 'Position', [100 100 900 450]);
errorbar(bin_s(v)*1e3, bin_angle_mean(v), bin_angle_std(v), ...
    '-o', 'LineWidth', 1.8, 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.2 0.4 0.8], 'Color', [0.2 0.4 0.8], 'CapSize', 5);
xlabel('Axial distance from tip  s [mm]', 'FontSize', 12);
ylabel('Angle from surface normal [°]', 'FontSize', 12);
title('B perpendicularity at closest air shell (P1 lower-half cone)', 'FontSize', 13);
grid on; xlim([0.9 6.2]);
yline(5, '--', '5°', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'FontSize', 11);
set(gca, 'FontSize', 11);

% --- Fig 2: |B| vs s ---
figure('Name', 'B Perpendicularity: |B| vs s', 'Position', [100 580 900 450]);
plot(bin_s(v)*1e3, bin_bmag_mean(v)*1e3, '-s', 'LineWidth', 1.8, 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.8 0.3 0.2], 'Color', [0.8 0.3 0.2]);
xlabel('Axial distance from tip  s [mm]', 'FontSize', 12);
ylabel('|B| [mT]', 'FontSize', 12);
title('|B| at closest air shell (P1 lower-half cone)', 'FontSize', 13);
grid on; xlim([0.9 6.2]);
set(gca, 'FontSize', 11);

fprintf('\nDone. 2 figures generated.\n');
