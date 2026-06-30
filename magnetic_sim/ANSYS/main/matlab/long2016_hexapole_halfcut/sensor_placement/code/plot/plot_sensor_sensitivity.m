%% plot_sensor_sensitivity.m
%  ------------------------------------------------------------------
%  Sensitivity curve: how does |V_out/V_in| diagonal change with
%  sensor placement distance? Demonstrates non-linearity.
%
%  Sources (.mat files saved by gen_Vout_Vin_*.m chain):
%    Vout_Vin_kA0p36.mat   (sensor 4.572 mm — baseline)
%    Vout_Vin_4p472.mat    (sensor 4.472 mm)
%    Vout_Vin_3p572.mat    (sensor 3.572 mm)
%    Vout_Vin_3p472.mat    (sensor 3.472 mm)
%    Vout_Vin_2p286.mat    (sensor 2.286 mm)
%
%  Output: magnetic_sim/ANSYS/main/figures/long2016_hexapole_halfcut/sensor_sensitivity.png
%  ------------------------------------------------------------------

clear; close all;

data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\bs_matrix\data';

% 4 new variants (loadable from .mat — V_VV stored in true V/V dimensionless)
new_variants = struct( ...
    'sensor_mm', {4.472, 3.572, 3.472, 2.286}, ...
    'file',      {'Vout_Vin_4p472.mat', 'Vout_Vin_3p572.mat', ...
                  'Vout_Vin_3p472.mat', 'Vout_Vin_2p286.mat'});

% Baseline 4.572 mm: hardcode from memory `project_long_fei_B_bar`
% (Vout_Vin_kA0p36.mat has a unit-label bug: stores values in mV/V scale
% but labels as "dimensionless"; safer to use the known-correct values)
baseline_sensor_mm = 4.572;
baseline_diag_mVV  = [1.0978, 1.5763, 1.0600, 1.6551, 1.5723, 1.2005];

N = 1 + length(new_variants);
sensor_mm   = zeros(N, 1);
diag_per_pole = zeros(N, 6);   % |V_VV(i,i)| in mV/V

% Baseline
sensor_mm(1)        = baseline_sensor_mm;
diag_per_pole(1, :) = baseline_diag_mVV;

% New variants (load .mat; V_VV is V/V dimensionless, * 1e3 -> mV/V)
for k = 1:length(new_variants)
    s = load(fullfile(data_dir, new_variants(k).file));
    sensor_mm(1+k)        = new_variants(k).sensor_mm;
    diag_per_pole(1+k, :) = abs(diag(s.V_VV)) * 1e3;
end

% Sort ascending by sensor distance for clean plot
[sensor_mm, idx]  = sort(sensor_mm);
diag_per_pole     = diag_per_pole(idx, :);
mean_diag         = mean(diag_per_pole, 2);

% Local sensitivity |dE/dΔ| between adjacent points
sensor_mid = 0.5 * (sensor_mm(1:end-1) + sensor_mm(2:end));
d_diag = diff(mean_diag);
d_sensor = diff(sensor_mm);
local_S = abs(d_diag ./ d_sensor);

% Chord-slope sensitivity from baseline (4.572)
baseline_idx  = find(abs(sensor_mm - 4.572) < 1e-6, 1);
baseline_mean = mean_diag(baseline_idx);
chord_S = zeros(N, 1);
for k = 1:N
    dx = sensor_mm(k) - sensor_mm(baseline_idx);
    if abs(dx) < 1e-6
        chord_S(k) = NaN;
    else
        chord_S(k) = abs(mean_diag(k) - baseline_mean) / abs(dx);
    end
end

% Echo to console
fprintf('Sensor (mm)   mean |V/V| (mV/V)   chord S   local S (midpt)\n');
for k = 1:N
    if isnan(chord_S(k))
        cs_str = '   —   ';
    else
        cs_str = sprintf('%6.3f', chord_S(k));
    end
    if k <= length(local_S)
        ls_str = sprintf('%6.3f @ %.3f', local_S(k), sensor_mid(k));
    else
        ls_str = '';
    end
    fprintf('  %6.3f       %7.3f           %s    %s\n', ...
            sensor_mm(k), mean_diag(k), cs_str, ls_str);
end

%% --- Figure: 2-panel layout ---
fig = figure('Color', 'w', 'Position', [80 80 1100 460]);

% --- Left: mean |V/V| vs sensor distance, with per-pole faint lines ---
ax1 = subplot(1, 2, 1);
hold(ax1, 'on');

pole_colors = lines(6);
for p = 1:6
    plot(ax1, sensor_mm, diag_per_pole(:, p), '-', ...
         'Color', [pole_colors(p,:), 0.30], 'LineWidth', 0.8, ...
         'HandleVisibility', 'off');
    plot(ax1, sensor_mm, diag_per_pole(:, p), 'o', ...
         'MarkerSize', 4, 'MarkerFaceColor', pole_colors(p,:), ...
         'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
end
h_mean = plot(ax1, sensor_mm, mean_diag, 'k-', 'LineWidth', 2.0, ...
              'DisplayName', 'mean over 6 poles');
plot(ax1, sensor_mm, mean_diag, 'ko', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'k', 'HandleVisibility', 'off');

% Baseline reference
xline(ax1, 4.572, '--', '4.572 mm baseline', ...
      'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
      'LabelHorizontalAlignment', 'left', ...
      'LabelVerticalAlignment', 'bottom');

% Linear-extrapolation reference (from baseline using small-Δ slope)
S_baseline = abs(diff(mean_diag(end-1:end))) / abs(diff(sensor_mm(end-1:end)));   % not great
% Use the chord between baseline and adjacent point (4.472) as "small-Δ" estimate
near_idx = find(sensor_mm < 4.572, 1, 'last');
S_lin = abs(mean_diag(near_idx) - baseline_mean) / abs(sensor_mm(near_idx) - 4.572);
x_lin = linspace(min(sensor_mm), 4.572, 50);
y_lin = baseline_mean - S_lin * (4.572 - x_lin);
plot(ax1, x_lin, y_lin, '--', 'Color', [0.85 0.30 0.30], 'LineWidth', 1.2, ...
     'DisplayName', sprintf('linear extrap (S=%.2f)', S_lin));

grid(ax1, 'on'); box(ax1, 'on');
xlabel(ax1, 'Sensor position (mm along surface from tip)', 'FontSize', 11);
ylabel(ax1, 'mean |V_{out}/V_{in}| diagonal  (mV/V)', 'FontSize', 11);
title(ax1, 'Signal magnitude vs sensor position', 'FontSize', 12);
legend(ax1, 'Location', 'northwest', 'FontSize', 10);
xlim(ax1, [min(sensor_mm)-0.2, max(sensor_mm)+0.2]);
ylim(ax1, [0, max(mean_diag)*1.15]);

% --- Right: local sensitivity (rate of change) ---
ax2 = subplot(1, 2, 2);
hold(ax2, 'on');

% Bar chart of local |dE/dΔ| over each segment
for k = 1:length(local_S)
    bar(ax2, sensor_mid(k), local_S(k), 0.18, ...
        'FaceColor', [0.20 0.40 0.70], 'EdgeColor', 'none');
    text(ax2, sensor_mid(k), local_S(k) + 0.03, ...
         sprintf('%.2f', local_S(k)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9, ...
         'Color', [0.15 0.30 0.55]);
end

% Chord-slope overlay (from baseline)
chord_keep = ~isnan(chord_S);
plot(ax2, sensor_mm(chord_keep), chord_S(chord_keep), 'r-s', ...
     'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'chord S (from 4.572)');

grid(ax2, 'on'); box(ax2, 'on');
xlabel(ax2, 'Sensor position (mm)', 'FontSize', 11);
ylabel(ax2, 'sensitivity  (mV/V per mm)', 'FontSize', 11);
title(ax2, 'Sensitivity: |\DeltaV/V/V| / |\Deltasensor|', 'FontSize', 12);
legend(ax2, {'local |dE/d\Delta| (per segment)', 'chord S (from baseline)'}, ...
       'Location', 'northwest', 'FontSize', 10);
xlim(ax2, [min(sensor_mm)-0.2, max(sensor_mm)+0.2]);
ylim(ax2, [0, max([local_S; chord_S(~isnan(chord_S))]) * 1.30]);

sgtitle(fig, 'Hall sensor placement sensitivity — Long Fei hexapole halfcut (FEM)', ...
        'FontSize', 13, 'FontWeight', 'bold');

%% --- Save ---
out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
out_path = fullfile(out_dir, 'sensor_sensitivity.png');
exportgraphics(fig, out_path, 'Resolution', 200);
fprintf('\nSaved: %s\n', out_path);
