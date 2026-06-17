%% match_target_matrix.m
%  ------------------------------------------------------------------
%  Compare a user-supplied target B-matrix against the 7 sensor variants
%  computed by gen_Vout_Vin_*.m, find the closest sensor distance.
%
%  3 distance metrics:
%   (1) raw element-wise RMSE (same units)
%   (2) shape similarity (after L2 normalisation; ignores absolute scale)
%   (3) best-fit scale factor + residual
%  ------------------------------------------------------------------

clear; clc;

%% --- Target matrix (from user image, 2026-05-25) ---
target = [ ...
   0.2816  -0.0052  -0.0406  -0.0305  -0.0316  -0.0505;
  -0.0111   0.2659  -0.0423  -0.0666  -0.0756  -0.0346;
  -0.0470  -0.0332   0.2795  -0.0064  -0.0237  -0.0550;
  -0.0308  -0.0646  -0.0129   0.2348  -0.0833  -0.0193;
  -0.0411  -0.0661  -0.0217  -0.0760   0.2537  -0.0118;
  -0.0477  -0.0208  -0.0402  -0.0250  -0.0056   0.2302 ];

fprintf('Target matrix:\n  diag mean = %.4f\n  off-diag mean |.| = %.4f\n\n', ...
        mean(diag(target)), ...
        mean(abs(target(~eye(6)))));

%% --- Load 7 sensor variants (V_VV is V/V dimensionless in new files) ---
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\bs_matrix';

% Baseline 4.572: hardcode from memory (Vout_Vin_kA0p36.mat has unit-label bug)
% diag in mV/V; full matrix from existing Vout_Vin_kA0p36.tex
baseline_full_mVV = [ ...
  +1.0978  -0.0230  -0.0395  -0.5845  -0.5789  -0.0402;
  -0.1182  +1.5763  -0.1917  -0.5672  -0.5696  -0.1928;
  -0.0455  -0.5634  +1.0600  -0.0199  -0.5609  -0.0452;
  -0.1992  -0.5844  -0.1288  +1.6551  -0.5891  -0.2092;
  -0.1905  -0.5632  -0.1798  -0.5711  +1.5723  -0.1166;
  -0.0284  -0.6245  -0.0274  -0.6166  -0.0327  +1.2005 ];

variants = struct( ...
    'sensor_mm', {4.572, 4.472, 3.572, 3.472, 2.572, 2.286, 1.572}, ...
    'file',      {'',    'Vout_Vin_4p472.mat', 'Vout_Vin_3p572.mat', ...
                  'Vout_Vin_3p472.mat', 'Vout_Vin_2p572.mat', ...
                  'Vout_Vin_2p286.mat', 'Vout_Vin_1p572.mat'});

N = length(variants);
mats_mVV = cell(N, 1);   % each matrix in mV/V scale, cosmetic applied

for k = 1:N
    if k == 1
        Vm_mVV = baseline_full_mVV;
    else
        s = load(fullfile(data_dir, variants(k).file));
        Vm_mVV = s.V_VV * 1e3;   % V/V -> mV/V
    end
    % Apply cosmetic: diag +|.|, off-diag -|.|
    for i = 1:6
        for j = 1:6
            if i == j
                Vm_mVV(i,j) =  abs(Vm_mVV(i,j));
            else
                Vm_mVV(i,j) = -abs(Vm_mVV(i,j));
            end
        end
    end
    mats_mVV{k} = Vm_mVV;
end

%% --- Compute 3 metrics for each variant ---
norm_target = norm(target, 'fro');
target_n    = target / norm_target;          % unit Frobenius

results = zeros(N, 5);   % [sensor, RMSE_raw, shape_dist, alpha, residual_rel]

for k = 1:N
    M = mats_mVV{k};

    % (1) Element-wise RMSE in mV/V (target assumed in mV/V for this metric)
    RMSE_raw = sqrt(mean((M(:) - target(:)).^2));

    % (2) Shape similarity: ||M/||M|| - target/||target|||| in Frobenius
    M_n        = M / norm(M, 'fro');
    shape_dist = norm(M_n - target_n, 'fro');

    % (3) Best-fit scalar α: minimize ||α·M - target||_F
    alpha      = sum(M(:).*target(:)) / sum(M(:).*M(:));
    residual   = norm(alpha*M - target, 'fro') / norm_target;

    results(k, :) = [variants(k).sensor_mm, RMSE_raw, shape_dist, alpha, residual];
end

%% --- Sort & print ---
fprintf('Comparison: 7 sensor variants vs target matrix\n');
fprintf('%s\n', repmat('=', 1, 88));
fprintf('%-12s %12s %12s %12s %12s\n', ...
        'sensor (mm)', 'RMSE raw', 'shape dist', 'best alpha', 'residual %');
fprintf('%-12s %12s %12s %12s %12s\n', ...
        '',            '(mV/V)',   '(0 = same)', '(scale fac)', '(after scale)');
fprintf('%s\n', repmat('-', 1, 88));
for k = 1:N
    fprintf('%-12.3f %12.4f %12.4f %12.4f %12.2f\n', ...
            results(k,1), results(k,2), results(k,3), results(k,4), results(k,5)*100);
end
fprintf('%s\n', repmat('=', 1, 88));

%% --- Winners by each metric ---
[~, idx_rmse]  = min(results(:, 2));
[~, idx_shape] = min(results(:, 3));
[~, idx_resid] = min(results(:, 5));

fprintf('\nWinners:\n');
fprintf('  (1) closest in raw value (RMSE)       : sensor = %.3f mm\n', results(idx_rmse, 1));
fprintf('  (2) closest in shape (scale-blind)    : sensor = %.3f mm (shape dist %.3f)\n', ...
        results(idx_shape, 1), results(idx_shape, 3));
fprintf('  (3) closest after best-fit rescale    : sensor = %.3f mm (residual %.1f%%, alpha=%.3f)\n', ...
        results(idx_resid, 1), results(idx_resid, 5)*100, results(idx_resid, 4));

%% --- Diagnose: is alpha consistent across variants? ---
fprintf('\nDiagnostic — best-fit alpha across 7 variants:\n');
fprintf('  range: [%.3f, %.3f]   mean: %.3f   std: %.3f\n', ...
        min(results(:,4)), max(results(:,4)), mean(results(:,4)), std(results(:,4)));
if std(results(:,4)) / mean(results(:,4)) < 0.20
    fprintf('  -> alpha STABLE: FEM shape matches target; absolute scale offset\n');
    fprintf('     is a uniform gain factor (different hardware sensitivity).\n');
else
    fprintf('  -> alpha NOT stable: FEM-vs-target mismatch is sensor-position-dependent;\n');
    fprintf('     no single rescale makes all variants fit equally well.\n');
end

%% --- Side-by-side diag check for top match ---
best_idx = idx_resid;
fprintf('\nDiagonal comparison for best match (sensor %.3f mm):\n', results(best_idx,1));
fprintf('  Pole    target     FEM (mV/V)    FEM x alpha    relative err %%\n');
M_best = mats_mVV{best_idx};
alpha_best = results(best_idx, 4);
for i = 1:6
    fem_val = M_best(i,i);
    fem_scaled = alpha_best * fem_val;
    rel_err = (fem_scaled - target(i,i)) / target(i,i) * 100;
    fprintf('  P%d     %+8.4f    %+8.4f      %+8.4f         %+6.2f %%\n', ...
            i, target(i,i), fem_val, fem_scaled, rel_err);
end
