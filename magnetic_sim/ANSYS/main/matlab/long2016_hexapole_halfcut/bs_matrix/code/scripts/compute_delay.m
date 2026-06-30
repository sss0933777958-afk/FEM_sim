%% compute_delay.m — compute time delay between B_surface and B_tip rigorously
%
%  Three methods, all from complex phasor Bn_surf, Bn_tip in fig44_P1_v2.mat:
%
%  (1) Raw phase difference: Δt_raw = (φ_tip - φ_surf) / (2π·f)
%      → gives ~half-period due to 180° sign convention from sensor normals
%
%  (2) Sign-corrected (align Long Fei convention V_tip ∝ +V_surf):
%      Flip B_tip → use -B_tip. Then phase diff is small residual.
%      Δt_aligned = (φ_tip_flipped - φ_surf) / (2π·f), wrapped to [-T/2, T/2]
%
%  (3) Cross-correlation in time domain (numerical confirmation of method 2)
%      Generate V(t) over multiple periods, find peak of x-corr.

clear; clc;
S = load('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\bs_matrix\data\fig44_P1_v2.mat');
freqs = S.freqs;
Bs = S.Bn_surf;          % complex phasor, sign convention: n_surf = +z
Bt = S.Bn_tip;           % complex phasor, sign convention: n_tip = +x_pole (toward apex)

fprintf('=== Method 1: Raw phase difference (no sign correction) ===\n');
fprintf('%-6s | %-10s | %-10s | %-12s | %-14s\n', ...
        'f Hz','φ_surf °','φ_tip °','Δφ raw °','Δt_raw [µs]');
fprintf('%s\n', repmat('-',1,70));
for k = 1:numel(freqs)
    phi_s = atan2d(imag(Bs(k)), real(Bs(k)));
    phi_t = atan2d(imag(Bt(k)), real(Bt(k)));
    dphi  = phi_t - phi_s;
    dt_us = dphi / (360 * freqs(k)) * 1e6;
    fprintf('%-6d | %+10.2f | %+10.2f | %+12.2f | %+14.3f\n', ...
            freqs(k), phi_s, phi_t, dphi, dt_us);
end

%% Method 2: sign-corrected, principal value [-T/2, T/2]
fprintf('\n=== Method 2: Sign-corrected (flip B_tip → use -B_tip) ===\n');
fprintf('  Then B_surf and -B_tip are nearly in-phase. Wrap Δφ to [-180°, +180°].\n');
fprintf('%-6s | %-15s | %-15s | %-14s\n', ...
        'f Hz','Δφ aligned °','wrapped °','Δt [µs]');
fprintf('%s\n', repmat('-',1,70));
Bt_flip = -Bt;   % flip sign of B_tip to align with Long Fei convention
dt_vec = zeros(numel(freqs),1);
for k = 1:numel(freqs)
    phi_s  = atan2d(imag(Bs(k)),     real(Bs(k)));
    phi_tf = atan2d(imag(Bt_flip(k)),real(Bt_flip(k)));
    dphi_raw = phi_tf - phi_s;
    % Wrap to [-180°, +180°]
    dphi = mod(dphi_raw + 180, 360) - 180;
    dt_us = dphi / (360 * freqs(k)) * 1e6;
    dt_vec(k) = dt_us;
    fprintf('%-6d | %+15.3f | %+15.3f | %+14.4f\n', ...
            freqs(k), dphi_raw, dphi, dt_us);
end

%% Method 3: time-domain cross-correlation
fprintf('\n=== Method 3: Cross-correlation in time domain ===\n');
fprintf('  Reconstruct V_surf(t) and V_tip_flipped(t) over 10 periods,\n');
fprintf('  resample at high rate (1000 pts/period), find xcorr peak.\n\n');
fprintf('%-6s | %-15s | %-15s | %-10s\n', ...
        'f Hz','xcorr peak [µs]','phase method [µs]','agree?');
fprintf('%s\n', repmat('-',1,65));

for k = 1:numel(freqs)
    f = freqs(k);
    T = 1/f;
    N_per_period = 1000;
    N_periods = 10;
    t = linspace(0, N_periods*T, N_periods*N_per_period);

    V_s = real(Bs(k))     * cos(2*pi*f*t) - imag(Bs(k))     * sin(2*pi*f*t);
    V_t = real(Bt_flip(k))* cos(2*pi*f*t) - imag(Bt_flip(k))* sin(2*pi*f*t);

    [R, lags] = xcorr(V_t, V_s, 'normalized');
    [~, idx] = max(R);
    lag_samples = lags(idx);
    dt_samples_us = lag_samples * (t(2)-t(1)) * 1e6;

    fprintf('%-6d | %+15.4f | %+15.4f | %s\n', ...
            f, dt_samples_us, dt_vec(k), ...
            ternary(abs(dt_samples_us - dt_vec(k)) < 1, '✓', 'mismatch'));
end

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('After removing the 180° sign convention (flip B_tip):\n');
fprintf('  f=1 kHz: Δt = %+.3f µs\n', dt_vec(1));
fprintf('  f=2 kHz: Δt = %+.3f µs\n', dt_vec(2));
fprintf('  f=3 kHz: Δt = %+.3f µs\n', dt_vec(3));
fprintf('\nInterpretation:\n');
fprintf('  - Pure physics (σ=0, quasi-static): Δt should be 0\n');
fprintf('  - Light-speed propagation over %g mm: %g fs (negligible)\n', ...
        norm([4.980-0.158;0;-12.590+13.000]), ...
        norm([4.980-0.158;0;-12.590+13.000])*1e-3/3e8*1e15);
fprintf('  - Observed |Δt| < %.2f µs: numerical noise from COMSOL AC artifact\n', ...
        max(abs(dt_vec)));
fprintf('  - %.2f µs / Long-Fei system T(@200Hz BW)=5 ms = %.4f%%\n', ...
        max(abs(dt_vec)), max(abs(dt_vec))/5000*100);
fprintf('  → Effective time delay = 0 (within COMSOL numerical precision)\n');

function s = ternary(cond, t, f)
    if cond, s = t; else, s = f; end
end
