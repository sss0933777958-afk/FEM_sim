%% main.m -- Maxwell single-pole calibration driver (current base)
%  sp_config -> sp_load (.fld -> tip frame, mT, iron mask) -> Nr ladder of on-axis samples ->
%  sp_fit -> NMAE on every air node of the ball -> convergence -> save .mat -> sp_emit PDF.
%  Frame: origin = pole tip apex, pole axis +x; l_hat = distance along x from the tip to the charge.
clearvars;  clc;

%% ---- per-run settings ------------------------------------------------------
MODEL    = 'NTU_hexapole';
GEOM     = 'single_pole';
USE_BIAS = true;            % false = charge on the axis (single), true = + e_y, e_z (eighteen)
l0       = 400e-6;          % initial l_hat [m]
I_actual = 1;               % drive current [A] (= FEM excitation: 50 A-turn / N_c 50)
NRMAX    = 200;             % ladder safety stop

%% ---- paths -----------------------------------------------------------------
here = fileparts(mfilename('fullpath'));                    % .../Maxwell/main/single_pole
CAL  = fileparts(fileparts(here));                          % .../Maxwell
addpath(fullfile(CAL, 'function', 'single_pole'));

%% ---- load ------------------------------------------------------------------
cfg = sp_config(MODEL, GEOM);
assert(cfg.calib.fix_ex, 'main: only e_x = 0 is implemented');
R   = cfg.calib.R;
ctr = cfg.calib.center(:).';
L   = sp_load(cfg);

% evaluation cloud: every air grid node inside the ball
[X, Y, Z] = ndgrid(L.xs, L.ys, L.zs);
in  = (X-ctr(1)).^2 + (Y-ctr(2)).^2 + (Z-ctr(3)).^2 <= R^2 & ~L.iron;
P_ev = [X(in), Y(in), Z(in)];
B_ev = reshape(L.B, [], 3);  B_ev = B_ev(in(:), :);
fprintf('[eval] air nodes inside R = %g um: %d\n', R*1e6, size(P_ev,1));

%% ---- Nr ladder: sample -> fit -> NMAE -> convergence -----------------------
%  samples x_k = -k*R/Nr (k = 1..Nr) on the axis (air side); a sample is dropped when its
%  trilinear stencil touches an iron node (its value would mix in the steel field).
NMAE = nan(NRMAX,1);  NPTS = nan(NRMAX,1);  LH = nan(NRMAX,1);  Nc = NaN;
for n = 1:NRMAX
    [P, b] = samples(n, R, ctr, L);
    fit = sp_fit(P, b, l0, USE_BIAS);
    NMAE(n) = sum(vecnorm(B_ev - fit.predict(P_ev), 2, 2)) / sum(vecnorm(B_ev, 2, 2)) * 100;
    NPTS(n) = size(P,1);  LH(n) = fit.l_hat*1e6;
    K = cfg.calib.conv_win;
    if n > K
        ch = abs(diff(NMAE(n-K:n))) ./ NMAE(n-K:n-1);       % last K relative steps
        if all(ch < cfg.calib.conv_tol), Nc = n - K;  break; end   % window start
    end
end
assert(isfinite(Nc), 'main: NMAE did not converge within Nr = %d', NRMAX);

[P, b] = samples(Nc, R, ctr, L);
fit = sp_fit(P, b, l0, USE_BIAS);
fprintf('[conv] Nr = %d -> %d samples (ladder scanned to %d)\n', Nc, size(P,1), n);

%% ---- record + save ---------------------------------------------------------
rec = struct('base','current', 'model',MODEL, 'GEOM',GEOM, 'USE_BIAS',USE_BIAS, ...
             'l_hat',fit.l_hat, 'e',fit.e, 'G',fit.G, 'gI_hat',fit.G/I_actual, 'J',fit.J, ...
             'npts',size(P,1), 'Nr',Nc, 'P',P, 'b',b, 'NMAE',NMAE(Nc), 'Np_eval',size(P_ev,1), ...
             'NMAE_ladder',NMAE(1:n), 'npts_ladder',NPTS(1:n), 'l_hat_ladder',LH(1:n), ...
             'R_select',R, 'center',ctr, 'l0',l0, 'I_actual',I_actual, ...
             'conv_tol',cfg.calib.conv_tol, 'conv_win',cfg.calib.conv_win, ...
             'WP',cfg.WP, 's_source',cfg.s_source, 'fld',L.fld);
tag    = 'single';  if USE_BIAS, tag = 'eighteen'; end
matdir = fullfile(CAL, 'data', MODEL, '.mat');
if ~exist(matdir, 'dir'), mkdir(matdir); end
matfile = fullfile(matdir, sprintf('calib_current_%s_N%d_R%03d_%s.mat', GEOM, rec.npts, round(R*1e6), tag));
save(matfile, '-struct', 'rec');
fprintf('saved %s\n', matfile);
fprintf('[current] l_hat = %.2f um  gI_hat = %.4f mT/A  e = (%.4f, %.4f, %.4f)  NMAE = %.4f %%  RMS = %.4g mT\n', ...
        rec.l_hat*1e6, rec.gI_hat, rec.e, rec.NMAE, sqrt(rec.J/(3*rec.npts)));

%% ---- PDF -------------------------------------------------------------------
sp_emit(matfile);

%% ---- local -----------------------------------------------------------------
function [P, b] = samples(n, R, ctr, L)
    xk = -(1:n).' * (R/n);
    P  = [xk, zeros(n,1), zeros(n,1)] + ctr;
    P  = P(L.air_stencil(P), :);
    b  = L.interp(P);
end
