% calib_noise.m -- noise effect on the flux calibration parameter l_hat (long2016 hexapole).
% Follows reference/Noise analyze/calib_noise_v2.pdf (sections 1-6), same symbols:
%   W = sum_i S_i^T S_i,  c_j = sum_i S_i^T b_ij,  g_j = W^-1 c_j,  eps_ij = b_ij - S_i g_j
%   J(l_hat) = sum_j sum_i eps_ij^T eps_ij
%   J' = -2 sum_j sum_i eps_ij^T S'_i g_j
%   g'_j = W^-1 sum_i (S'_i^T eps_ij - S_i^T S'_i g_j)
%   eps'_ij = -S'_i g_j - S_i g'_j
%   J'' = 2 sum_j sum_i ||eps'_ij||^2 - 2 sum_j sum_i eps_ij^T S''_i g_j - 4 sum_j sum_i eps_ij^T S'_i g'_j
%   noise (sec. 1): b~_ij = b_ij + n_ij,  delta_l_hat ~= -J~'(l_hat_0) / J~''(l_hat_0)
% Hexapole, single (e = 0): S_i = [s_i1 ... s_i6] (3x6), s_ik = d_ik / ||d_ik||^3,
%   d_ik = p_bar_i - p_bar_ck, p_bar_i = p_i / l_hat; current, axes shells N43 (Nr 7), R 150 um, N_I = 6.
% s'_ik and s''_ik: each element of s_ik is differentiated symbolically with respect to l_hat, the element
% functions are evaluated at every point and assembled back into S'_i and S''_i.
% Lengths in um, b in mT (per 1 A).
% Checks: J(l_hat_0) against the reference fit; s' against its closed form; J', J'' against finite
% differences of J; J'(l_hat_0) ~ 0.
% Noise: n_ij ~ N(0, SIG2) independent per component; NMC realizations.  REFIT = true also refits each
% noisy data set (true minimizer) to compare with the linearized delta_l_hat.
% BIAS = true: eighteen calibration (axes shells N55, Nr 9, R 150 um); e_hat is FIXED at its
%   calibrated value, p_bar_ck = Pc_base + E(e_hat), and only l_hat is perturbed (same 1-D formulas).
%   l_hat_0 = minimizer of J over l_hat with e_hat fixed; refit = 1-D fminbnd over l_hat, e_hat fixed.
% NTU single / pair pole: utils/scripts/NTU_hexapole/calib_noise.m.
% Result (sec. 5-6): linearized in the noise, with the noise-free eps'_ij at l_hat_0,
%   delta_l_hat ~= -(2 / J''(l_hat_0)) sum_j sum_i eps'_ij^T n_ij = EPS^T N,
%   EPS = -(2 / J'') [eps'_11; eps'_21; ...; eps'_NcNI],  N = [n_11; n_21; ...; n_NcNI]  (3 Nc NI x 1)
%   sigma_l_hat = 2 sqrt(sum_j sum_i ||eps'_ij||^2) / J''(l_hat_0) * sigma  (E[n_ij] = 0, independent, Cov = sigma^2 E3).
%   SIG_U = sigma_l_hat / sigma [um/mT]; dl_res(m) = EPS^T N for each noise realization.
% Params: BIAS, SIG2 [mT^2] (noise variance; [] = checks only, no noise part), NMC, SEED, REFIT,
%   LADDER (Nr list; [] = off). LADDER set: single and eighteen calibration (BIAS ignored), sigma = sqrt(SIG2)
%   (default 0.0128 mT). For each Nr (N = 6 Nr + 1 points, R 150 um) and each model: calibration fit,
%   then e fixed (e = 0 / e_hat), l_hat_0 = 1-D minimizer, closed-form sigma_l_hat (sec. 6). No Monte Carlo.
%   Convergence: first window of KW = 10 consecutive step changes |s_k - s_k-1| / s_k-1 < TOLS = 0.01 %.
% Output: utils/data/long2016_hexapole_halfcut/calib_noise[_bias].mat; LADDER: calib_noise_ladder.mat
if ~exist('BIAS','var'), BIAS = false; end
if ~exist('SIG2','var'), SIG2 = []; end
if ~exist('NMC','var'), NMC = 1000; end
if ~exist('SEED','var'), SEED = 0; end
if ~exist('REFIT','var'), REFIT = true; end
if ~exist('LADDER','var'), LADDER = []; end
clearvars -except BIAS SIG2 NMC SEED REFIT LADDER;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
addpath(fullfile(CALROOT, 'function', 'hexapole'), fullfile(CALROOT, 'common_path'));

%% ---- LADDER: sigma_l_hat vs number of calibration points (eighteen) ----
if ~isempty(LADDER)
    MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';  VARIANT = 'maxwell';
    R_SEL = 150e-6;  TOLS = 0.01;  KW = 10;
    if isempty(SIG2), SIG2 = 0.0128^2; end
    cfg = model_config(MODEL, GEOM);
    opt = struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, 'frame','actuator', ...
                 'raw',extract_maxwell_data(cfg, 'all', VARIANT));
    SF = sym_S();
    % columns: 1 = single parameter (e = 0), 2 = eighteen parameters (e_hat calibrated, then fixed)
    Nr = LADDER(:);  nq = numel(Nr);
    N = nan(nq, 1);  l0 = nan(nq, 2);  J0 = nan(nq, 2);  J2 = nan(nq, 2);  SIG_U = nan(nq, 2);
    e_hat = nan(17, nq);
    lg = [0.5e-3, 0.5e-3];  t0 = tic;
    for q = 1:nq
        evalc('[Pm, B] = conv_design_ws(Nr(q), R_SEL, opt);');           % silence the iron-exclusion report
        N(q) = size(Pm, 1);
        P  = Pm * 1e6;
        for m = 1:2
            [eq, lq] = fitting(Pm, B, cfg.Pc_base, lg(m), m == 2);       % calibration
            lg(m) = lq;
            if m == 2, e_hat(:,q) = eq; end
            Pc = make_Pc(eq, cfg.Pc_base);                               % e fixed (e = 0 for single)
            l0(q,m) = fitl(P, Pc, B, lq*1e6, SF);
            [J0(q,m), ~, J2(q,m), EPS1] = jcost(P, Pc, B, l0(q,m), SF);
            SIG_U(q,m) = 2 * norm(EPS1) / J2(q,m);
        end
        fprintf('Nr=%3d N=%4d | sigma_l_hat single %.6g, eighteen %.6g um | %.0f s\n', ...
                Nr(q), N(q), SIG_U(q,:)*sqrt(SIG2), toc(t0));
    end
    SIG_L = SIG_U * sqrt(SIG2);
    ch = [nan(1, 2); abs(diff(SIG_L)) ./ SIG_L(1:end-1,:) * 100];        % step change [%]
    iw = nan(1, 2);
    for m = 1:2
        for i = 2:nq-KW+1
            if all(ch(i:i+KW-1,m) < TOLS), iw(m) = i - 1;  break; end    % window start = rung before the first step
        end
        if isnan(iw(m))
            fprintf('model %d: not converged up to N = %d (min step change %.4g %%)\n', m, N(end), min(ch(:,m)));
        else
            fprintf('model %d: converged at N = %d, sigma_l_hat %.6g um\n', m, N(iw(m)), SIG_L(iw(m),m));
        end
    end
    DAT = fullfile(CALROOT, 'utils', 'data', MODEL);
    save(fullfile(DAT, 'calib_noise_ladder.mat'), 'MODEL', 'R_SEL', 'SIG2', 'TOLS', 'KW', ...
         'Nr', 'N', 'l0', 'e_hat', 'J0', 'J2', 'SIG_U', 'SIG_L', 'ch', 'iw');
    fprintf('saved %s\n', fullfile(DAT, 'calib_noise_ladder.mat'));
    return
end

%% ---- data, symbolic element functions, noiseless minimizer l_hat_0 ----
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';  VARIANT = 'maxwell';
if BIAS
    NR = 9;  CTAG = 'axshN55_R150_eighteen';              % axes shells: N = 6*NR + 1 = 55 points
else
    NR = 7;  CTAG = 'axshN43_R150_single';                % 43 points
end
R_SEL = 150e-6;
CALIB = fullfile(CALROOT, 'data', MODEL, '.mat', ['calib_current_maxwell_' CTAG '.mat']);
ref = load(CALIB, 'l_hat', 'e', 'J');
cfg = model_config(MODEL, GEOM);
[Pm, B] = conv_design_ws(NR, R_SEL, struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, 'frame','actuator'));
P  = Pm * 1e6;                                            % Np x 3 [um], actuator frame
Pc = cfg.Pc_base;                                         % 3 x 6 charge lattice p_bar_ck (e = 0)
e_hat = zeros(17, 1);
if BIAS, e_hat = ref.e;  Pc = make_Pc(e_hat, Pc); end     % p_bar_ck = Pc_base + E(e_hat), fixed
SF = sym_S();
Np = size(P, 1);
if BIAS
    l0 = fitl(P, Pc, B, ref.l_hat*1e6, SF);               % e_hat fixed: J over l_hat only
    Jref = ref.J;
    fprintf('l_hat_0 = %.4f um (saved %.4f)   J = %.6g mT^2 (saved %.6g), e_hat fixed\n', ...
            l0, ref.l_hat*1e6, jcost(P, Pc, B, l0, SF), ref.J);
else
    [~, l0m, Jref] = fitting(Pm, B, Pc, 0.5e-3, false);
    l0 = l0m * 1e6;                                       % [um]
    fprintf('l_hat_0 = %.4f um (saved %.4f)   J = %.6g mT^2 (saved %.6g)\n', l0, ref.l_hat*1e6, Jref, ref.J);
end

%% ---- checks at l_hat_0 ----
[S, S1] = build_S(P, Pc, l0, SF);
S1c = zeros(size(S1));                                    % closed form: -1/(l ||d||^3) (E3 - 3 d d^T/||d||^2) p_bar_i
for i = 1:Np
    pb = P(i,:).' / l0;
    for k = 1:6
        dk = pb - Pc(:,k);  nd = norm(dk);
        S1c(3*i-2:3*i, k) = -1/(l0*nd^3) * (eye(3) - 3*(dk*dk.')/nd^2) * pb;
    end
end
fprintf('s'' symbolic vs closed form: max |diff| = %.3g (max |s''| = %.3g)\n', max(abs(S1(:)-S1c(:))), max(abs(S1(:))));
[J0, J1, J2, EPS1] = jcost(P, Pc, B, l0, SF);
h = 1e-4 * l0;                                            % finite-difference step [um]
[Jp, J1p] = jcost(P, Pc, B, l0 + h, SF);
[Jm, J1m] = jcost(P, Pc, B, l0 - h, SF);
fprintf('J(l_hat_0)  = %.6g mT^2 (reference fit %.6g)\n', J0, Jref);
fprintf('J''(l_hat_0) = %.4g mT^2/um   (central diff of J  %.4g)\n', J1, (Jp - Jm)/(2*h));
fprintf('J''''(l_hat_0)= %.6g mT^2/um^2 (central diff of J'' %.6g, of J %.6g)\n', ...
        J2, (J1p - J1m)/(2*h), (Jp - 2*J0 + Jm)/h^2);

%% ---- sec. 5-6: EPS and sigma_l_hat (closed form, noise-free eps' at l_hat_0) ----
EPS   = -2 / J2 * EPS1;                                 % delta_l_hat = EPS^T N
SIG_U = 2 * norm(EPS1) / J2;                            % sigma_l_hat / sigma [um/mT]
SIG_L = [];
if ~isempty(SIG2), SIG_L = SIG_U * sqrt(SIG2); end
fprintf('sigma_l_hat / sigma = 2 sqrt(sum ||eps''||^2) / J'''' = %.6g um/mT\n', SIG_U);

%% ---- sec. 1: noise, delta_l_hat ~= -J~'/J~'' ----
dl_lin = [];  dl_fit = [];  dl_res = [];
if isempty(SIG2)
    fprintf('SIG2 not set: noise part skipped.\n');
else
    rng(SEED);
    dl_lin = zeros(NMC, 1);  dl_fit = nan(NMC, 1);  dl_res = zeros(NMC, 1);
    for m = 1:NMC
        nij = sqrt(SIG2) * randn(size(B));                % n_ij
        Bn  = B + nij;                                    % b~_ij = b_ij + n_ij
        [~, Jt1, Jt2] = jcost(P, Pc, Bn, l0, SF);         % J~', J~'' at l_hat_0
        dl_lin(m) = -Jt1 / Jt2;
        N = nij(:);                                       % N = [n_11; n_21; ...]
        dl_res(m) = EPS.' * N;                            % EPS^T N
        if REFIT
            if BIAS
                lf = fitl(P, Pc, Bn, l0, SF);
            else
                [~, lf] = fitting(Pm, Bn, Pc, l0m, false);  lf = lf*1e6;
            end
            dl_fit(m) = lf - l0;
        end
    end
    fprintf('SIG2 = %g mT^2, NMC = %d\n', SIG2, NMC);
    fprintf('  result    : sigma_l_hat %.4g um, var %.4g um^2 (closed form) | EPS^T N: std %.4g um, max |EPS^T N - lin| %.3g um\n', ...
            SIG_L, SIG_L^2, std(dl_res), max(abs(dl_res - dl_lin)));
    fprintf('  linearized: mean %.4g um, std %.4g um, var %.4g um^2\n', mean(dl_lin), std(dl_lin), var(dl_lin));
    if REFIT
        fprintf('  refit     : mean %.4g um, std %.4g um, var %.4g um^2, max |lin - refit| %.3g um\n', ...
                mean(dl_fit), std(dl_fit), var(dl_fit), max(abs(dl_lin - dl_fit)));
    end
end

DAT = fullfile(CALROOT, 'utils', 'data', MODEL);
OUTN = 'calib_noise.mat';
if BIAS, OUTN = 'calib_noise_bias.mat'; end
save(fullfile(DAT, OUTN), 'BIAS', 'MODEL', 'P', 'B', 'Pc', 'e_hat', ...
     'l0', 'J0', 'J1', 'J2', 'EPS', 'SIG_U', 'SIG_L', 'SIG2', 'NMC', 'SEED', 'REFIT', 'dl_lin', 'dl_fit', 'dl_res');
fprintf('saved %s\n', fullfile(DAT, OUTN));

%% ---- local functions ----
function SF = sym_S()
% Element functions of s_ik, s'_ik, s''_ik (symbolic in l_hat): SF{order+1, component}.
syms l x y z cx cy cz real
d  = [x; y; z] / l - [cx; cy; cz];
s  = d / (d.' * d)^(3/2);                                 % s_ik (3x1)
s1 = diff(s, l);                                          % s'
s2 = diff(s, l, 2);                                       % s''
VARS = {x, y, z, cx, cy, cz, l};
SF = cell(3, 3);
for a = 1:3
    SF{1,a} = matlabFunction(s(a),  'Vars', VARS);
    SF{2,a} = matlabFunction(s1(a), 'Vars', VARS);
    SF{3,a} = matlabFunction(s2(a), 'Vars', VARS);
end
end

function [S, S1, S2] = build_S(P, Pc, l, SF)
% S, S', S'' stacked over points: 3Np x 6, rows [x; y; z] per point (same layout as fitting.m).
Np = size(P, 1);
S = zeros(3*Np, 6);  S1 = S;  S2 = S;
for k = 1:6
    for a = 1:3
        rows = a:3:3*Np;
        args = {P(:,1), P(:,2), P(:,3), Pc(1,k), Pc(2,k), Pc(3,k), l};
        S(rows, k)  = SF{1,a}(args{:});
        S1(rows, k) = SF{2,a}(args{:});
        S2(rows, k) = SF{3,a}(args{:});
    end
end
end

function [J, J1, J2, EPS1] = jcost(P, Pc, B, l, SF)
% J, J', J'' (sec. 1, 2, 4); sums over i are inner products of the stacked vectors.
% EPS1 = eps'_ij stacked over i, then j (same order as B(:)).
[S, S1, S2] = build_S(P, Pc, l, SF);
W = S.' * S;
J = 0;  J1 = 0;  J2 = 0;  EPS1 = zeros(numel(B), 1);  nr = size(B, 1);
for j = 1:size(B, 2)
    b    = B(:, j);
    g    = W \ (S.' * b);                                 % g_j = W^-1 c_j
    eps0 = b - S * g;                                     % eps_ij
    g1   = W \ (S1.' * eps0 - S.' * S1 * g);              % g'_j
    eps1 = -S1 * g - S * g1;                              % eps'_ij
    EPS1((j-1)*nr+1 : j*nr) = eps1;
    J  = J  + eps0.' * eps0;
    J1 = J1 - 2 * eps0.' * S1 * g;
    J2 = J2 + 2 * (eps1.' * eps1) - 2 * eps0.' * (S2 * g) - 4 * eps0.' * (S1 * g1);
end
end

function l = fitl(P, Pc, B, lg, SF)
% e_hat fixed: minimizer of J over l_hat in [lg - 50, lg + 50] um (fminbnd).
l = fminbnd(@(q) jcost(P, Pc, B, q, SF), lg - 50, lg + 50, optimset('TolX', 1e-9));
assert(abs(l - lg) < 49, 'calib_noise:edge', 'minimum at the search edge (guess %.1f um)', lg);
end

function Pc = make_Pc(e17, Pc_base)
% Pc_bar = Pc_base + E(e), with the e6z constraint (same as fitting.m).
E = zeros(3, 6);
E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
E(:,5) = e17(13:15);
E(1,6) = e17(16);      E(2,6) = e17(17);
E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
Pc = Pc_base + E;
end
