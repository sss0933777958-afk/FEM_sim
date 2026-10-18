% calib_noise.m -- noise effect on the flux calibration parameter l_hat (NTU single / pair pole).
% Follows reference/Noise analyze/calib_noise_v2.pdf (sections 1-6), same symbols:
%   W = sum_i S_i^T S_i,  c_j = sum_i S_i^T b_ij,  g_j = W^-1 c_j,  eps_ij = b_ij - S_i g_j
%   J(l_hat) = sum_j sum_i eps_ij^T eps_ij
%   J' = -2 sum_j sum_i eps_ij^T S'_i g_j
%   g'_j = W^-1 sum_i (S'_i^T eps_ij - S_i^T S'_i g_j)
%   eps'_ij = -S'_i g_j - S_i g'_j
%   J'' = 2 sum_j sum_i ||eps'_ij||^2 - 2 sum_j sum_i eps_ij^T S''_i g_j - 4 sum_j sum_i eps_ij^T S'_i g'_j
%   noise (sec. 1): b~_ij = b_ij + n_ij,  delta_l_hat ~= -J~'(l_hat_0) / J~''(l_hat_0)
% CASE = 'single_pole' | 'pair_pole' (P1 excitation): 1-D model b_i = a I / (x_i - l_hat)^2 on the pole
%   axis, i.e. S_i = s_i = I / (x_i - l_hat)^2 (scalar), g = a, N_I = 1; data |b_x| at the 17 points
%   x = -20 ... -340 um in front of the tip (origin = tip, l_hat > 0 into the pole).
% s'_i and s''_i: s_i is differentiated symbolically with respect to l_hat and evaluated at every point.
% Lengths in um, b in mT (per 1 A).
% Checks: J(l_hat_0) against the reference fit; s' against its closed form; J', J'' against finite
% differences of J; J'(l_hat_0) ~ 0.
% Noise: n_i ~ N(0, SIG2) independent; NMC realizations.  REFIT = true also refits each noisy data set
% (true minimizer) to compare with the linearized delta_l_hat.
% long2016 hexapole: utils/scripts/long2016_hexapole_halfcut/calib_noise.m.
% Result (sec. 5-6): linearized in the noise, with the noise-free eps'_ij at l_hat_0,
%   delta_l_hat ~= -(2 / J''(l_hat_0)) sum_j sum_i eps'_ij^T n_ij = EPS^T N,
%   EPS = -(2 / J'') [eps'_11; eps'_21; ...; eps'_NcNI],  N = [n_11; n_21; ...; n_NcNI]  (Nc NI x 1 here, 1-D)
%   sigma_l_hat = 2 sqrt(sum_j sum_i ||eps'_ij||^2) / J''(l_hat_0) * sigma  (E[n_ij] = 0, independent, Var = sigma^2).
%   SIG_U = sigma_l_hat / sigma [um/mT]; dl_res(m) = EPS^T N for each noise realization.
% Params: CASE, SIG2 [mT^2] (noise variance; [] = checks only, no noise part), NMC, SEED, REFIT.
% Output: utils/data/NTU_hexapole/calib_noise_<CASE>.mat
if ~exist('CASE','var'), CASE = 'single_pole'; end
if ~exist('SIG2','var'), SIG2 = []; end
if ~exist('NMC','var'), NMC = 1000; end
if ~exist('SEED','var'), SEED = 0; end
if ~exist('REFIT','var'), REFIT = true; end
clearvars -except CASE SIG2 NMC SEED REFIT;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell

%% ---- data, symbolic element function, noiseless minimizer l_hat_0 ----
PRE = struct('single_pole','sp', 'pair_pole','pp');  PRE = PRE.(CASE);
MODEL = 'NTU_hexapole';
DAT = fullfile(CALROOT, 'utils', 'data', MODEL);
A = load(fullfile(DAT, [PRE '_axis.mat']));
P  = A.x(:) * 1e6;                                        % Np x 1 [um], tip frame (negative, in front of the tip)
B  = abs(A.b(:,1));                                       % |b_x| [mT], P1 excitation
Pc = A.I;                                                 % drive current [A] (enters s_i)
syms l x I real
s  = I / (x - l)^2;                                       % s_i
SF = {matlabFunction(s,           'Vars', {x, I, l}), ... % SF{order+1}
      matlabFunction(diff(s, l),    'Vars', {x, I, l}), ...
      matlabFunction(diff(s, l, 2), 'Vars', {x, I, l})};
l0 = fit1d(P, Pc, B, A.xc, SF);
Jref = jcost(P, Pc, B, l0, SF);
fprintf('l_hat_0 = %.4f um (saved %.4f)   J = %.6g mT^2 (saved %.6g)\n', l0, A.xc, Jref, A.J);

%% ---- checks at l_hat_0 ----
[~, S1] = build_S(P, Pc, l0, SF);
S1c = 2 * Pc ./ (P - l0).^3;                              % closed form of s'
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
            dl_fit(m) = fit1d(P, Pc, Bn, l0, SF) - l0;
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

OUTN = ['calib_noise_' CASE '.mat'];
save(fullfile(DAT, OUTN), 'CASE', 'MODEL', 'P', 'B', 'Pc', ...
     'l0', 'J0', 'J1', 'J2', 'EPS', 'SIG_U', 'SIG_L', 'SIG2', 'NMC', 'SEED', 'REFIT', 'dl_lin', 'dl_fit', 'dl_res');
fprintf('saved %s\n', fullfile(DAT, OUTN));

%% ---- local functions ----
function [S, S1, S2] = build_S(P, Pc, l, SF)
% S, S', S'' over the points: Np x 1.
S = SF{1}(P, Pc, l);  S1 = SF{2}(P, Pc, l);  S2 = SF{3}(P, Pc, l);
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

function l = fit1d(P, Pc, B, lg, SF)
% 1-D minimizer of J(l_hat) near lg: coarse scan, then fminbnd (as in sp_noise_mc.m).
lgrid = (max(max(P) + 0.5, lg - 100) : 0.5 : lg + 100).';
Jg = arrayfun(@(q) jcost(P, Pc, B, q, SF), lgrid);
[~, im] = min(Jg);
assert(im > 1 && im < numel(lgrid), 'calib_noise:edge', 'minimum on the scan edge (guess %.1f um)', lg);
l = fminbnd(@(q) jcost(P, Pc, B, q, SF), lgrid(im-1), lgrid(im+1), optimset('TolX', 1e-9));
end
