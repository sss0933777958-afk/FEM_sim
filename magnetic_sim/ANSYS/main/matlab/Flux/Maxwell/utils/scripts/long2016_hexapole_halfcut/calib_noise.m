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
% Document (reference/Noise analyze/VarianceCalculation.pdf, its symbols in doc_var / doc_dl): Var_u = Var(delta_l_hat) / sigma_n^2,
%   Var_dl = Var(delta_l_hat) at sigma_n^2 = SIG2, dl_doc = its delta_l_hat - E[delta_l_hat] per realization;
%   *_add = the same with the omitted delta_l_hat (S'_o,i)^T n_ij term put back (document's denominator kept).
%   *_v2 = report/406/10.21/10.pptx, Version 2 (second-order expansion, without omission), its final sigma_l_hat.
% Bias (draft, sec. 4 "Bias and Random Fluctuation of delta_l_hat", see doc_bias): E_u = E[delta_l_hat] / sigma_n^2,
%   E_dl = E[delta_l_hat] at sigma_n^2 = SIG2, dl_doc_b = N_1/(2A) + (N_2 + N_1 B)/(2A) per realization;
%   *_add, *_v2: the same steps on the cost of VarianceCalculation..pdf sec. 2 / 3.
% Params: BIAS, SIG2 [mT^2] (noise variance; [] = checks only, no noise part), NMC, SEED, REFIT,
%   REUSE (true = keep the saved refits, regenerate the same noise, add the document's results),
%   NPAR (0 = serial loop; > 0 = split over NPAR separate MATLAB processes, chunks of 1e4 realizations, chunk c
%   drawn from substream c of Threefry(SEED); output calib_noise[_bias]_<NMC, e.g. 4M>.mat; no REUSE),
%   PART (NPAR > 0: 1..NPAR = compute chunks PART, PART + NPAR, ... -> <output>_p<PART>.mat and stop;
%   0 = merge the part files and continue as the serial run),
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
if ~exist('REUSE','var'), REUSE = false; end
if ~exist('NPAR','var'), NPAR = 0; end
if ~exist('PART','var'), PART = 0; end
clearvars -except BIAS SIG2 NMC SEED REFIT LADDER REUSE NPAR PART;  close all;
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
    % REUSE = true: Nr, l_hat_0 and e_hat from the saved calib_noise_ladder.mat (no calibration fits).
    if REUSE
        old = load(fullfile(CALROOT, 'utils', 'data', MODEL, 'calib_noise_ladder.mat'), 'Nr', 'l0', 'e_hat', 'SIG2');
        LADDER = old.Nr;  SIG2 = old.SIG2;
    end
    Nr = LADDER(:);  nq = numel(Nr);
    N = nan(nq, 1);  l0 = nan(nq, 2);  J0 = nan(nq, 2);  J2 = nan(nq, 2);  SIG_U = nan(nq, 2);  Var_u = nan(nq, 2);
    Var_u_add = nan(nq, 2);  Var_u_v2 = nan(nq, 2);
    e_hat = nan(17, nq);
    lg = [0.5e-3, 0.5e-3];  t0 = tic;
    for q = 1:nq
        evalc('[Pm, B] = conv_design_ws(Nr(q), R_SEL, opt);');           % silence the iron-exclusion report
        N(q) = size(Pm, 1);
        P  = Pm * 1e6;
        for m = 1:2
            if REUSE
                eq = zeros(17, 1);  if m == 2, eq = old.e_hat(:,q); end
                Pc = make_Pc(eq, cfg.Pc_base);
                l0(q,m) = old.l0(q,m);
            else
                [eq, lq] = fitting(Pm, B, cfg.Pc_base, lg(m), m == 2);   % calibration
                lg(m) = lq;
                Pc = make_Pc(eq, cfg.Pc_base);                           % e fixed (e = 0 for single)
                l0(q,m) = fitl(P, Pc, B, lq*1e6, SF);
            end
            if m == 2, e_hat(:,q) = eq; end
            [J0(q,m), ~, J2(q,m), EPS1] = jcost(P, Pc, B, l0(q,m), SF);
            SIG_U(q,m) = 2 * norm(EPS1) / J2(q,m);
            [S, S1, S2] = build_S(P, Pc, l0(q,m), SF);
            Var_u(q,m) = doc_var(S, S1, B);                              % VarianceCalculation.pdf, sec. 1
            Var_u_add(q,m) = doc_var(S, S1, B, 'add');                   % sec. 2
            Var_u_v2(q,m) = doc_var(S, S1, B, 'v2', S2);                 % sec. 3
        end
        fprintf('Nr=%3d N=%4d | eta single %.6g (doc %.6g / %.6g / %.6g), eighteen %.6g (doc %.6g) um/mT | %.0f s\n', ...
                Nr(q), N(q), SIG_U(q,1), sqrt(Var_u(q,1)), sqrt(Var_u_add(q,1)), sqrt(Var_u_v2(q,1)), ...
                SIG_U(q,2), sqrt(Var_u(q,2)), toc(t0));
    end
    SIG_L = SIG_U * sqrt(SIG2);  Var_dl = Var_u * SIG2;  Var_dl_add = Var_u_add * SIG2;  Var_dl_v2 = Var_u_v2 * SIG2;
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
         'Nr', 'N', 'l0', 'e_hat', 'J0', 'J2', 'SIG_U', 'SIG_L', 'ch', 'iw', 'Var_u', 'Var_dl', ...
         'Var_u_add', 'Var_dl_add', 'Var_u_v2', 'Var_dl_v2');
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
OUTN = 'calib_noise.mat';
if BIAS, OUTN = 'calib_noise_bias.mat'; end
if NPAR > 0                                               % parallel run: own file, the serial results stay
    if NMC >= 1e6, NTAG = sprintf('%gM', NMC/1e6); else, NTAG = sprintf('%gk', NMC/1e3); end
    OUTN = strrep(OUTN, '.mat', ['_' NTAG '.mat']);
end
CALIB =fullfile(CALROOT, 'data', MODEL, '.mat', ['calib_current_maxwell_' CTAG '.mat']);
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

%% ---- reference/Noise analyze/VarianceCalculation.pdf: Var(delta_l_hat), the document's symbols ----
% S_o,i = S(p_i / l_hat_o),  S'_o,i = dS(p_i / l_hat)/dl_hat at l_hat_o  (S, S1 above),  sigma_n^2 = SIG2.
% Var_u = Var(delta_l_hat) / sigma_n^2 [um^2/mT^2];  Var_dl = Var(delta_l_hat) [um^2]  (see doc_var).
% *_add: the omitted delta_l_hat (S'_o,i)^T n_ij term put back (doc_var / doc_dl with VER = 'add').
% *_v2: report/406/10.21/10.pptx Version 2 (VER = 'v2', needs S''_o,i).
Var_u  = doc_var(S, S1, B);
Var_u_add = doc_var(S, S1, B, 'add');
[~, ~, S2] = build_S(P, Pc, l0, SF);                     % S''_o,i
Var_u_v2 = doc_var(S, S1, B, 'v2', S2);                  % 10.pptx, Version 2
Var_dl = [];  Var_dl_add = [];  Var_dl_v2 = [];
if ~isempty(SIG2), Var_dl = Var_u * SIG2;  Var_dl_add = Var_u_add * SIG2;  Var_dl_v2 = Var_u_v2 * SIG2; end
E_u = doc_bias(S, S1, B);                                % draft sec. 4: E[delta_l_hat] / sigma_n^2 (sec. 1 cost)
E_u_add = doc_bias(S, S1, B, [], 'add');                 % sec. 2 cost
E_u_v2  = doc_bias(S, S1, B, [], 'v2', S2);              % sec. 3 cost
E_dl = [];  E_dl_add = [];  E_dl_v2 = [];
if ~isempty(SIG2), E_dl = E_u * SIG2;  E_dl_add = E_u_add * SIG2;  E_dl_v2 = E_u_v2 * SIG2; end
fprintf('E[delta_l_hat] / sigma_n^2 = %.6g (sec. 1), %.6g (sec. 2), %.6g (sec. 3) um/mT^2\n', E_u, E_u_add, E_u_v2);
fprintf('document: sqrt(Var(delta_l_hat)) / sigma_n = %.6g um/mT, omitted term added %.6g, Version 2 %.6g (this script %.6g)\n', ...
        sqrt(Var_u), sqrt(Var_u_add), sqrt(Var_u_v2), SIG_U);

%% ---- sec. 1: noise, delta_l_hat ~= -J~'/J~'' ----
% REUSE = true: keep the saved refits (dl_fit) and regenerate the same noise (same SEED, same draw order)
%   to add the document's delta_l_hat (dl_doc); the regenerated EPS^T N is checked against the saved dl_res.
dl_lin = [];  dl_fit = [];  dl_res = [];  dl_doc = [];  dl_doc_add = [];  dl_doc_v2 = [];  dl_doc_b = [];  dl_doc_b_add = [];  dl_doc_b_v2 = [];
if REUSE
    old = load(fullfile(CALROOT, 'utils', 'data', MODEL, OUTN), 'SIG2', 'NMC', 'SEED', 'REFIT', 'dl_fit', 'dl_res');
    SIG2 = old.SIG2;  NMC = old.NMC;  SEED = old.SEED;  REFIT = old.REFIT;
    SIG_L = SIG_U * sqrt(SIG2);  Var_dl = Var_u * SIG2;  Var_dl_add = Var_u_add * SIG2;  Var_dl_v2 = Var_u_v2 * SIG2;
    E_dl = E_u * SIG2;  E_dl_add = E_u_add * SIG2;  E_dl_v2 = E_u_v2 * SIG2;
end
if isempty(SIG2)
    fprintf('SIG2 not set: noise part skipped.\n');
else
    if NPAR > 0
    assert(~REUSE, 'calib_noise:par', 'REUSE needs the serial loop (NPAR = 0)');
    NK = 1e4;  NCH = ceil(NMC / NK);  NMC = NCH * NK;            % chunks of NK realizations
    PFN = @(k) fullfile(CALROOT, 'utils', 'data', MODEL, strrep(OUTN, '.mat', sprintf('_p%d.mat', k)));
    if PART > 0                                           % process PART of NPAR: chunks PART, PART + NPAR, ...
        Pm_ = []; l0m_ = [];  if ~BIAS, Pm_ = Pm;  l0m_ = l0m; end
        CH = PART:NPAR:NCH;  RP = zeros(NK, 10, numel(CH));  t0 = tic;
        for q = 1:numel(CH)
            rs = RandStream('threefry4x64_20', 'Seed', SEED);  rs.Substream = CH(q);
            RP(:,:,q) = mc_chunk(rs, NK, SIG2, P, Pc, B, l0, SF, EPS, S, S1, S2, BIAS, REFIT, Pm_, l0m_);
            fprintf('part %d: chunk %d (%d / %d), %.0f s\n', PART, CH(q), q, numel(CH), toc(t0));
        end
        save(PFN(PART), 'CH', 'RP', 'NMC', 'SEED', 'SIG2');
        fprintf('saved %s\n', PFN(PART));
        return
    end
    R = zeros(NK, 10, NCH);                               % PART = 0: merge the NPAR part files
    for k = 1:NPAR
        d = load(PFN(k), 'CH', 'RP', 'NMC', 'SEED', 'SIG2');
        assert(d.NMC == NMC && d.SEED == SEED && d.SIG2 == SIG2, 'calib_noise:part', 'part %d settings differ', k);
        R(:,:,d.CH) = d.RP;
    end
    assert(all(any(R(:,1,:) ~= 0, 1), 'all'), 'calib_noise:part', 'missing chunks');
    R = reshape(permute(R, [1 3 2]), NMC, 10);
    dl_lin = R(:,1);  dl_res = R(:,2);  dl_doc = R(:,3);  dl_doc_add = R(:,4);  dl_doc_v2 = R(:,5);
    dl_doc_b = R(:,6);  dl_doc_b_add = R(:,7);  dl_doc_b_v2 = R(:,8);  dl_fit = R(:,9);
    clear R
    else
    rng(SEED);
    dl_lin = zeros(NMC, 1);  dl_fit = nan(NMC, 1);  dl_res = zeros(NMC, 1);  dl_doc = zeros(NMC, 1);  dl_doc_add = zeros(NMC, 1);  dl_doc_v2 = zeros(NMC, 1);  dl_doc_b = zeros(NMC, 1);  dl_doc_b_add = zeros(NMC, 1);  dl_doc_b_v2 = zeros(NMC, 1);
    for m = 1:NMC
        nij = sqrt(SIG2) * randn(size(B));                % n_ij
        Bn  = B + nij;                                    % b~_ij = b_ij + n_ij
        [~, Jt1, Jt2] = jcost(P, Pc, Bn, l0, SF);         % J~', J~'' at l_hat_0
        dl_lin(m) = -Jt1 / Jt2;
        N = nij(:);                                       % N = [n_11; n_21; ...]
        dl_res(m) = EPS.' * N;                            % EPS^T N
        dl_doc(m) = doc_dl(S, S1, B, nij);                % document: delta_l_hat - E[delta_l_hat]
        dl_doc_add(m) = doc_dl(S, S1, B, nij, 'add');     % document + omitted (S'_o,i)^T n_ij term
        dl_doc_v2(m) = doc_dl(S, S1, B, nij, 'v2', S2);   % 10.pptx, Version 2
        [~, dl_doc_b(m)] = doc_bias(S, S1, B, nij);       % draft sec. 4: N_1/(2A) + (N_2 + N_1 B)/(2A)
        [~, dl_doc_b_add(m)] = doc_bias(S, S1, B, nij, 'add');
        [~, dl_doc_b_v2(m)]  = doc_bias(S, S1, B, nij, 'v2', S2);
        if REUSE
            dl_fit = old.dl_fit;
        elseif REFIT
            if BIAS
                lf = fitl(P, Pc, Bn, l0, SF);
            else
                [~, lf] = fitting(Pm, Bn, Pc, l0m, false);  lf = lf*1e6;
            end
            dl_fit(m) = lf - l0;
        end
    end
    end
    fprintf('SIG2 = %g mT^2, NMC = %d\n', SIG2, NMC);
    fprintf('  result    : sigma_l_hat %.4g um, var %.4g um^2 (closed form) | EPS^T N: std %.4g um, max |EPS^T N - lin| %.3g um\n', ...
            SIG_L, SIG_L^2, std(dl_res), max(abs(dl_res - dl_lin)));
    fprintf('  linearized: mean %.4g um, std %.4g um, var %.4g um^2\n', mean(dl_lin), std(dl_lin), var(dl_lin));
    fprintf('  document  : Var(delta_l_hat) %.4g um^2, sqrt %.4g um (closed form) | per realization: mean %.4g um, std %.4g um\n', ...
            Var_dl, sqrt(Var_dl), mean(dl_doc), std(dl_doc));
    fprintf('  doc + term: Var(delta_l_hat) %.4g um^2, sqrt %.4g um (closed form) | per realization: mean %.4g um, std %.4g um\n', ...
            Var_dl_add, sqrt(Var_dl_add), mean(dl_doc_add), std(dl_doc_add));
    fprintf('  Version 2 : Var(delta_l_hat) %.4g um^2, sigma_l_hat %.4g um (closed form) | per realization: mean %.4g um, std %.4g um\n', ...
            Var_dl_v2, sqrt(Var_dl_v2), mean(dl_doc_v2), std(dl_doc_v2));
    fprintf('  bias      : E[delta_l_hat] %.4g um (closed form) | per realization (N_1 + N_2 + N_1 B)/(2A): mean %.4g um\n', ...
            E_dl, mean(dl_doc_b));
    fprintf('  bias sec.2: %.4g um (closed form) | per realization mean %.4g um;  sec.3: %.4g um | %.4g um\n', ...
            E_dl_add, mean(dl_doc_b_add), E_dl_v2, mean(dl_doc_b_v2));
    if REUSE, fprintf('  REUSE: max |EPS^T N - saved| %.3g um (same noise)\n', max(abs(dl_res - old.dl_res))); end
    if REFIT
        fprintf('  refit     : mean %.4g um, std %.4g um, var %.4g um^2, max |lin - refit| %.3g um\n', ...
                mean(dl_fit), std(dl_fit), var(dl_fit), max(abs(dl_lin - dl_fit)));
    end
end

DAT = fullfile(CALROOT, 'utils', 'data', MODEL);
save(fullfile(DAT, OUTN), 'BIAS', 'MODEL', 'P', 'B', 'Pc', 'e_hat', ...
     'l0', 'J0', 'J1', 'J2', 'EPS', 'SIG_U', 'SIG_L', 'SIG2', 'NMC', 'SEED', 'REFIT', 'dl_lin', 'dl_fit', 'dl_res', ...
     'Var_u', 'Var_dl', 'dl_doc', 'Var_u_add', 'Var_dl_add', 'dl_doc_add', 'Var_u_v2', 'Var_dl_v2', 'dl_doc_v2', ...
     'E_u', 'E_dl', 'dl_doc_b', 'E_u_add', 'E_dl_add', 'dl_doc_b_add', 'E_u_v2', 'E_dl_v2', 'dl_doc_b_v2');
fprintf('saved %s\n', fullfile(DAT, OUTN));
if NPAR > 0 && REFIT
    % simulated statistics only (the full file exceeds the GitHub size limit): <output>_stat.mat
    E_sim  = mean(dl_fit);                                % simulated bias = mean(l_hat) - l_hat_0 [um]
    SE_sim = std(dl_fit) / sqrt(NMC);                     % its standard error [um]
    std_fit = std(dl_fit);                                % std of the refitted l_hat [um]
    save(fullfile(DAT, strrep(OUTN, '.mat', '_stat.mat')), 'MODEL', 'l0', 'SIG2', 'NMC', 'SEED', ...
         'E_sim', 'SE_sim', 'std_fit', 'SIG_L', 'Var_dl', 'Var_dl_add', 'Var_dl_v2', 'E_dl', 'E_dl_add', 'E_dl_v2');
    fprintf('saved %s\n', fullfile(DAT, strrep(OUTN, '.mat', '_stat.mat')));
end

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

function Var_u = doc_var(S_o, Sp_o, b, VER, Spp_o)
% VarianceCalculation.pdf, last equation: Var(delta_l_hat) / sigma_n^2 [um^2/mT^2].
% S_o, Sp_o: S_o,i and S'_o,i stacked over i (3Nc x 6); b: b_ij stacked over i, column j (3Nc x N_I).
% ADD = true: the omitted term delta_l_hat (S'_o,i)^T n_ij is kept in the expansion (s'_nj = sum_i (S'_o,i)^T n_ij),
%   which adds - s_oj^T W_o^-1 s'_nj to the numerator of delta_l_hat; with E[s_nj s_nj^T] = sigma_n^2 W_o,
%   E[s'_nj s'_nj^T] = sigma_n^2 sum_i (S'_o,i)^T S'_o,i, E[s_nj s'_nj^T] = sigma_n^2 sum_i S_o,i^T S'_o,i:
%   Var(delta_l_hat) ~= sum_j (a_j^T W_o a_j - 2 a_j^T (sum_i S_o,i^T S'_o,i) c_j + c_j^T (sum_i (S'_o,i)^T S'_o,i) c_j)
%                       / (sum_j [...])^2 * sigma_n^2,   a_j = D_o s_oj - W_o^-1 s'_oj,  c_j = W_o^-1 s_oj.
%   The denominator is the document's.
% VER: 'v1' (default) = the document | 'add' = the omitted term put back | 'v2' = report/406/10.21/10.pptx,
%   Version 2 (second-order expansion, without omission), last slide, Spp_o = S''_o,i stacked like S_o:
%   sigma_l_hat = sqrt(sum_j (s_oj^T D_o W_o D_o s_oj - 2 s'_oj^T D_o s_oj + s'_oj^T W_o^-1 s'_oj
%                  - 2 (s_oj^T D_o - s'_oj^T W_o^-1) V_o W_o^-1 s_oj + s_oj^T W_o^-1 U_o W_o^-1 s_oj))
%                 / sum_j (s'_oj^T W_o^-1 s'_oj - 2 s'_oj^T D_o s_oj + s''_oj^T W_o^-1 s_oj + 1/2 s_oj^T Z_o s_oj) * sigma_n
%   V_o = sum_i S_o,i^T S'_o,i,  U_o = sum_i (S'_o,i)^T S'_o,i,  Z_o = 2 W_o^-1 W'_o W_o^-1 W'_o W_o^-1 - W_o^-1 W''_o W_o^-1,
%   W''_o = sum_i [(S''_o,i)^T S_o,i + 2 (S'_o,i)^T S'_o,i + S_o,i^T S''_o,i].
if nargin < 4, VER = 'v1'; end
ADD = any(strcmp(VER, {'add', 'v2'}));
W_o     = S_o.' * S_o;                                    % W_o = sum_i S_o,i^T S_o,i
W_o_inv = inv(W_o);                                       % W_o^-1
Wp_o    = Sp_o.' * S_o + S_o.' * Sp_o;                    % W'_o = sum_i [(S'_o,i)^T S_o,i + S_o,i^T S'_o,i]
D_o     = W_o_inv * Wp_o * W_o_inv;                       % D_o = W_o^-1 W'_o W_o^-1
STSp_o  = S_o.' * Sp_o;                                   % sum_i S_o,i^T S'_o,i
SpTSp_o = Sp_o.' * Sp_o;                                  % sum_i (S'_o,i)^T S'_o,i
if strcmp(VER, 'v2')
    V_o   = STSp_o;                                       % V_o = sum_i S_o,i^T S'_o,i
    U_o   = SpTSp_o;                                      % U_o = sum_i (S'_o,i)^T S'_o,i
    Wpp_o = Spp_o.' * S_o + 2 * (Sp_o.' * Sp_o) + S_o.' * Spp_o;     % W''_o
    Z_o   = 2 * W_o_inv * Wp_o * W_o_inv * Wp_o * W_o_inv - W_o_inv * Wpp_o * W_o_inv;
end
num = 0;  den = 0;
for j = 1:size(b, 2)
    s_oj  = S_o.'  * b(:, j);                             % s_oj  = sum_i S_o,i^T b_ij
    sp_oj = Sp_o.' * b(:, j);                             % s'_oj = sum_i (S'_o,i)^T b_ij
    num = num + (s_oj.' * D_o * W_o * D_o * s_oj - 2 * sp_oj.' * D_o * s_oj + sp_oj.' * W_o_inv * sp_oj);
    if strcmp(VER, 'v2')
        spp_oj = Spp_o.' * b(:, j);                       % s''_oj = sum_i (S''_o,i)^T b_ij
        num = num + (-2 * (s_oj.' * D_o - sp_oj.' * W_o_inv) * V_o * W_o_inv * s_oj + s_oj.' * W_o_inv * U_o * W_o_inv * s_oj);
        den = den + (spp_oj.' * W_o_inv * s_oj + 0.5 * s_oj.' * Z_o * s_oj);
    elseif ADD
        a_j = D_o * s_oj - W_o_inv * sp_oj;
        c_j = W_o_inv * s_oj;
        num = num - 2 * a_j.' * STSp_o * c_j + c_j.' * SpTSp_o * c_j;
    end
    den = den + (sp_oj.' * W_o_inv * sp_oj - 2 * sp_oj.' * D_o * s_oj);
end
Var_u = num / den^2;
end

function dl = doc_dl(S_o, Sp_o, b, n, VER, Spp_o)
% VarianceCalculation.pdf: delta_l_hat - E[delta_l_hat] for one noise realization n (n_ij stacked like b) [um].
% ADD = true: numerator also gets - s_oj^T W_o^-1 s'_nj (the omitted delta_l_hat (S'_o,i)^T n_ij term).
% VER = 'v2': 10.pptx Version 2, delta_l_hat - E[delta_l_hat] = sum_j (s_oj^T D_o s_nj - s'_oj^T W_o^-1 s_nj
%   - s'_nj^T W_o^-1 s_oj) / sum_j (s'_oj^T W_o^-1 s'_oj - 2 s'_oj^T D_o s_oj + s''_oj^T W_o^-1 s_oj + 1/2 s_oj^T Z_o s_oj).
if nargin < 5, VER = 'v1'; end
ADD = any(strcmp(VER, {'add', 'v2'}));
W_o     = S_o.' * S_o;
W_o_inv = inv(W_o);
Wp_o    = Sp_o.' * S_o + S_o.' * Sp_o;
D_o     = W_o_inv * Wp_o * W_o_inv;
if strcmp(VER, 'v2')
    Wpp_o = Spp_o.' * S_o + 2 * (Sp_o.' * Sp_o) + S_o.' * Spp_o;
    Z_o   = 2 * W_o_inv * Wp_o * W_o_inv * Wp_o * W_o_inv - W_o_inv * Wpp_o * W_o_inv;
end
num = 0;  den = 0;
for j = 1:size(b, 2)
    s_oj  = S_o.'  * b(:, j);
    sp_oj = Sp_o.' * b(:, j);
    s_nj  = S_o.'  * n(:, j);                             % s_nj  = sum_i S_o,i^T n_ij
    num = num + (s_oj.' * D_o - sp_oj.' * W_o_inv) * s_nj;
    if ADD
        sp_nj = Sp_o.' * n(:, j);                         % s'_nj = sum_i (S'_o,i)^T n_ij
        num = num - s_oj.' * W_o_inv * sp_nj;
    end
    den = den + (sp_oj.' * W_o_inv * sp_oj - 2 * sp_oj.' * D_o * s_oj);
    if strcmp(VER, 'v2')
        spp_oj = Spp_o.' * b(:, j);
        den = den + (spp_oj.' * W_o_inv * s_oj + 0.5 * s_oj.' * Z_o * s_oj);
    end
end
dl = num / den;
end

function [E_u, dl] = doc_bias(S_o, Sp_o, b, n, VER, Spp_o)
% Draft sec. 4, Bias and Random Fluctuation of delta_l_hat:
%   N_1 = 2 sum_j (s_oj^T D_o - (s'_oj)^T W_o^-1) s_nj,   N_2 = sum_j s_nj^T D_o s_nj,
%   A = sum_j [(s'_oj)^T W_o^-1 s'_oj - 2 (s'_oj)^T D_o s_oj],   B = 2/A sum_j (s'_oj)^T D_o s_nj,
%   delta_l_hat ~= N_1/(2A) + (N_2 + N_1 B)/(2A),   E[delta_l_hat] ~= E[N_2 + N_1 B]/(2A).
% With E[s_nj s_nk^T] = sigma_n^2 W_o (k = j), 0 (k ~= j):
%   E[N_2]   = sigma_n^2 sum_j tr(D_o W_o)
%   E[N_1 B] = 4/A sigma_n^2 sum_j (s_oj^T D_o - (s'_oj)^T W_o^-1) W_o D_o s'_oj
%   E_u = E[delta_l_hat] / sigma_n^2 [um/mT^2].
% n given: dl = N_1/(2A) + (N_2 + N_1 B)/(2A) for that realization (n_ij stacked like b) [um].
% VER = 'add' / 'v2' (VarianceCalculation..pdf sec. 2 / 3): the same steps on that section's cost, i.e.
%   s_oj + s_nj -> a_j,  s'_oj -> c_j = s'_oj + s'_nj,  (sec. 3) e_j = s''_oj + s''_nj,  W_o^-1 - dl D_o (+ dl^2/2 Z_o):
%   dJ/d(dl) = 0 to first order in dl gives dl ~= Num / Den,
%   Num = sum_j (a_j^T D_o a_j - 2 c_j^T W_o^-1 a_j)                     (its noise-free part is -J'(l_hat_o)/2 = 0)
%   Den = sum_j (2 c_j^T W_o^-1 c_j - 4 c_j^T D_o a_j [+ 2 e_j^T W_o^-1 a_j + a_j^T Z_o a_j])
%   Num = g^T n + n^T Q n,  Den = Den_o + h^T n + ...,  and to second order in the noise
%   E[dl] ~= sigma_n^2 (tr(Q) / Den_o - g^T h / Den_o^2);  VER = 'v1' reproduces the draft's N_1, N_2, A, B
%   (g^T n = N_1, n^T Q n = N_2, Den_o = 2A, h^T n = -2AB).
if nargin < 4, n = []; end
if nargin < 5, VER = 'v1'; end
ADD = any(strcmp(VER, {'add', 'v2'}));  V2 = strcmp(VER, 'v2');
W_o     = S_o.' * S_o;
W_o_inv = inv(W_o);
Wp_o    = Sp_o.' * S_o + S_o.' * Sp_o;
D_o     = W_o_inv * Wp_o * W_o_inv;
if V2
    Wpp_o = Spp_o.' * S_o + 2 * (Sp_o.' * Sp_o) + S_o.' * Spp_o;
    Z_o   = 2 * W_o_inv * Wp_o * W_o_inv * Wp_o * W_o_inv - W_o_inv * Wpp_o * W_o_inv;
end
Den_o = 0;  trQ = 0;  gh = 0;  Num_1 = 0;  Num_2 = 0;  Den_1 = 0;
for j = 1:size(b, 2)
    s_oj  = S_o.'  * b(:, j);
    sp_oj = Sp_o.' * b(:, j);
    Den_o = Den_o + 2 * sp_oj.' * W_o_inv * sp_oj - 4 * sp_oj.' * D_o * s_oj;
    g = 2 * S_o * (D_o * s_oj - W_o_inv * sp_oj);         % Num, linear:    g^T n_j
    h = -4 * S_o * D_o * sp_oj;                           % Den, linear:    h^T n_j
    Q = S_o * D_o * S_o.';                                % Num, quadratic: n_j^T Q n_j
    if ADD
        g = g - 2 * Sp_o * W_o_inv * s_oj;
        h = h + 4 * Sp_o * W_o_inv * sp_oj - 4 * Sp_o * D_o * s_oj;
        Q = Q - (Sp_o * W_o_inv * S_o.' + S_o * W_o_inv * Sp_o.');
    end
    if V2
        spp_oj = Spp_o.' * b(:, j);
        Den_o = Den_o + 2 * spp_oj.' * W_o_inv * s_oj + s_oj.' * Z_o * s_oj;
        h = h + 2 * Spp_o * W_o_inv * s_oj + 2 * S_o * W_o_inv * spp_oj + 2 * S_o * Z_o * s_oj;
    end
    trQ = trQ + trace(Q);  gh = gh + g.' * h;
    if ~isempty(n)
        Num_1 = Num_1 + g.' * n(:, j);  Num_2 = Num_2 + n(:, j).' * Q * n(:, j);  Den_1 = Den_1 + h.' * n(:, j);
    end
end
E_u = trQ / Den_o - gh / Den_o^2;
dl = [];
if ~isempty(n), dl = (Num_1 + Num_2) / Den_o - Num_1 * Den_1 / Den_o^2; end
end

function R = mc_chunk(rs, NK, SIG2, P, Pc, B, l0, SF, EPS, S, S1, S2, BIAS, REFIT, Pm, l0m)
% NK realizations of the serial loop body from stream rs; columns: dl_lin, dl_res, dl_doc, dl_doc_add,
% dl_doc_v2, dl_doc_b, dl_doc_b_add, dl_doc_b_v2, dl_fit, (spare).
R = nan(NK, 10);
for m = 1:NK
    nij = sqrt(SIG2) * randn(rs, size(B));                % n_ij
    Bn  = B + nij;
    [~, Jt1, Jt2] = jcost(P, Pc, Bn, l0, SF);
    R(m,1) = -Jt1 / Jt2;
    R(m,2) = EPS.' * nij(:);
    R(m,3) = doc_dl(S, S1, B, nij);
    R(m,4) = doc_dl(S, S1, B, nij, 'add');
    R(m,5) = doc_dl(S, S1, B, nij, 'v2', S2);
    [~, R(m,6)] = doc_bias(S, S1, B, nij);
    [~, R(m,7)] = doc_bias(S, S1, B, nij, 'add');
    [~, R(m,8)] = doc_bias(S, S1, B, nij, 'v2', S2);
    if REFIT
        if BIAS
            lf = fitl(P, Pc, Bn, l0, SF);
        else
            [~, lf] = fitting(Pm, Bn, Pc, l0m, false);  lf = lf*1e6;
        end
        R(m,9) = lf - l0;
    end
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
