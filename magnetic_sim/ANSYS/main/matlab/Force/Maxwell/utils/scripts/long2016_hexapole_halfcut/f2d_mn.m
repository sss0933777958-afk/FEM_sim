%% f2d_mn.m -- summary of the 2D Fourier direction-function sweep M, N = 0..5 (fourier2d.pdf).
%  Fits: physics loss, LAM = 0, single model, Phi_hat 6 dim (PH6), 300 blue points x 6 deg (ls_fit.m, S0M 'lab', MAXI 3000).
%  Tests: TEST S25, the 343 red points (nn_grid_sph_r45.mat) x 1742 directions every 6 deg.
%  Reads  utils/data/long2016_hexapole_halfcut/fit_f2d6_single_<tag>.mat  (Z2: pole-condition null space)
%         utils/data/long2016_hexapole_halfcut/test_f2d6_single_<tag>_all.mat  (S25.e N_j x N_p [%], S25.z 6 x N_j x N_p)
%  tag = lam0_M<M>N<N>, except the eight older pairs (names below).
%  NPC = true: the fits WITHOUT pole conditions (ls_fit.m PC = false, tag lam0_M<M>N<N>_npc) for M >= 2; M = 0, 1 have no
%  pole conditions, so their entries are the ordinary fits -> f2d_mn_npc.mat
if ~exist('NPC','var'), NPC = false; end
%  Output utils/data/long2016_hexapole_halfcut/f2d_mn.mat, all 6 x 6 (row M+1, column N+1):
%    L        basis terms before the pole conditions, L = (N+1)(2M+1) - 2 ceil(M/2)
%    L_free   terms of the fit (columns of Z2, after the pole conditions), per element of each D_k
%    e_mean, e_max            force error over all red points and directions [%]
%    Phi_mean, Phi_med, Phi_max, Phi_min   ||Phi_hat|| over all red points and directions
%  GRD = true: instead, the grid relay of 2026-11-19 (M, N = 1..3, LM mu >= 1e-6, LAM 0, no pole conditions): (1,1) from the
%  5-deg labels (fit_..._lam0_M1N1_lmmu6), every other (M,N) from its upstream (M-1,N) / (M,N-1) with the new terms zero,
%  the upstream of smaller final loss -> f2d_grid.mat, all 3 x 3 (row M, column N):
%    T        terms per element of each D_k;  src  initial value ('labels' or '(m,n)');  loss  final training loss [pN^2];
%    it       LM iterations;  e_blue  mean force error, 300 blue points x 1742 directions (6 deg) [%];  e_red  same, 343 red
%    E_blue   cell, force error per sample (1742 directions x 300 points, single) [%]
%    A, B, C  cell, matrices A_h,mn / B_h,mn of fourier2d_new.pdf Sec. 4.3 (4 x 4 x (M+1) x (N+1) x 6, phi = azimuth,
%             theta from +z) and coefficients (10 x T x 6) of each fit;  DEF = 'new'
if ~exist('GRD','var'), GRD = false; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
DAT  = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut');
if GRD
    fn = @(m, n) sprintf('f2d6_single_lam0_M%dN%d_lmg', m, n);  f11 = 'f2d6_single_lam0_M1N1_lmmu6';
    [T, loss, it, e_blue, e_red] = deal(nan(3));  src = cell(3);  [E_blue, A, B, C] = deal(cell(3));
    for M = 1:3
        for N = 1:3
            tag = fn(M, N);  if M == 1 && N == 1, tag = f11; end
            F = load(fullfile(DAT, ['fit_' tag '.mat']), 'M', 'N', 'C', 'A', 'B', 'DEF', 'fv', 'LH', 'e_q', 'phi_j');
            assert(F.M == M && F.N == N, 'f2d_mn:grd', '%s is M = %d, N = %d', tag, F.M, F.N);
            assert(strcmp(F.DEF, 'new'), 'f2d_mn:def', '%s: not fourier2d_new (convert with ls_fit.m CNV)', tag);
            S = load(fullfile(DAT, ['test_' tag '_all.mat']), 'S25');
            Nj = numel(F.phi_j);  Eb = F.e_q(1:Nj, :);                   % grid directions only (the rest are cell centres)
            T(M,N) = size(F.C, 2);  loss(M,N) = F.fv;  it(M,N) = numel(F.LH) - 1;
            e_blue(M,N) = mean(Eb(:));  e_red(M,N) = mean(S.S25.e(:));  E_blue{M,N} = single(Eb);
            A{M,N} = F.A;  B{M,N} = F.B;  C{M,N} = F.C;
            src{M,N} = 'labels';  fb = inf;                              % the upstream of smaller final loss
            for u = {[M-1 N], [M N-1]}
                m_ = u{1}(1);  n_ = u{1}(2);  if m_ < 1 || n_ < 1, continue, end
                t_ = fn(m_, n_);  if m_ == 1 && n_ == 1, t_ = f11; end
                U = load(fullfile(DAT, ['fit_' t_ '.mat']), 'fv');
                if U.fv < fb, fb = U.fv;  src{M,N} = sprintf('(%d,%d)', m_, n_); end
            end
        end
    end
    M_ = 1:3;  N_ = 1:3;  LAM = 0;  CAL = 'single';  DEF = 'new';
    save(fullfile(DAT, 'f2d_grid.mat'), 'M_', 'N_', 'T', 'src', 'loss', 'it', 'e_blue', 'e_red', 'E_blue', 'A', 'B', 'C', 'LAM', 'CAL', 'DEF');
    disp('e_blue [%] (rows M = 1..3, columns N = 1..3):');  disp(e_blue);
    fprintf('-> %s\n', fullfile(DAT, 'f2d_grid.mat'));
    return
end
OLD  = containers.Map({'4 4', '1 4', '2 4', '3 4', '3 3', '3 2', '3 1', '1 1'}, ...
                      {'lam0', 'lam0_M1', 'lam0_M2', 'lam0_M3', 'lam0_M3_N3', 'lam0_M3_N2', 'lam0_M3_N1', 'lam0_M1N1'});
[L, L_free, e_mean, e_max, Phi_mean, Phi_med, Phi_max, Phi_min] = deal(nan(6));
for M = 0:5
    for N = 0:5
        k = sprintf('%d %d', M, N);
        if isKey(OLD, k), tag = OLD(k); else, tag = sprintf('lam0_M%dN%d', M, N); end
        if NPC && M >= 2, tag = sprintf('lam0_M%dN%d_npc', M, N); end
        F = load(fullfile(DAT, ['fit_f2d6_single_' tag '.mat']), 'Z2', 'M', 'N');
        assert(F.M == M && F.N == N, 'f2d_mn:mn', '%s is M = %d, N = %d', tag, F.M, F.N);
        T = load(fullfile(DAT, ['test_f2d6_single_' tag '_all.mat']), 'S25');
        e = T.S25.e;  nPhi = squeeze(vecnorm(T.S25.z, 2, 1));
        L(M+1,N+1) = size(F.Z2, 1);  L_free(M+1,N+1) = size(F.Z2, 2);
        e_mean(M+1,N+1) = mean(e(:));  e_max(M+1,N+1) = max(e(:));
        Phi_mean(M+1,N+1) = mean(nPhi(:));  Phi_med(M+1,N+1) = median(nPhi(:));  Phi_max(M+1,N+1) = max(nPhi(:));  Phi_min(M+1,N+1) = min(nPhi(:));
    end
end
M_ = 0:5;  N_ = 0:5;  LAM = 0;  CAL = 'single';
OUT = 'f2d_mn.mat';  if NPC, OUT = 'f2d_mn_npc.mat'; end
save(fullfile(DAT, OUT), 'M_', 'N_', 'L', 'L_free', 'e_mean', 'e_max', 'Phi_mean', 'Phi_med', 'Phi_max', 'Phi_min', 'LAM', 'CAL', 'NPC');
disp('L_free (rows M = 0..5, columns N = 0..5):');  disp(L_free);
fprintf('-> %s\n', fullfile(DAT, OUT));
