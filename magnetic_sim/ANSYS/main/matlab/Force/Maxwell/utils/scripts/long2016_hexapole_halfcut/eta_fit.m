%% eta_fit.m -- flux-norm ratio of a fitted direction function to the global optimum,
%  eta = || Phi_hat || / || Phi_hat_opt ||, at every position x direction of an ls_fit.m test file.
%  Inputs (set before running):
%    TSF  test file of ls_fit.m (no .mat), S25.z = Phi_hat (6 x N_q x N_p), e.g. 'test_f2d6_single_lam0_M1N1_lmmu6_all_b6'
%         (the (1,1) fit at the 300 blue points x 6 deg)
%    LBF  global-optimum label file (no .mat) on the same positions, IOPT (6 x NPH x NTH x N_p), default
%         'nn_label_sph_r45_mid_flux_single' (300 blue points x 6 deg, 6-dim Phi_hat)
%    ETN  output tag -> utils/data/long2016_hexapole_halfcut/eta_<ETN>.mat
%  Output: eta (N_q x N_p), phi, theta (N_q x 1 [deg]), P (N_p x 3 [um]), nPhi, nOpt (N_q x N_p), TSF, LBF
if ~exist('TSF','var'), TSF = 'test_f2d6_single_lam0_M1N1_lmmu6_all_b6'; end
if ~exist('LBF','var'), LBF = 'nn_label_sph_r45_mid_flux_single'; end
if ~exist('ETN','var'), ETN = 'f2d11_b6'; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
DAT  = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut');
T = load(fullfile(DAT, [TSF '.mat']), 'S25');  S = T.S25;
L = load(fullfile(DAT, [LBF '.mat']), 'IOPT', 'P_act', 'ph', 'th');
assert(max(abs(L.P_act - S.P), [], 'all') < 1e-9, 'eta_fit:pos', 'positions of the test and label files differ');
N_p  = size(S.P, 1);
nPhi = squeeze(vecnorm(S.z, 2, 1));                            % N_q x N_p
nOpt = nan(size(nPhi));
for q = 1:numel(S.phi)                                         % the label at the same (phi, theta)
    a = abs(L.ph - S.phi(q)) < 1e-9;  b = abs(L.th - S.theta(q)) < 1e-9;
    nOpt(q,:) = vecnorm(reshape(L.IOPT(:, a, b, :), 6, N_p), 2, 1);
end
eta = nPhi ./ nOpt;  phi = S.phi;  theta = S.theta;  P = S.P;
save(fullfile(DAT, ['eta_' ETN '.mat']), 'eta', 'phi', 'theta', 'P', 'nPhi', 'nOpt', 'TSF', 'LBF');
fprintf('eta_%s: %d x %d | mean %.4f, std %.4f, min %.4f, max %.4f -> %s\n', ETN, size(eta), mean(eta(:)), std(eta(:)), ...
        min(eta(:)), max(eta(:)), fullfile(DAT, ['eta_' ETN '.mat']));
