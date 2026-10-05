%% eps_shell.m -- shell-averaged calibration residual vs radius (long2016, current, Nr = 200)
%  e_i      = sum_j ||S_i g_j - b_ij||^2          (point i, 6 excitations)       [mT^2]
%  eps(r_k) = (1/6) * sum of e_i over the 6 axis points of shell k            [mT^2]
%  Input : data/long2016_hexapole_halfcut/.mat/calib_current_maxwell_axshN1201_R150_{single,eighteen}.mat
%  Output: utils/data/<MODEL>/eps_shell.mat  (rk, EPS single/eighteen, e_centre)
%  Plot  : plot/long2016_hexapole_halfcut/current/eps_shell.m

MODEL = 'long2016_hexapole_halfcut';   GEOM = 'tip40um';   VARIANT = 'maxwell';
R = 150e-6;   NR = 200;
TAG = {'single', 'eighteen'};

here = fileparts(mfilename('fullpath'));   CAL = fileparts(fileparts(fileparts(here)));
addpath(fullfile(CAL, 'function'), fullfile(CAL, 'common_path'));
cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
[P, Bs] = conv_design_ws(NR, R, struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, 'frame','actuator', 'raw',raw));
r  = vecnorm(P, 2, 2);   sh = round(r / (R/NR));            % shell index, 0 = centre
rk = (1:NR).' * R / NR;                                      % [m]

EPS = nan(NR, 2);   e_centre = nan(1, 2);
for m = 1:2
    c   = load(fullfile(CAL, 'data', MODEL, '.mat', ['calib_current_maxwell_axshN1201_R150_' TAG{m} '.mat']));
    res = kernel(c.l_hat, pc(c.e, ad.Pc_base), P) * c.G - Bs;
    ei  = sum(reshape(sum(res.^2, 2), 3, []), 1).';
    for k = 1:NR
        s = sh == k;   assert(nnz(s) == 6, 'shell %d has %d points', k, nnz(s));
        EPS(k, m) = sum(ei(s)) / 6;
    end
    e_centre(m) = ei(sh == 0);
    assert(abs((e_centre(m) + 6*sum(EPS(:,m))) / c.J - 1) < 1e-10, 'shell sum does not reproduce J');
end

out = fullfile(CAL, 'utils', 'data', MODEL, 'eps_shell.mat');
save(out, 'MODEL', 'R', 'NR', 'rk', 'EPS', 'e_centre', 'TAG');
fprintf('saved %s\n', out);

% ---- same charge grid and kernel as function/fitting.m ----------------------
function Pc = pc(e17, Pc_base)
    E = zeros(3, 6);
    E(:,1) = e17(1:3);   E(:,2) = e17(4:6);   E(:,3) = e17(7:9);   E(:,4) = e17(10:12);
    E(:,5) = e17(13:15); E(1,6) = e17(16);    E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
function S = kernel(l_hat, Pc, P)
    Np = size(P, 1);   pbar = P / l_hat;   S = zeros(3*Np, 6);
    for k = 1:6
        d = pbar - Pc(:,k).';   r3 = sum(d.^2, 2).^1.5;
        S(:,k) = reshape((d ./ r3).', 3*Np, 1);
    end
end
