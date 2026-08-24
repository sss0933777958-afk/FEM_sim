%% svd_sigma_hung_ws.m -- 純後處理：讀 hung R300/R700 current 校正結果 → 空間 SVD → 存 σ1/σ2/σ3 .mat
% =========================================================================
%  分工（見 calibration-shared-structure 規則）：
%    校正 = main.m + function/（產 calib_current_<V>_R150_eighteen.mat 到 data/）。
%    本腳本（utils/<model>/）= 純後處理：**讀** calib .mat → 逐節點解 T=S(p)·Ĥ_I 的 SVD →
%       **存** σ1/σ2/σ3 到 data/。**不重跑校正**。
%    畫圖 = plot/hung_hexapole/current/plot_sigma_hist_hung_ws.m（讀本腳本產的 .mat）。
%
%  每節點 T = S(p)·Ĥ_I 為 3×6 → svd 得 σ1≥σ2≥σ3（單位 mT/A）。
%    S(p) = (p/ℓ̂ − Pc) ./ |·|³（無因次 Coulomb kernel，Pc = biased 電荷格 = Pc_base+E(e)）。
%    Ĥ_I  = ĝ_I·K̄_I [mT/A]。取樣 = R≤150µm 球內均勻 3D 網格（非 FEM 節點，避免密度偏差）。
%  輸出：data/hung_hexapole/.mat/svd_sigma_hung_ws_R150.mat（自描述，供 plot 讀）。
% =========================================================================

clear; clc;

R_EVAL   = 150e-6;                 % SVD 取樣球半徑 [m]（= 校正擬合區，兩變體可比）
VARIANTS = {'R300','R700'};
MODEL    = 'hung_hexapole';

%% ---- paths（相對自身定位）----
here  = fileparts(mfilename('fullpath'));            % utils/hung_hexapole
CAL   = fileparts(fileparts(here));                  % Calibration_using_FEM_modeling
addpath(fullfile(CAL,'function'));                   % model_config 等
addpath(fullfile(CAL,'common_path'));
matdir = fullfile(CAL,'data',MODEL,'.mat');

cfg = model_config(MODEL, '');                       % 取 Pc_base（幾何常數）
Pc_base = cfg.Pc_base;

%% ---- 逐變體：讀 calib → ball_grid → 逐節點 svd → 收 σ1/σ2/σ3 ----
sig = struct('variants',{VARIANTS}, 'R_eval_um',round(R_EVAL*1e6), 'unit','mT/A');
P = ball_grid(R_EVAL, 32);                           % 兩變體同取樣點雲（WP 原點）
Np = size(P,1);
for c = 1:numel(VARIANTS)
    V = VARIANTS{c};
    src = fullfile(matdir, sprintf('calib_current_%s_R%03d_eighteen.mat', V, round(R_EVAL*1e6)));
    assert(exist(src,'file')==2, '缺校正 .mat：%s（先跑 main.m，MODEL=hung、VARIANT=%s）', src, V);
    rec = load(src, 'l_hat','e','KI_bar','gI_hat','C_mean','kappa_mean');

    Pc    = make_Pc_local(rec.e, Pc_base);           % 3×6 biased 電荷格
    Hhat  = rec.gI_hat * rec.KI_bar;                 % 6×6 Ĥ_I [mT/A]
    ell_m = rec.l_hat;

    s1 = zeros(Np,1);  s2 = zeros(Np,1);  s3 = zeros(Np,1);
    for i = 1:Np
        Dk = P(i,:).'/ell_m - Pc;                    % 3×6
        sv = svd((Dk ./ (vecnorm(Dk).^3)) * Hhat);   % σ1≥σ2≥σ3
        s1(i) = sv(1);  s2(i) = sv(2);  s3(i) = sv(3);
    end

    % sanity：∏σ ≈ rec.C_mean、σ3/σ1 ≈ rec.kappa_mean（同區、取樣法略異，應同數量級）
    Cchk = mean(s1.*s2.*s3);  Kchk = mean(s3./s1);
    fprintf('[%s] N=%d | mean σ=[%.3f %.3f %.3f] mT/A | ∏σ=%.4g (calib %.4g) | σ3/σ1=%.4f (calib %.4f)\n', ...
            V, Np, mean(s1), mean(s2), mean(s3), Cchk, rec.C_mean, Kchk, rec.kappa_mean);

    sig.(V) = struct('sigma1',s1, 'sigma2',s2, 'sigma3',s3, 'npts',Np, ...
                     'variant',V, 'ell_um',ell_m*1e6, 'gI_hat',rec.gI_hat, ...
                     'C_mean_svd',Cchk, 'kappa_svd',Kchk, 'src',src);
end

%% ---- 存 σ .mat（自描述，供 plot 讀）----
if ~exist(matdir,'dir'), mkdir(matdir); end
outmat = fullfile(matdir, sprintf('svd_sigma_hung_ws_R%d.mat', round(R_EVAL*1e6)));
save(outmat, '-struct', 'sig');
fprintf('saved %s\n', outmat);

%% ======================= local functions ==================================
function Pc = make_Pc_local(e17, Pc_base)
% 電荷格 Pc = Pc_base + E(e)（含 e6z 約束；忠實複製 function/fitting.m 的 make_Pc）
    E = zeros(3,6);
    E(:,1)=e17(1:3); E(:,2)=e17(4:6); E(:,3)=e17(7:9); E(:,4)=e17(10:12); E(:,5)=e17(13:15);
    E(1,6)=e17(16);  E(2,6)=e17(17);  E(3,6)=e17(1)-e17(4)+e17(8)-e17(11)+e17(15);
    Pc = Pc_base + E;
end

function P = ball_grid(R, nr)
% 球內均勻 3D 網格（WP 原點）；點數只由 nr 決定、與 R 無關（避免密度偏差）
    h = R/nr;  v = -R:h:R;
    [X,Y,Z] = ndgrid(v,v,v);
    in = (X.^2 + Y.^2 + Z.^2) <= R^2;
    P = [X(in), Y(in), Z(in)];
end
