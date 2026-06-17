%% verify_d_vs_kmQS.m
%  驗證:用 d 解出的 sensor 場  b = g_H·S·diag(d)·V_j
%        是否等於 hw.pdf eq(1) 的  b = (k_m/ℓ̂²)·S·Q̂_j ,
%        其中前面已求的磁荷  Q̂_j = (1/μ0)·diag(d)·V_j 。
%  關係:g_H = k_m/(ℓ̂²·μ0)  ⇒  兩式應逐元素相同(機器精度)。
%  全部用 R=150um 下的 fitting 值(calib_sensor_d.mat)。
clear; clc;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
matf = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\calib_sensor_d.mat';

cnst = mt_constants();
S = load(matf);                       % d(=d_final), gH, Vmat(all-source), ell_hat, R_select
d = S.d; Vmat = S.Vmat; gH = S.gH; ell = S.ell_hat;

%% (1) 橋接常數:g_H  vs  k_m/(ℓ̂²·μ0)
gH_formula = cnst.k_m / (ell^2 * cnst.mu_0);
fprintf('=== (1) gain 常數 ===\n');
fprintf('  g_H (stored)        = %.6e  1/m^2\n', gH);
fprintf('  k_m/(l^2 mu0)       = %.6e  1/m^2\n', gH_formula);
fprintf('  相對差              = %.3e\n\n', abs(gH-gH_formula)/gH);

%% (2) 前面已求的磁荷  Q̂ = (1/μ0)·diag(d)·V  (6x6,欄 j = 第 j 次模擬)
Q = (1/cnst.mu_0) * diag(d) * Vmat;
fprintf('=== (2) 磁荷 Q̂ = (1/mu0) diag(d) Vmat  [Wb·m  等效] ===\n');
for i=1:6, fprintf('   % .4e % .4e % .4e % .4e % .4e % .4e\n', Q(i,:)); end
fprintf('\n');

%% (3) 在 R<=150um 真實節點上比兩種寫法
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);                          % 3x6 極尖單位方向(= build_S 的 Pc)
d1   = import_ansys_data(fullfile(results_root, 'coil1', 'standard'),'wp','coil1');
air  = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
zwp  = d1.z - cnst.SPH_OFST;
Pall = [d1.x(air), d1.y(air), zwp(air)];
insel = sum(Pall.^2,2) <= S.R_select^2;
P = Pall(insel,:); Np = size(P,1);

km_l2 = cnst.k_m / ell^2;                             % k_m/ℓ̂²
maxabs = 0; sse = 0; ssb = 0;
for i = 1:Np
    Di = P(i,:).'/ell - dhat;                        % 3x6  (pbar - pbar_ck)
    Si = Di ./ (vecnorm(Di).^3);                     % 3x6  空間函數 S_i
    for j = 1:6
        b_d = gH    * Si * (Vmat(:,j).*d);           % d 形式: g_H S diag(d) V_j
        b_Q = km_l2 * Si * Q(:,j);                   % eq(1):  (k_m/ℓ̂²) S Q̂_j
        df  = b_d - b_Q;
        maxabs = max(maxabs, max(abs(df)));
        sse = sse + sum(df.^2); ssb = ssb + sum(b_Q.^2);
    end
end
fprintf('=== (3) 在 R<=150um 全 %d 節點 × 6 sim 比較 ===\n', Np);
fprintf('  max|b_d - b_Q|              = %.3e T\n', maxabs);
fprintf('  相對 ||b_d-b_Q||/||b_Q||    = %.3e\n', sqrt(sse/ssb));
fprintf('\n結論: 兩式%s相同(d 解出的場 = k_m·Q̂·S/ℓ̂²)。\n', ...
        ternary(maxabs < 1e-12, '在機器精度內', '不'));

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
