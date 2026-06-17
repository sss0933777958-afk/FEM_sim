%% main.m -- Hall-sensor-based hexapole model: per-pole constant d (求 d 一條龍)
%  Documents: "Calibration using FEM modeling_V2.pdf" (page-1 步驟1-9) + "hw.pdf"/"derivation.pdf"
%  Pipeline (求 d → 後處理 → 驗證 → LaTeX, all in this driver):
%    page-1  載 6-coil FEM 場 → 選 R 球內節點 → 擬合特徵長度 ℓ̂(fminbnd) → 在 ℓ̂ 建 M, c
%    page-2  建 sensor 幾何 → 抽 Vmat(all-source) → 解 d(含增益 g_H) → 殘差 nrmse_sensor
%    存     calib_sensor_d.mat(同原路徑,供舊版下游/驗證沿用)
%    LaTeX  d_v2 / d_final（write_d_tex）、KH_v2 / KH_final（compute_KH + write_KH_tex）→ results/sensor_d/
%  Model    : b_ij = g_H · S_i V_j d ;  g_H = k_m/(ℓ̂²μ0) = 1/(4πℓ̂²)
%  Current  : I = 1 A = FEM excitation (per fit-current-matches-sim rule).
%  Sign     : 物理 signed B·n+,all-source(翻下極激發 P1/P3/P6)(per charge-model-source-convention).
%  All model math lives in code\function\ ; this file is just the driver.
%  改寫自 fixl_fit/code/scripts/calib_fem.m(page-1+page-2)+ sensor_d gen_*_latex(後處理/LaTeX)。

clear; clc;

%% ---- config ----------------------------------------------------------------
R_select   = 150e-6;           % 取點半徑 [m]:只保留 |p|<=R_select 的 air 節點(fit ℓ̂ 用)
I_actual   = 1;                % drive current [A] = FEM excitation (1 A)
S_hall     = 130;              % Hall 靈敏度 [V/T](EQ-730L:130 V/T,非舊版 130e-3)
N_I        = 6;                % FEM 模擬次數 = 6 個單線圈解

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir'];
% mt_constants / import_ansys_data / filter_iron_nodes 在 hexapole-long2016 的 analysis 夾
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
addpath(fullfile(TREE,'code','function'));                                       % 模型輔助函式
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
mat_out      = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\' ...
                'long2016_hexapole_halfcut\charge_fit\calib_sensor_d.mat'];
tex_dir      = fullfile(TREE,'results');
if ~exist(tex_dir,'dir'); mkdir(tex_dir); end
if ~exist(fileparts(mat_out),'dir'); mkdir(fileparts(mat_out)); end

%% ---- constants + pole-tip directions ---------------------------------------
cnst = mt_constants();                                 % 幾何常數(R_norm、極尖、SPH_OFST、k_m、mu_0、N_c ...)
apdl_to_paper_idx = [1,3,6,5,2,4];                     % APDL coil j → paper pole 索引
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];  % 3x6 極尖位置(WP 座標)
dhat = tip ./ vecnorm(tip);                            % 3x6 極尖單位方向
Pc   = dhat;                                           % 電荷(極)位置,放單位距離(build_S 用)

F = zeros(6, N_I);                                     % 6x6 電流矩陣(pole 空間,置換)
for j = 1:N_I, F(apdl_to_paper_idx(j), j) = 1; end     % 第 j 次:paper pole 通 1A
assert(rank(F) == 6, '電流矩陣 F 必須 rank = 6');

%% ===========================================================================
%  PAGE 1 - 載 6-coil FEM 場、選球內節點、擬合 ℓ̂、在 ℓ̂ 建 M 與 c
%  ===========================================================================
% 步驟1: coil1 air 節點當工作點 p_i,取 R_select 球內(6 顆 coil 共用同一 mesh/順序)
d1   = import_ansys_data(fullfile(results_root, 'coil1', 'standard'),'wp','coil1');
air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
zwp1 = d1.z - cnst.SPH_OFST;
P_all = [d1.x(air1), d1.y(air1), zwp1(air1)];
insel = sum(P_all.^2,2) <= R_select^2;
P     = P_all(insel,:);  Np = size(P,1);
fprintf('PAGE1 步驟1: R_select = %g um 內選出 N_p = %d 個點\n', R_select*1e6, Np);

% 步驟3: 記錄 b_ij = -B^FEM(all-source 慣例,整體變號);存 B(Np,3,N_I)
B = zeros(Np,3,N_I);
for j = 1:N_I
    if j == 1, dj = d1; airj = air1;
    else
        cn = sprintf('coil%d', j);
        dj = import_ansys_data(fullfile(results_root, cn, 'standard'),'wp',cn);
        airj = filter_iron_nodes(dj.x,dj.y,dj.z,cnst,struct('visualize',false));
    end
    Bj_all = -[dj.bx(airj), dj.by(airj), dj.bz(airj)];
    B(:,:,j) = Bj_all(insel,:);
end
fprintf('PAGE1 步驟3: 已記錄 b_ij(%d 模擬 × %d 點)\n', N_I, Np);

% 步驟5-8: 擬合特徵長度 ell_hat(g_j 已剖面化解掉 → 一維 fminbnd)
obj    = @(x) cost_J(x, P, B, Pc);
ell_lo = 0.2e-3;  ell_hi = 2.0e-3;
fopt   = optimset('TolX',1e-9,'Display','iter');
[ell_hat, Jmin] = fminbnd(obj, ell_lo, ell_hi, fopt);
fprintf('PAGE1 步驟8: ell_hat = %.4f mm  (Jmin = %.6e)\n', ell_hat*1e3, Jmin);

% 步驟9: 在 ell_hat 重建法矩陣 M = Σ S^T S 與右端 c_j = Σ S^T b_ij(page-2 解 d 重用)
M = zeros(6,6);  c = zeros(6,N_I);
for i = 1:Np
    Si = build_S(P(i,:), ell_hat, Pc);
    M  = M + Si.' * Si;
    for j = 1:N_I
        c(:,j) = c(:,j) + Si.' * squeeze(B(i,:,j)).';
    end
end

%% ===========================================================================
%  PAGE 2 - Hall-sensor 模型 + 求 d
%  ===========================================================================
[sensor_pos, sensor_n, disc_u, disc_v, disc_local, Ndisc] = build_sensor_geometry(cnst);
[Vmat, exc_sign] = extract_Vmat(results_root, cnst, apdl_to_paper_idx, ...
                                sensor_pos, sensor_n, disc_u, disc_v, disc_local, Ndisc, S_hall);
[d, gH]      = solve_d(Vmat, exc_sign, M, c, ell_hat, cnst, N_I);
nrmse_sensor = sensor_residual(P, B, Vmat, exc_sign, ell_hat, Pc, d, gH, N_I);

fprintf('\n=============== PAGE 2: Hall-sensor 模型 (ell_hat=%.3f mm) ===============\n', ell_hat*1e3);
fprintf('  sensor 電壓 Vmat [V](列=sensor 極 P1..P6,欄=激發 coil1..6):\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat(i,:)); end
fprintf('  d (6x1, P1..P6, 含增益 d_final):\n'); fprintf('   % .4e\n', d);
fprintf('  sensor 模型在 R<=%dum 的相對 RMSE = %.2f%%\n', round(R_select*1e6), nrmse_sensor);
fprintf('========================================================================\n');

%% ---- 存 calib_sensor_d.mat(同原路徑,供 verify_* / 舊下游沿用)----
save(mat_out, 'd','gH','Vmat','ell_hat','nrmse_sensor','S_hall', ...
              'sensor_pos','sensor_n','apdl_to_paper_idx','R_select');
fprintf('已存 %s\n', mat_out);

%% ===========================================================================
%  後處理 + LaTeX:最終解 d_final（write_d_tex）+ KH_final（compute_KH）
%  ===========================================================================
SS = struct('ell_hat',ell_hat,'R_select',R_select,'S_hall',S_hall, ...
            'nrmse_sensor',nrmse_sensor,'gH',gH,'Vmat',Vmat, ...
            'apdl_to_paper_idx',apdl_to_paper_idx);            % 給 writers 的設定 struct

d_final = d;                                                  % 最終解(含增益 g_H)
Vp = zeros(6,6);                                               % 欄重排:激發 coil → 激發 paper pole
for kc = 1:6, Vp(:, apdl_to_paper_idx(kc)) = Vmat(:, kc); end

write_d_tex(fullfile(tex_dir,'d_final.tex'), 'final', d_final, Vp, gH, SS);

mu0 = cnst.mu_0;  Nc = cnst.N_c;                               % derivation.pdf 的 K_H 流程
[gF, KH, Ra] = compute_KH(d_final, Vmat, F, mu0, Nc);
fprintf('[final] g_F = %.4e,  R_a = %.4e,  K_H(1,1) = %.4f\n', gF, Ra, KH(1,1));
write_KH_tex(fullfile(tex_dir,'KH_final.tex'), 'final', 'hw.pdf (gain $d$, $b=g_H S V d$)', d_final, gF, KH, Ra, SS);

fprintf('done: 2 .tex (d_final, KH_final) in %s\n', tex_dir);
