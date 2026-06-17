%% main.m -- Hall-sensor-based hexapole model: per-pole d, 18-param bias 版 (no_fix_dir)
%  跟 Hall_sensor_base_fix_dir 一樣求 sensor per-pole 常數 d,但 page-1 校正改用
%  18 參數 bias 模型(no_fix_dir,actuator 框、電荷可離軸 Pc_18)而非 fix-ℓ。
%
%  Pipeline:
%    page-1  載入既有 calib_bias.mat(no_fix_dir 校正結果:ell_hat, Pc_18, R, F)
%            → 重載 6-coil FEM 場 → R_select 近場球 → 旋進 actuator 框 → build_A → M, c
%    page-2  build_sensor_geometry → extract_Vmat(all-source) → solve_d(含增益 g_H)
%            → sensor_residual_bias(actuator 框殘差)
%    存     calib_sensor_d_no_fix_dir.mat(不蓋 fix_dir 版的 calib_sensor_d.mat)
%    LaTeX  d_v2 / d_final（write_d_tex）、KH_v2 / KH_final（compute_KH+write_KH_tex）→ results/
%
%  Model    : b_ij = g_H · S_i V_j d ;  g_H = 1/(4πℓ̂²);  charges at actuator-frame Pc_18.
%  Current  : I = 1 A = FEM excitation (per fit-current-matches-sim rule).
%  Sign     : 物理 signed B·n+,all-source(翻下極激發 P1/P3/P6).
%  All model math lives in code\function\ ; this file is just the driver.

clear; clc;

%% ---- config ----------------------------------------------------------------
R_select = 150e-6;             % 近場取點半徑 [m](= calib_bias.mat 的 R_select)
I_actual = 1;                  % drive current [A] = FEM excitation (1 A)
S_hall   = 130;               % Hall 靈敏度 [V/T](EQ-730L)
N_I      = 6;                  % FEM 模擬次數 = 6 個單線圈解
dataset  = 'wp';              % page-1 建 M,c 用近場 'wp'(=R_select 球內,與 calib_bias 一致)

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_no_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));                                                      % 模型輔助函式
results_root  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
charge_dir    = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit';
calib_bias_in = fullfile(charge_dir, 'calib_bias.mat');
mat_out       = fullfile(charge_dir, 'calib_sensor_d_no_fix_dir.mat');
tex_dir       = fullfile(TREE,'results');
if ~exist(tex_dir,'dir'); mkdir(tex_dir); end

%% ---- 載入 18-param 校正結果(no_fix_dir) ------------------------------------
assert(exist(calib_bias_in,'file')==2, 'calib_bias.mat 不存在,請先跑 no_fix_dir 校正');
CB = load(calib_bias_in);                              % R, Pc_18, ell_hat, F, apdl_to_paper_idx, R_select, ...
ell_hat = CB.ell_hat;                                  % 18-param ℓ̂ [m]
Pc      = CB.Pc_18;                                    % 3x6 bias 電荷格點(actuator 框)
Rrot    = CB.R;                                        % 3x3 measure→actuator 旋轉
F       = CB.F;                                        % 6x6 電流矩陣
apdl_to_paper_idx = CB.apdl_to_paper_idx;             % coil→paper pole
if isfield(CB,'R_select'), R_select = CB.R_select; end % 與校正一致
fprintf('載入 calib_bias.mat:ell_hat=%.4f mm,Pc_18(actuator),R(det=%.4f)\n', ell_hat*1e3, det(Rrot));

cnst = mt_constants();

%% ===========================================================================
%  PAGE 1 - 重載 6-coil FEM、選近場球、旋進 actuator 框、建 M 與 c
%  ===========================================================================
d1   = import_ansys_data(fullfile(results_root, 'coil1', 'standard'), dataset, 'coil1');
air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
zwp1 = d1.z - cnst.SPH_OFST;
P_all = [d1.x(air1), d1.y(air1), zwp1(air1)];          % measure 框
insel = sum(P_all.^2,2) <= R_select^2;                % 近場球(measure 框)
P     = (Rrot * P_all(insel,:).').';                  % 旋進 actuator 框 → Np x 3
Np    = size(P,1);
fprintf('PAGE1: R_select=%g um 內選出 N_p=%d 點(actuator 框)\n', R_select*1e6, Np);

B = zeros(Np,3,N_I);                                   % b_ij(actuator 框,= -FEM all-source)
for j = 1:N_I
    if j == 1, dj = d1; airj = air1;
    else
        cn = sprintf('coil%d', j);
        dj = import_ansys_data(fullfile(results_root, cn, 'standard'), dataset, cn);
        airj = filter_iron_nodes(dj.x,dj.y,dj.z,cnst,struct('visualize',false));
    end
    Bj_all = -[dj.bx(airj), dj.by(airj), dj.bz(airj)]; % all-source 整體變號
    B(:,:,j) = (Rrot * Bj_all(insel,:).').';          % 旋進 actuator 框
end
Bstack = zeros(3*Np, N_I);                            % 堆疊(點-major、xyz 交錯)
for j = 1:N_I, Bstack(:,j) = reshape(B(:,:,j).', 3*Np, 1); end

A = build_A(ell_hat, Pc, P);                          % 3Np x 6 空間函數矩陣(actuator)
M = A.' * A;                                          % 6x6
c = A.' * Bstack;                                     % 6xN_I(各欄 c_j = A^T b_j)
fprintf('PAGE1: 已建 A(%dx6)、M(6x6)、c(6x%d)\n', 3*Np, N_I);

%% ===========================================================================
%  PAGE 2 - Hall-sensor 模型 + 求 d
%  ===========================================================================
[sensor_pos, sensor_n, disc_u, disc_v, disc_local, Ndisc] = build_sensor_geometry(cnst);
[Vmat, exc_sign] = extract_Vmat(results_root, cnst, apdl_to_paper_idx, ...
                                sensor_pos, sensor_n, disc_u, disc_v, disc_local, Ndisc, S_hall);
[d, gH]      = solve_d(Vmat, exc_sign, M, c, ell_hat, cnst, N_I);
nrmse_sensor = sensor_residual_bias(A, Bstack, Vmat, exc_sign, d, gH, N_I);

fprintf('\n=============== PAGE 2: Hall-sensor 模型 (18-param bias, ell_hat=%.3f mm) ===============\n', ell_hat*1e3);
fprintf('  sensor 電壓 Vmat [V](列=sensor 極 P1..P6,欄=激發 coil1..6):\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat(i,:)); end
fprintf('  d (6x1, P1..P6, 含增益 d_final):\n'); fprintf('   % .4e\n', d);
fprintf('  sensor 模型在 R<=%dum 的相對 RMSE = %.2f%%\n', round(R_select*1e6), nrmse_sensor);
fprintf('========================================================================\n');

%% ---- 存 calib_sensor_d_no_fix_dir.mat(新檔名,不蓋 fix_dir 版)----
save(mat_out, 'd','gH','Vmat','ell_hat','Pc','Rrot','nrmse_sensor','S_hall', ...
              'sensor_pos','sensor_n','apdl_to_paper_idx','R_select');
fprintf('已存 %s\n', mat_out);

%% ===========================================================================
%  後處理 + LaTeX:最終解 d_final + KH_final
%  ===========================================================================
SS = struct('ell_hat',ell_hat,'R_select',R_select,'S_hall',S_hall, ...
            'nrmse_sensor',nrmse_sensor,'gH',gH,'Vmat',Vmat, ...
            'apdl_to_paper_idx',apdl_to_paper_idx);

d_final = d;                                                 % 最終解(含增益 g_H)
Vp = zeros(6,6);                                              % 欄重排:激發 coil → 激發 paper pole
for kc = 1:6, Vp(:, apdl_to_paper_idx(kc)) = Vmat(:, kc); end

write_d_tex(fullfile(tex_dir,'d_final.tex'), 'final', d_final, Vp, gH, SS);

mu0 = cnst.mu_0;  Nc = cnst.N_c;
[gF, KH, Ra] = compute_KH(d_final, Vmat, F, mu0, Nc);
fprintf('[final] g_F = %.4e,  R_a = %.4e,  K_H(1,1) = %.4f\n', gF, Ra, KH(1,1));
write_KH_tex(fullfile(tex_dir,'KH_final.tex'), 'final', 'hw.pdf (gain $d$, $b=g_H S V d$, 18-param bias)', d_final, gF, KH, Ra, SS);

fprintf('done: 2 .tex (d_final, KH_final) in %s\n', tex_dir);
