%% main.m -- Hall-sensor-based hexapole model: per-pole d, 18-param bias 版 (no_fix_dir)
%  跟 Hall_sensor_base_fix_dir 一樣求 sensor per-pole 常數 d,但 page-1 校正改用
%  18 參數 bias 模型(no_fix_dir,actuator 框、電荷可離軸 Pc_18)而非 fix-ℓ。
%
%  Pipeline:
%    page-1  載入既有 calib_bias.mat(no_fix_dir 校正結果:ell_hat, Pc_18, R, F)
%            → 重載 6-coil FEM 場 → R_select 近場球 → 旋進 actuator 框 → build_S 逐點 → M, c
%    page-2  build_sensor_geometry → extract_Vmat(all-source,真實節點) → solve_d(no-gain)
%            → sensor_residual_bias(actuator 框 cost J)
%    存     fitting_d/calib_sensor_d_no_fix_dir.mat(不蓋 fix_dir 版的 calib_sensor_d.mat)
%
%  Model    : b_ij = S_i V_j d  (no-gain,無 g_H / K_H);  charges at actuator-frame Pc_18.
%  Output   : 與 Hall_sensor_base_fix_dir 對齊 —— 只回 d、Vmat、cost J = Σ‖ε‖²(無 RMSE)。
%  Current  : I = 1 A = FEM excitation (per fit-current-matches-sim rule).
%  Sign     : 物理 signed B·n+,all-source(翻下極激發 P1/P3/P6).
%  All model math lives in code\function\ ; this file is just the driver.

clear; clc;

%% ---- config ----------------------------------------------------------------
R_select = 150e-6;             % 近場取點半徑 [m](= calib_bias.mat 的 R_select)
I_actual = 1;                  % drive current [A] = FEM excitation (1 A)
S_hall   = 130;               % Hall 靈敏度 [V/T](EQ-730L)
N_I      = 6;                  % FEM 模擬次數 = 6 個單線圈解
VARIANT  = 'standard';        % 讀哪個 FEM 變體子夾：'standard'(baseline) | 'sensref'(sensor 加密)；與 fix_dir 對齊
dataset  = 'wp';              % page-1 建 M,c 用近場 'wp'(=R_select 球內,與 calib_bias 一致)

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_no_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));                                                      % 模型輔助函式
results_root  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
charge_dir    = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit';
calib_bias_in = fullfile(charge_dir, 'calibration', 'calib_bias.mat');
mat_out       = fullfile(charge_dir, 'fitting_d', 'calib_sensor_d_no_fix_dir.mat');  % 歸到 charge_fit/fitting_d/
if ~exist(fileparts(mat_out),'dir'); mkdir(fileparts(mat_out)); end

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
% 建 M、c（用 build_S 逐點，取代 build_A 堆疊版；M=Σ S_iᵀS_i、c_j=Σ S_iᵀ·頁1 負號版場）
M = zeros(6,6); c = zeros(6, N_I);
for i = 1:Np
    Si = build_S(P(i,:), ell_hat, Pc);                % 3×6 空間核（actuator 框、Pc_18）
    M  = M + Si.' * Si;                               % Σ_i S_iᵀ S_i
    for j = 1:N_I
        c(:,j) = c(:,j) + Si.' * squeeze(B(i,:,j)).'; % c_j += S_iᵀ·B(i,:,j)（B = −FEM）
    end
end
fprintf('PAGE1: 已建 M(6x6)、c(6x%d)（build_S 逐點）\n', N_I);

%% ===========================================================================
%  PAGE 2 - Hall-sensor 模型 + 求 d
%  ===========================================================================
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);
[Vmat, exc_sign] = extract_Vmat(results_root, cnst, apdl_to_paper_idx, ...
                                sensor_pos, sensor_n, S_hall, VARIANT);
d = solve_d(Vmat, exc_sign, M, c, N_I);
J = sensor_residual_bias(P, B, Vmat, exc_sign, ell_hat, Pc, d, N_I);

fprintf('\n=============== PAGE 2: Hall-sensor 模型 (18-param bias, ell_hat=%.3f mm) ===============\n', ell_hat*1e3);
fprintf('  sensor 電壓 Vmat [V](列=sensor 極 P1..P6,欄=激發 coil1..6):\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat(i,:)); end
fprintf('  d (6x1, P1..P6, no-gain):\n'); fprintf('   % .4e\n', d);
fprintf('  cost J = %.6e  [T^2]\n', J);
fprintf('========================================================================\n');

%% ---- 存 fitting_d/calib_sensor_d_no_fix_dir.mat(新檔名,不蓋 fix_dir 版)----
%  輸出對齊 Hall_sensor_base_fix_dir：只存 d / Vmat / cost J（+ bias 專屬 Pc/Rrot）；無 g_H / K_H / RMSE。
save(mat_out, 'd','Vmat','exc_sign','ell_hat','Pc','Rrot','J','S_hall', ...
              'sensor_pos','sensor_n','apdl_to_paper_idx','R_select');
fprintf('已存 %s\n', mat_out);
