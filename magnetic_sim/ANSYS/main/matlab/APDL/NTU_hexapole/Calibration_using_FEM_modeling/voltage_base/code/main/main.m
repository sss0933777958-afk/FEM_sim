%% main.m -- NTU Hall-sensor 校正 driver（voltage_base）：輸出 V̄, D̄, ^Bĝ_V
%  NTU full_assembly 6-coil：電荷模型 bias 擬合 WP 場 → 6 S2 sensor 1000 點內插 V̄ → D̄/^Bĝ_V。
%  ⚠ NTU 幾何（扁平板、WP 在 (0,0,9.30mm)、無 magic-angle）：mt_constants/load_coils_actuator/
%    build_sensor_geometry/build_V_matrix 已改成 NTU 版；charge fit bias 殘差 ~4%（sanity）。
% =========================================================================
%  解 V̄（sensor 電壓矩陣）、D̄、^Bĝ_V（電壓側）。Ĥ_V 內部算來得 D̄，但不輸出 PDF。
%  USE_BIAS 統一 fix（在軸）/ no_fix（18-param bias）：幾何擬合沿用 current_base 的 fitting
%  （變數投影，逐激發 profile 電荷 g_j），USE_BIAS 只決定 ê 是否在非線性參數集。
%  記號統一論文《Lumped-Parameter…》step 9：b_ij = S(p_i/ℓ̄[,ê])·G_j，G profiled（code Dv=G=D^v）。
%    Recover：Ĥ_V = G·Vᵀ(VVᵀ)⁻¹ → ^Bĝ_V=(6/5)h₁₁、D̄=(5/(6h₁₁))Ĥ_V（solve_D_bar_gain）。
%  Sign  : all-source（flip-sink：只翻下極 sink P1/P3/P6；上極不翻）→ D^v 對角全正。
%  Current: I = 1 A = FEM 激發（per fit-current-matches-sim）。
%
%  管線（driver 只 call 函式）：
%    ① 載參數  ② 沒疊代 bias 就 ê=0
%    ③ build_S_matrix（每點 S_matrix，沿用 current_base）
%    ④ fitting（沿用 current_base，回 ℓ̄/ê/Pc/J）→ profile G=D^v
%    ⑤ build_V_matrix（sensor 電壓 V̄）→ solve_D_bar_gain（D̄、^Bĝ_V、Ĥ_V）
%    ⑥ save .mat + emit_model_results（+18-param 加印 ê）+ gen_B_matrix（另出 B 矩陣 PDF）
% =========================================================================

clear; clc;

%% ---- ① 參數 --------------------------------------------------------------
USE_BIAS   = true;             % NTU：用 18-param bias（single-param 對平板殘差 ~23%，bias ~4%）
VARIANT    = 'full_assembly';  % NTU FEM 變體：data/full_assembly/coilN/
N_I        = 6;                % FEM 模擬次數 = 6 個單線圈解
R_select   = 300e-6;           % 取點半徑 [m]（NTU WP：150µm 只 945 節點、300µm ~7000 且 fit 較穩）
I_actual   = 1;                % 驅動電流 [A] = FEM 激發（1 A）
TURNS      = 50;               % 線圈匝數（記錄用；per-Ampere 模型不引用）
S_hall     = 130;              % Hall 靈敏度 [mV/mT]
ell0       = 0.5e-3;           % ℓ̄ 初值 [m]（fitting 在 SI）
n_uniform  = 1000;             % sensor 圓柱內內插取樣點數（使用者指定 1000）
AXIAL_TOL  = 0.10e-3;          % sensor 圓柱厚度 [m]：canonical 0.10e-3
htag       = '';  if abs(AXIAL_TOL-0.10e-3)>1e-12, htag = sprintf('_h%gum', AXIAL_TOL*1e6); end
dataset    = 'all';

%% ---- paths ---------------------------------------------------------------
CAL  = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'NTU_hexapole\Calibration_using_FEM_modeling'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');     % mt_constants/import_ansys_data/filter_iron_nodes
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');  % ansys_path
addpath(fullfile(CAL,'current_base','code','main_function'));     % load_coils_actuator/select_ball/build_S_matrix/fitting
addpath(fullfile(CAL,'voltage_base','code','main_function'));     % build_sensor_geometry/build_V_matrix/solve_D_bar_gain/emit_model_results/gen_B_matrix
model        = 'NTU_hexapole';
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';
mesh_csv_dir = fullfile(fileparts(results_root),'csv');   % sensor-local tet CSV（basegap 網格）
if ~strcmp(VARIANT,'gap_200um') && ~strcmp(VARIANT,'no_gap'), mesh_csv_dir = fullfile(mesh_csv_dir,VARIANT); end
TREE         = fullfile(CAL,'voltage_base');
mat_dir      = fullfile(TREE,'data');      if ~exist(mat_dir,'dir'); mkdir(mat_dir); end

%% ---- 常數 + 慣例 ---------------------------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1,3,6,5,2,4];

%% ---- ③ 載 6-coil FEM（actuator 框）→ 選 R 球 → all-source（② ê 由 USE_BIAS 決定）
fprintf('載 6 coils（variant ''%s'', dataset ''%s'', USE_BIAS=%d）...\n', VARIANT, dataset, USE_BIAS);
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, dataset, VARIANT);
[P, Bstack, npts] = select_ball(D, R_select);
fprintf('選出 N_p = %d 個 air 節點（R ≤ %g µm 球內）\n', npts, R_select*1e6);
s_sink = ones(1,6);                                          % all-source（只翻下極 sink P1/P3/P6）
for j = 1:6, if ismember(apdl_to_paper_idx(j), [1 3 6]), s_sink(j) = -1; end; end
Bstack = (-Bstack) .* s_sink;                                % 還原 raw 後只翻下極 → all-source（Bstack 已 mT）

%% ---- ④ fitting：擬合 ℓ̄（+ê 若 USE_BIAS）→ profile G=D^v ------------------
[ell_hat, e_hat, Pc, J] = fitting(P, Bstack, D.Pc_base, ell0, USE_BIAS);
fprintf('擬合：ℓ̄ = %.2f µm（USE_BIAS=%d）| J = %.6e | ‖ê‖ = %.4e\n', ell_hat*1e6, USE_BIAS, J, norm(e_hat));
A   = build_S_matrix(ell_hat, Pc, P);   % 3Np×6 stacked S_matrix
Dv  = (A.' * A) \ (A.' * Bstack);        % G=D^v：Dv_j=(ΣSᵀS)⁻¹(ΣSᵀb_ij)（6×N_I）
ell_hat = ell_hat * 1e6;                 % m → µm（此後 print/save/PDF 用 µm）

%% ---- ⑤ sensor 電壓 V̄（build_V_matrix）→ 解 D̄、^Bĝ_V（solve_D_bar_gain）----
[sensor_pos, sensor_n] = build_sensor_geometry(cnst);
[Vmat, exc_sign] = build_V_matrix(results_root, cnst, apdl_to_paper_idx, ...
                       sensor_pos, sensor_n, S_hall, mesh_csv_dir, n_uniform, [], AXIAL_TOL, VARIANT);
% 自激發（對角）為正：翻下極激發欄 P1/P3/P6。G(Dv) 與 V 同翻 → H_V=G·Vᵀ(VVᵀ)⁻¹ 不變（D̄ 全正不動），
%   G 與 V 的自激發對角同時轉正（激發極自己的電壓/電荷為正）。
col_sign = ones(1,6);
for j = 1:6, if ismember(apdl_to_paper_idx(j), [1 3 6]), col_sign(j) = -1; end; end
Dv   = Dv   .* col_sign;
Vmat = Vmat .* col_sign;
[D_bar, ghat_V_B, H_V] = solve_D_bar_gain(Dv, Vmat);   % Ĥ_V 內部算、不印（D̄ 對 col_sign 不變）
[U_hv, S_hv, W_hv] = svd(H_V);
kappa = S_hv(1,1) / S_hv(6,6);

%% ---- 區域場誤差 + 重建檢查 -----------------------------------------------
errpct = 100 * sqrt(J / sum(Bstack(:).^2));
recon  = norm(H_V*Vmat - Dv,'fro') / norm(Dv,'fro');

%% ---- 激發欄重排成 paper P1..P6 + ê(3×6)（存檔/下游用）--------------------
[~, paper_to_apdl] = sort(apdl_to_paper_idx);   % = [1 5 2 6 4 3]
Vmat_p = Vmat(:, paper_to_apdl);
Dv_p   = Dv(:,   paper_to_apdl);
E36    = Pc - D.Pc_base;                         % 3×6 ê 偏移（含 e6z 約束；fix 時全 0）

%% ---- 印結果（主輸出 V̄, D̄, ^Bĝ_V）---------------------------------------
fprintf('\n========= Hall-sensor（USE_BIAS=%d, variant=%s, l_hat=%.2f µm）=========\n', USE_BIAS, VARIANT, ell_hat);
fprintf('  N_p=%d | J=%.4e | region err=%.3f%% | recon ||H_V*V-G||/||G||=%.2e\n', npts, J, errpct, recon);
fprintf('  V [mV]（列=sensor P1..P6，欄=激發 P1..P6）：\n');
for i=1:6, fprintf('   % .3e % .3e % .3e % .3e % .3e % .3e\n', Vmat_p(i,:)); end
fprintf('  D_bar（gauge D_bar(1,1)=5/6）：\n');
for i=1:6, fprintf('   % .4f % .4f % .4f % .4f % .4f % .4f\n', D_bar(i,:)); end
fprintf('  ghat_V_B (^B g_V) = %.4f x10^-3 mT/mV\n', ghat_V_B*1e3);
fprintf('  SVD(H_V) σ = [% .4e % .4e % .4e % .4e % .4e % .4e]，κ=σ1/σ6=%.3g\n', diag(S_hv), kappa);
fprintf('===================================================================\n');

%% ---- 存 .mat 到本包 data/（含 ê/E36 供下游；PDF 前先存）------------------
mat_out = fullfile(mat_dir, sprintf('calib_D_%s%s.mat', VARIANT, htag));
Dmat = H_V;  g_V = ghat_V_B;     % alias：維持 .mat field 名 'Dmat'/'g_V'
save(mat_out, 'Dmat','D_bar','g_V','Dv_p','Vmat_p','exc_sign','ell_hat','e_hat','E36','Pc', ...
              'J','errpct','recon','S_hall','R_select','npts','VARIANT','USE_BIAS','TURNS', ...
              'apdl_to_paper_idx','paper_to_apdl','sensor_pos','sensor_n','n_uniform', ...
              'U_hv','S_hv','W_hv','kappa');
fprintf('已存 %s\n', mat_out);

%% ---- ⑥ 打包 PDF -> results/{single_param|eighteen_param}/ ----------------
emit_model_results(VARIANT, USE_BIAS);    % 主結果 PDF（D̄/G/V/ℓ̄/^Bĝ_V，18-param 加印 ê）
% gen_B_matrix：NTU 無參考 B_ref（V_out/V_in 對照矩陣）→ 暫不出（主結果已在 emit_model_results）。
