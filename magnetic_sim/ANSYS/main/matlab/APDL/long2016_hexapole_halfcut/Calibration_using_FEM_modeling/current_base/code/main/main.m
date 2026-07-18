%% main.m -- hexapole point-charge model calibration driver (current_base)
%  解出 K̄_I、G、^Bĝ_I（電流側）。USE_BIAS 統一 fix（在軸）/ no_fix（18-param bias）：
%    USE_BIAS = false -> fix   : 電荷在極軸（ê=0），只擬合 ℓ̄。
%    USE_BIAS = true  -> no_fix: 18-param bias（ℓ̄ + ê(17)），電荷離軸。
%  Model : B(p) = ^Bg_I * sum_i (sum_j K_bar_ij I_j) * S(p); 電荷位置 ℓ̄*(Pc_base+E(ê))。
%  Gauge : solve_KI_bar_gain 以 K_bar(1,1)=5/6 -> ^Bg_I=ghat_I_B (mT/A)、K_bar (6x6)。
%  Current : I = 1 A = FEM 激發電流（per fit-current-matches-sim rule）。
%
%  管線（driver 只 call 函式，一步一支）：
%    ① 載參數（ℓ̄ 初值、電流、取點半徑、匝數）        ② 沒疊代 bias 就 ê=0
%    ③ build_S_matrix（每點 S_matrix，fitting/solve 共用）
%    ④ fitting            -- 變數投影疊代解 ℓ̄（+ê），回傳 (ℓ̄, ê, Pc, J)
%    ⑤ solve_KI_bar_gain  -- 解 K̄_I、G、^Bĝ_I（+ calc_ellipsoid 出 𝒞/κ）
%    ⑥ save .mat + emit_model_results -> results/{single_param|eighteen_param}/*.pdf

clear; clc;

%% ---- ① 參數 --------------------------------------------------------------
USE_BIAS    = true;            % false = fix (在軸); true = no_fix (18-param bias)
MODE        = 'single';        % 'single' | 'sweep'
R_single_um = 150;             % single-mode 取點球半徑 [µm]（= PDF 報告半徑）
R_start_um  = 50;              % sweep start [µm]
R_step_um   = 5;               % sweep step  [µm]
R_end_um    = 500;             % sweep end   [µm]
I_actual    = 1;               % 驅動電流 [A] = FEM 激發（1 A）
TURNS       = 50;              % 線圈匝數（記錄用；FEM 場已含匝數，per-Ampere 模型不再引用）
SHAPE       = 'ball';          % 取點區域：ball ||p|| <= R about WP
ell0        = 0.5e-3;          % ℓ̄ 初值 [m]（= ℓ_design；SI，尺度佳）
dataset     = 'all';           % standard-mesh dataset
VARIANT     = 'gap_calibrate'; % FEM 變體子夾（gap_calibrate = 逐極單盤氣隙校準）

%% ---- paths ---------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\current_base'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');       % ansys_path
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');                  % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','main_function'));                                     % 全部模型函式（唯一函式夾）
model = 'long2016_hexapole_halfcut';

%% ---- 常數 + 慣例 ---------------------------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];                 % coil k excites this paper pole

%% ---- 載 FEM 一次（6 coils，actuator 框）+ ② ê 由 USE_BIAS 決定 ----------
fprintf('loading 6 coils (variant ''%s'', dataset ''%s'', USE_BIAS=%d) ...\n', VARIANT, dataset, USE_BIAS);
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, dataset, VARIANT);
vtag = ''; if ~strcmp(VARIANT,'standard'), vtag = ['_' VARIANT]; end
tag  = 'fixl'; if USE_BIAS, tag = 'bias'; end            % .mat 名稱前綴：fix -> fixl、bias -> bias

%% ---- 依 MODE 決定半徑清單 ------------------------------------------------
switch lower(MODE)
    case 'single', R_um_list = R_single_um;
    case 'sweep',  R_um_list = R_start_um:R_step_um:R_end_um;
    otherwise,     error('MODE must be ''single'' or ''sweep''.');
end

%% ---- 逐 R 擬合 + 存 .mat -------------------------------------------------
cal_dir = fullfile(TREE,'data'); if ~exist(cal_dir,'dir'); mkdir(cal_dir); end
for R_um = R_um_list
    [P, Bstack, npts] = select_ball(D, R_um*1e-6);                     % ③ P [m]、Bstack [mT]（S_matrix 由 build_S_matrix 提供）
    [ell, e_hat, Pc, J]  = fitting(P, Bstack, D.Pc_base, ell0, USE_BIAS);   % ④ 變數投影，回 (ℓ̄,ê,Pc,J)
    [K_bar, ghat_I_B, G] = solve_KI_bar_gain(ell, Pc, P, Bstack, D.F);      % ⑤ K̄_I、^Bg_I=ghat_I_B [mT/A]、G [mT]
    ell               = ell * 1e6;                                     % m → µm（此後 print/save 用 µm）
    gB = ghat_I_B;  Khat = K_bar;                                      % alias（.mat field 名沿用 gB/Khat）
    % 控制範圍性能量（singular_value.pdf）：C_mean=mean ∏σ_k、kappa_mean=mean σ₃/σ₁（球內真實節點）
    rm = calc_ellipsoid(P, gB*Khat, ell*1e-6, Pc);
    C_mean = rm.C_mean;  kappa_mean = rm.kappa_mean;  C_min = rm.C_min;  kappa_worst = rm.kappa_worst;  Np_range = rm.Np;
    save(fullfile(cal_dir, sprintf('fit_%s_R%03dum%s.mat', tag, R_um, vtag)), ...
         'ell','gB','Khat','G','e_hat','Pc','J','R_um','I_actual','TURNS','SHAPE','VARIANT','USE_BIAS', ...
         'C_mean','kappa_mean','C_min','kappa_worst','Np_range');
    fprintf('R=%3d um | npts=%6d | ell=%.2f um | ^Bg_I=%.4e mT/A | J=%.4e | C_mean=%.4g (mT/A)^3 | kappa_mean=%.4f (Np=%d)\n', ...
            R_um, npts, ell, ghat_I_B, J, C_mean, kappa_mean, Np_range);
end
fprintf('done (%s mode, USE_BIAS=%d, variant=%s): %d 個 .mat 存到 %s\n', ...
        MODE, USE_BIAS, VARIANT, numel(R_um_list), cal_dir);

%% ---- ⑥ 打包 PDF -> results/{single_param|eighteen_param}/ ---------------
emit_model_results(USE_BIAS, VARIANT, R_single_um);
