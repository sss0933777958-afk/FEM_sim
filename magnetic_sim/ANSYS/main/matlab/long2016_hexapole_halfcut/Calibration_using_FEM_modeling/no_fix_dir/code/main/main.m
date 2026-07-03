%% main.m -- hexapole point-charge model (no_fix_l, 18-param bias) calibration driver
%  Document : "no_fix_l.pdf"
%  Model    : charges at pc = ell*(Pc_base + E(e_hat)), actuator frame; 18 params
%             = characteristic length ell + 1x17 bias e (e6z constrained).
%             Per simulation the 6 charge magnitudes g_j are profiled out (LS) -> G=D^v;
%             gauge with K_bar(1,1)=5/6 -> ^Bg_I = (6/5)h11, K_bar = (5/(6 g11))H_I  (paper step 8).
%  Current  : I = 1 A = FEM excitation current (per fit-current-matches-sim rule).
%
%  MODE switch:
%    'single' -> one fit at R = R_single_um
%    'sweep'  -> fits at R = R_start_um : R_step_um : R_end_um
%  Output: per-R fit .mat to data\ + console summary (PDF via code\main_function\emit_model_results.m).
%  All model math lives in code\main_function\ ; this file is just the driver.
%
%  Pipeline (top-to-bottom call order):
%    1) load_coils_actuator -- load 6-coil FEM, rotate to actuator frame  (load data)
%    2) select_ball         -- keep nodes inside the sampling ball R       (pick region)
%    3) fit_bias            -- lsqnonlin fit {ell, e_hat(17)}              (fit)
%    4) gauge_KI            -- profile g_j, gauge -> ^Bg_I, K_bar          (gauge)
%    5) region_field_err    -- relative RMS field error over region        (accuracy)
%    6) save fit .mat       -- console summary; PDF via main_function/emit_model_results.m

clear; clc;

%% ---- config ----------------------------------------------------------------
MODE        = 'single';        % 'single' | 'sweep'
R_single_um = 150;             % single-mode sampling-ball radius [um]
R_start_um  = 50;              % sweep start [um]
R_step_um   = 5;               % sweep step  [um]
R_end_um    = 500;             % sweep end   [um]
I_actual    = 1;               % drive current [A] = FEM excitation (1 A)
SHAPE       = 'ball';          % sampling region: ball ||p|| <= R about WP
ell0        = 0.5e-3;         % ell initial guess [m] (= ell_design; fit_bias 在 SI、well-scaled)
dataset     = 'all';           % standard-mesh dataset
VARIANT     = 'gap200um_mueq'; % [MODIFIED] FEM 變體子夾（'standard' = baseline；'gap200um_mueq' = gap200 2 段式 μ_eff）

%% ---- paths -----------------------------------------------------------------
TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
        'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\no_fix_dir'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');         % ansys_path
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants/import_ansys_data/filter_iron_nodes
addpath(fullfile(TREE,'code','function'));                                          % model helpers (last -> precedence)
addpath(fullfile(TREE,'code','main_function'));                                          % model helpers (last -> precedence)
model   = 'long2016_hexapole_halfcut';

%% ---- constants + conventions -----------------------------------------------
cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];                 % coil k excites this paper pole
coil_sign = [1 -1 1 -1 -1 1];                           % all-source display flip (upper P2/P4/P5)

%% ---- load FEM once (all 6 coils, actuator frame) ---------------------------
fprintf('loading 6 coils (variant ''%s'', dataset ''%s'') ...\n', VARIANT, dataset);   % [MODIFIED]
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, dataset, VARIANT);              % [MODIFIED] variant
vtag = '';                                                                              % [ADDED] output suffix
if ~strcmp(VARIANT,'standard'), vtag = ['_' VARIANT]; end

%% ---- radius list by mode ---------------------------------------------------
switch lower(MODE)
    case 'single', R_um_list = R_single_um;
    case 'sweep',  R_um_list = R_start_um:R_step_um:R_end_um;
    otherwise,     error('MODE must be ''single'' or ''sweep''.');
end

%% ---- fit + emit LaTeX result per R -----------------------------------------
for R_um = R_um_list
    [P, Bstack, npts] = select_ball(D, R_um*1e-6);                  % P [m]、Bstack [mT]
    [ell, e_hat, J]   = fit_bias(P, Bstack, D.Pc_base, ell0);       % 在 SI 公尺擬合（well-scaled）→ ell [m]
    Pc                = make_Pc(e_hat, D.Pc_base);
    [K_bar, ghat_I_B] = gauge_KI(ell, Pc, P, Bstack, D.F);          % K_bar、^Bg_I=ghat_I_B [mT/A]（ell/P 均公尺）
    ell               = ell * 1e6;                                 % m → µm（此後 write/print 用 µm）
    errpct            = region_field_err(Bstack, J);
    % PDF 輸出已分離到 code/main_function/emit_model_results.m（功能分開：main 只算+存 .mat + console）
    gB = ghat_I_B;  Khat = K_bar;                                   % alias（.mat field 名沿用 gB/Khat）
    % [ADDED] 控制範圍（R≤R_um 球）性能量（singular_value.pdf）：C_mean=mean ∏σ_k、kappa_mean=mean σ₃/σ₁
    %   P=actuator 節點、Pc=actuator bias 電荷（同框直接配對，不需 R_act）；C/kappa 框無關。
    rm = calc_range_metrics(P, gB*Khat, ell*1e-6, Pc);
    C_mean = rm.C_mean;  kappa_mean = rm.kappa_mean;  C_min = rm.C_min;  kappa_worst = rm.kappa_worst;  Np_range = rm.Np;
    cal_dir = fullfile(TREE,'data'); if ~exist(cal_dir,'dir'); mkdir(cal_dir); end
    save(fullfile(cal_dir, sprintf('fit_bias_R%03dum%s.mat', R_um, vtag)), ...
         'ell','gB','Khat','e_hat','J','errpct','R_um','I_actual','SHAPE','VARIANT', ...
         'C_mean','kappa_mean','C_min','kappa_worst','Np_range');            % [MODIFIED] + C_mean/kappa_mean
    fprintf('R=%3d um | npts=%6d | ell=%.2f µm | ^Bg_I=%.4e mT/A | err=%.2f%% | C_mean=%.4g (mT/A)^3 | kappa_mean=%.4f (Np=%d, C_min=%.4g, kappa_worst=%.4f)\n', ...
            R_um, npts, ell, ghat_I_B, errpct, C_mean, kappa_mean, Np_range, C_min, kappa_worst);
end
fprintf('done (%s mode, variant=%s): %d 個 .mat 存到 %s（PDF 由 code/main_function/emit_model_results.m 產生）\n', MODE, VARIANT, numel(R_um_list), fullfile(TREE,'data'));
