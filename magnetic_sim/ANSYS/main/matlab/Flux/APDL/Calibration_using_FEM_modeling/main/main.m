%% main.m -- 共用校正 driver（current + voltage）
%  cfg=model_config → extract_ansys_data(純讀) → build_actuator_data(轉 actuator+all-source) →
%  select_ball → fitting → [current: solve_current(F)] / [voltage: build_V_matrix→solve_voltage(V)] →
%  存 .mat(結果+設定) → emit_results 出 PDF。改頂部 per-run 調參即可切模型/變體/base/半徑。

clear; clc;

%% ---- per-run 調參（不進 config）---------------------------------------------
MODEL    = 'long2016_hexapole_halfcut';
GEOM     = 'tip40um';       % config 幾何變體：long2016 用 tip40um|tip400um；hung/NTU 用 ''（flat config）
VARIANT  = '';              % '' = 用該 geom 的 default_variant（data 變體）
DATASET  = 'all';
BASE     = 'voltage';       % 'current' | 'voltage'
USE_BIAS = true;            % e 開關：false=fix(single)、true=18-param(eighteen)
R_select = 150e-6;          % 取點球半徑 [m]
l0       = 0.5e-3;          % l_hat 初值 [m]
I_actual = 1;               % 驅動電流 [A]（= FEM 激發）
% 通用化旗標（'' → 由 cfg 提供預設；特例幾何自動吃自己的設定，不用手動改）
INTERP_TO = '';             % '' 正常；否則把本 variant 場內插到此參考 geom 的 R≤R_select 點雲（公平比較）
V_METHOD  = '';             % '' → cfg.v_method（'csv-tet' 預設 / 'scattered'）
% voltage-only 取樣調參
SOFF_upper = 4.572e-3;      % 上極 sensor 沿錐面距極尖 [m]
n_uniform  = 500;           % 每 sensor 圓柱撒點數（[MODIFIED 2026-08-06] 使用者拍板統一 500；
                            %   取樣標準誤 ≈ (σ/μ)/√n ≈ 0.11%。舊結果用 1e4（0.024%））
sensor_r   = 0.15e-3;       % sensor 圓柱半徑 [m]
axial_tol  = 0.10e-3;       % sensor 圓柱高（沿 n+）[m]

%% ---- paths（自身相對定位）---------------------------------------------------
here  = fileparts(mfilename('fullpath'));                    % .../Calibration_using_FEM_modeling/main
CAL   = fileparts(here);                                     % .../Calibration_using_FEM_modeling
APDL  = fileparts(CAL);                                      % .../main/matlab/APDL
ANSYS = fileparts(fileparts(fileparts(APDL)));               % .../ANSYS
addpath(fullfile(CAL, 'function'));                                     % pipeline + import_ansys_data/filter_iron_nodes/interp_field_to_points
addpath(fullfile(CAL, 'common_path'));                                  % ansys_path（共用路徑 resolver）
addpath(fullfile(CAL, 'utils'));                                          % [ADDED 2026-08-08] pole_sensor_geometry（sensor 幾何唯一來源）

%% ---- pipeline 前段（共用）---------------------------------------------------
cfg = model_config(MODEL, GEOM);
if isempty(VARIANT),   VARIANT   = cfg.default_variant;              end
if isempty(INTERP_TO), INTERP_TO = getdef(cfg, 'interp_to', '');     end
if isempty(V_METHOD),  V_METHOD  = getdef(cfg, 'v_method', 'csv-tet'); end

raw = extract_ansys_data(cfg, DATASET, VARIANT);            % 純讀 .dat（measure, T）
if isempty(INTERP_TO)
    ad = build_actuator_data(raw, cfg);                     % → actuator frame, mT, all-source
    [P, Bstack, npts] = cfg.select_ball(ad, R_select);
    Pc_base = ad.Pc_base;
else
    % 特例：把本 variant 場內插到「參考 geom」的點雲（同數量同分布，公平比較，如 tip400→tip40）
    cfg_ref = model_config(MODEL, INTERP_TO);
    raw_ref = extract_ansys_data(cfg_ref, DATASET, cfg_ref.default_variant);
    ad_ref  = build_actuator_data(raw_ref, cfg_ref);
    [P, ~, npts] = cfg_ref.select_ball(ad_ref, R_select);   % 參考點雲
    Bstack  = interp_field_to_points(cfg, VARIANT, cfg.R_act, P, getdef(cfg,'r_loc',0.6e-3));
    Pc_base = cfg.Pc_base;                                   % 本 variant 的電荷格（非 canonical）
end
[e, l_hat, J] = fitting(P, Bstack, Pc_base, l0, USE_BIAS);

%% ---- 共用結果紀錄（rec = 結果 + 設定，自描述）------------------------------
rec = struct('base',BASE, 'l_hat',l_hat, 'e',e, 'J',J, 'npts',npts, ...
             'model',MODEL, 'GEOM',GEOM, 'VARIANT',VARIANT, 'DATASET',DATASET, 'USE_BIAS',USE_BIAS, ...
             'INTERP_TO',INTERP_TO, 'V_METHOD',V_METHOD, ...
             'R_select',R_select, 'l0',l0, 'I_actual',I_actual, ...
             'apdl_to_paper_idx',cfg.apdl_to_paper_idx, 's_source',cfg.s_source);

%% ---- 後段（base 分兩路）-----------------------------------------------------
switch BASE
    case 'current'
        F = zeros(6, cfg.N_I);
        for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
        [KI_bar, gI_hat, G, rm] = solve_current(l_hat, e, Pc_base, P, Bstack, F);
        rec.KI_bar = KI_bar;  rec.gI_hat = gI_hat;  rec.G = G;  rec.Fmap = F;
    case 'voltage'
        [V, ~] = build_V_matrix(cfg, VARIANT, raw, cfg.S_hall, SOFF_upper, n_uniform, sensor_r, axial_tol, [], V_METHOD);
        [D_bar, gV_hat, G, rm] = solve_voltage(l_hat, e, Pc_base, P, Bstack, V);
        rec.D_bar = D_bar;  rec.gV_hat = gV_hat;  rec.G = G;  rec.V = V;
        % [ADDED 2026-08-06] sensor 取樣參數記進 .mat（自描述；與 Maxwell 分支同步）
        rec.SOFF_upper = SOFF_upper;  rec.n_uniform = n_uniform;
        rec.sensor_r = sensor_r;  rec.axial_tol = axial_tol;
    otherwise
        error('BASE 必為 ''current'' | ''voltage''');
end
rec.C_mean = rm.C_mean;  rec.kappa_mean = rm.kappa_mean;  rec.C_min = rm.C_min;  rec.kappa_worst = rm.kappa_worst;
if isfield(rm,'RMSPE'), rec.RMSPE = rm.RMSPE; end
if isfield(rm,'NMAE'),  rec.NMAE  = rm.NMAE;  end     % [ADDED 2026-08-15] 擬合 NMAE [%]（PDF 印這個）

%% ---- 存 .mat（自描述）------------------------------------------------------
tag = 'single'; if USE_BIAS, tag = 'eighteen'; end
matdir = fullfile(CAL, 'data', MODEL, '.mat');
if ~exist(matdir, 'dir'), mkdir(matdir); end
matfile = fullfile(matdir, sprintf('calib_%s_%s_R%03d_%s.mat', BASE, VARIANT, round(R_select*1e6), tag));
save(matfile, '-struct', 'rec');
fprintf('saved %s\n', matfile);
if strcmp(BASE,'current')
    fprintf('[current] l_hat=%.1f um  gI_hat=%.4e mT/A  K(1,1)=%.4f  C_mean=%.4g  kappa=%.4f\n', ...
            rec.l_hat*1e6, rec.gI_hat, rec.KI_bar(1,1), rec.C_mean, rec.kappa_mean);
else
    fprintf('[voltage] l_hat=%.1f um  gV_hat=%.4e mT/mV  D(1,1)=%.4f  C_mean=%.4g  kappa=%.4f\n', ...
            rec.l_hat*1e6, rec.gV_hat, rec.D_bar(1,1), rec.C_mean, rec.kappa_mean);
end

%% ---- 輸出 PDF ---------------------------------------------------------------
emit_results(matfile);

%% ---- local ------------------------------------------------------------------
function v = getdef(s, f, d)   % cfg 有欄位 f 且非空就用，否則預設 d
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
