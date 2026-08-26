%% main.m -- Maxwell 校正 driver（current + voltage）
%  cfg=model_config → extract_maxwell_data(純讀 .fld) → build_actuator_data(轉 actuator+all-source) →
%  select_ball → fitting → [current: solve_current(F)] / [voltage: build_V_matrix→solve_voltage(V)] →
%  存 .mat(結果+設定) → emit_results 出 PDF。改頂部 per-run 調參即可切模型/變體/base/半徑。

clear; clc;

%% ---- per-run 調參（不進 config）---------------------------------------------
MODEL    = 'zhi_peng';
GEOM     = 'R500';          % config 幾何變體
VARIANT  = 'maxwell_split';  % 志鵬現行模型（磁極切開重解版）；default 'maxwell' 是已被取代的舊匯出
DATASET  = 'all';
BASE     = 'current';       % 'current' | 'voltage'
USE_BIAS = false;           % e 開關：false=fix(single)、true=18-param(eighteen)
R_select = 150e-6;          % 取點球半徑 [m]
l0       = 0.5e-3;          % l_hat 初值 [m]
I_actual = 1;               % 驅動電流 [A]（= FEM 激發）
% [ADDED 2026-08-10] K̄(2,2) 物理約束（規則 calibration-transfer-matrix-output.md 附則二）。
%   []      = 自由擬合（預設，行為與先前完全相同）
%   0.8340  = 依「上磁極全錐應強於半切下磁極」設定；只改 G(2,2)、ℓ̂/e/其餘元素不動，
%             輸出變體名自動加 _k22_0p8340（不覆蓋自由擬合版）。
%   ⚠ 一旦設定，**current 與 voltage 兩邊都要用同一個值重跑**（兩邊必須共用同一顆 G）。
K22_SET  = [];
% [MODIFIED 2026-08-23 使用者拍板] 工作空間的取樣設計改由**收斂迴圈自動決定**
%   （原本是手填三元組）。三種給法：
%     'auto'    跑收斂階梯自動決定 (N_r,N_phi,N_theta)  <- 定案做法
%     []        R 球內**全部** .fld 格點（R150 -> 1791 點；全取樣基準，不迴圈）
%     [a b c]   手動指定設計（重現舊結果用，不迴圈）
%   R=150um 既有收斂設計：long2016 single (3,8,22)=528 / eighteen (2,4,10)=80；
%   zhi_peng R500 maxwell_split single (1,3,8)=24 / eighteen (1,2,5)=10。
%   輸出檔名自動加 _convN<點數>，不覆蓋全格點版。
GRID_NRPT = [6 6 6];
% ---- 收斂判準（工作空間球）：後 KWIN 步「l_hat 與 g_I」變化率都 < TOL，且 K_I_bar 合物理
%   ⚠ 判準序列固定用 solve_current 的 [l_hat, g_I]（**兩個 base 都是**），與 2026-08-23
%     之前的 conv_design 同一把尺 —— 這樣既有記錄的 N_c 值才可比。voltage 的
%     solve_voltage 在收斂之後才跑一次。
CONV_SEED    = [1 2 3];     % 階梯種子
CONV_TOL     = 0.005;       % 變化率門檻（0.5%）
CONV_KWIN    = 10;          % 連續穩定步數
CONV_NDMAX   = 150;         % 最多掃幾級
CONV_KI_GATE = false;       % false = 放寬「K_I_bar 非對角全負」（六極不等強，如 zhi_peng）
CONV_KI_REQ  = true;        % false = K_I_bar 完全不參與判準，只看 l_hat + g_I
% 通用化旗標（'' → 由 cfg 提供預設；特例幾何自動吃自己的設定，不用手動改）
INTERP_TO = '';             % '' 正常；否則把本 variant 場內插到此參考 geom 的 R≤R_select 點雲（公平比較）
V_METHOD  = '';             % '' → cfg.v_method（現行 'grid' = 規則格三線性；另有 'scattered' / 'csv-tet'）
% voltage-only 取樣調參
SOFF_upper = 4.572e-3;        % 上極 sensor 沿錐面距極尖 [m]（定案值 4.572e-3）
SOFF_lower = 4.572e-3;        % [ADDED] 下極 sensor 沿錐面距極尖 [m]（原寫死 4.572e-3；兩層要一起移才與示意圖一致）
% [REMOVED 2026-08-23 使用者拍板] n_uniform **廢除** —— 每 sensor 的取樣設計不再手填，
%   改由「conv_design_sensor 產點 <-> build_V_matrix 組 V」的收斂迴圈自動決定（見下）。
%   舊定案手填值是等測度網格 (3,20,3)=180 點/sensor（比舊亂數 500 點準 11 倍）；
%   收斂迴圈實測會停在 (3,13,2)=78 點，下游 g_V 只差 0.009%、NMAE 完全相同。
SEN_SEED = [1 3 2];         % sensor 階梯種子（地板 N_theta>=3、N_z>=2 由 conv_design_sensor 內建）
SEN_TOL  = 0.001;           % 判準：V 的 36 個元素**全部**離密參考 < 0.1%
SEN_KWIN = 10;              % ...且要**連續 SEN_KWIN 級**都過（見下方迴圈註解）
SEN_REF  = [10 100 10];     % 密參考設計（10000 點/sensor；自身誤差 0.0042%）
sensor_r   = 0.15e-3;       % sensor 圓柱半徑 [m]
axial_tol  = 0.10e-3;       % sensor 圓柱高（沿 n+）[m]

%% ---- paths（自身相對定位）---------------------------------------------------
here  = fileparts(mfilename('fullpath'));                    % .../Maxwell/main
CAL   = fileparts(here);                                     % .../Maxwell
SOLV  = fileparts(CAL);                                      % .../main/matlab
ANSYS = fileparts(fileparts(fileparts(SOLV)));               % .../ANSYS
addpath(fullfile(CAL, 'function'));                                     % pipeline + import_maxwell_fld / filter_iron_nodes / conv_design_ws
addpath(fullfile(CAL, 'common_path'));                                  % 共用路徑 resolver（如有）
addpath(fullfile(CAL, 'utils'));                                          % [ADDED 2026-08-08] pole_sensor_geometry（sensor 幾何唯一來源）

%% ---- pipeline 前段（共用）---------------------------------------------------
cfg = model_config(MODEL, GEOM);
if isempty(VARIANT),   VARIANT   = cfg.default_variant;              end
if isempty(INTERP_TO), INTERP_TO = getdef(cfg, 'interp_to', '');     end
if isempty(V_METHOD),  V_METHOD  = getdef(cfg, 'v_method', 'csv-tet'); end

raw = extract_maxwell_data(cfg, DATASET, VARIANT);          % 純讀 .fld（Maxwell frame, T）
if ~isempty(INTERP_TO)
    % [REMOVED 2026-08-23 使用者拍板] 「內插到參考 geom 點雲」這條路已停用。
    %   它靠 function/interp_field_to_points.m，而那支呼叫的 extract_ansys_data
    %   **只存在於 APDL 分支**（Maxwell 樹下沒有）-> 在純 Maxwell path 上必定
    %   Undefined function。該檔已刪除。
    error('main:interpToRemoved', ...
          ['INTERP_TO=''%s'' 已停用：interp_field_to_points 於 2026-08-23 刪除' char(10) ...
           '（它相依 APDL 分支的 extract_ansys_data，在 Maxwell 樹下無法運作）。'], INTERP_TO);
end
ad      = build_actuator_data(raw, cfg);                    % → actuator frame, mT, all-source
Pc_base = ad.Pc_base;
% [ADDED 2026-08-26 使用者拍板] 誤差指標的**評估點雲** = R 內全部 .fld 格點。
%   solve_* 仍用取樣點 P 擬合出 G，但 NMAE / RMSPE 改在這組點上算 —— 也就是
%   「用 N 個取樣點校正出來的模型，拿去擬合全格點的場」（out-of-sample）。
%   ⚠ 舊行為是 in-sample（在取樣點自己身上算），對點數有系統性樂觀：實測 35 點的
%     in-sample NMAE 2.60% 比全格點擬合的 2.81% 還「低」，純粹因為點少好配；同一顆
%     模型拿去配全格點其實是 2.81%。故所有帶 _convN 的舊 PDF/.mat 的 NMAE 不可與新值並列。
%   ⚠ GRID_NRPT=[]（全格點取樣）時評估點 = 取樣點，新舊值逐位相同。
[P_ev, B_ev, npts_ev] = cfg.select_ball(ad, R_select);
fprintf('[eval] NMAE 評估點雲：R<=%g um 內全部格點 %d 點\n', R_select*1e6, npts_ev);
F = zeros(6, cfg.N_I);
for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

%% ---- ① sensor 取樣設計收斂（voltage only）----------------------------------
%  [ADDED 2026-08-23 使用者拍板] conv_design_sensor **只產點** -> build_V_matrix
%  內插組 V -> 這裡判「離密參考的偏差」。收斂後 V 就定住，不參與 ② 的迭代。
V = [];   sen_tri = [];   dev = NaN;
if strcmp(BASE, 'voltage')
    % Maxwell 的 sensor 場是**另一組匯出**（WP 細格框 ±0.6mm 涵蓋不到 WP 外 ~4.5mm 的
    % sensor）→ 另載 dataset='voltage' 的粗格，只餵 build_V_matrix；上面的電荷擬合
    % 仍用 WP 細格 raw。兩者同一個 Maxwell 座標框，可直接並用。
    raw_v = extract_maxwell_data(cfg, 'voltage', VARIANT);
    fprintf('[voltage] sensor 場格點 %d（WP 擬合場格點 %d）\n', numel(raw_v.x), numel(raw.x));
    bV = @(t) build_V_matrix(cfg, VARIANT, raw_v, cfg.S_hall, SOFF_upper, t, ...
                             sensor_r, axial_tol, [], V_METHOD, SOFF_lower);
    % ⚠ info 是**第 6 個**輸出（三線性併入後多了 B）；少一個 ~ 會拿到 tri 而非 info。
    [~,~,~,~,~,si] = conv_design_sensor(SEN_SEED(1), SEN_SEED(2), SEN_SEED(3), ...
        struct('ladder',CONV_NDMAX, 'sensor_r',sensor_r, 'axial_tol',axial_tol));
    Vref = bV(SEN_REF);
    fprintf('[sensor] 密參考 (%d,%d,%d) = %d 點/sensor\n', SEN_REF, prod(SEN_REF));
    % [MODIFIED 2026-08-23 使用者拍板] 判準由「第一次低於門檻就停」改成
    %   **連續 SEN_KWIN 級都低於門檻**，與 ws 那條的穩健度看齊（原本 ws 要連續 10 步、
    %   sensor 只要一次，兩者不對等；實測第一次過關是 0.0985% 對門檻 0.100%，
    %   只差 1.5%，靠單一級決定太脆弱）。
    %   ⚠ 這裡**維持 'ref' 判準**（離密參考的偏差），只是多要求連續 KWIN 級 ——
    %     **不可**改成 ws 那種「相對前一步的變化率」：實測 sensor 端那條會連兩次
    %     假收斂（5 點 / 8 點，實際誤差都是 0.6%），因為階梯在角向緩慢爬升時
    %     「變化率小」與「離真值近」是兩回事。
    %   收斂設計取**連續窗的第一級**（與 ws 的 q_hit 同慣例）。
    sen_hit = false;
    sen_ok  = false(1, size(si.ladder,1));
    sen_dev = nan(1, size(si.ladder,1));
    for q = 1:size(si.ladder,1)
        V   = bV(si.ladder(q,:));
        dev = max(abs(V(:) - Vref(:)) ./ abs(Vref(:)));      % 36 個元素全部要過
        sen_ok(q) = dev < SEN_TOL;   sen_dev(q) = dev;
        i0 = ftrue(sen_ok(1:q), SEN_KWIN);
        if ~isnan(i0)
            sen_tri = si.ladder(i0,:);
            V = bV(sen_tri);                                 % 回到窗首那級重算 V
            fprintf(['[sensor] 收斂設計 (%d,%d,%d)、%d 點/sensor' ...
                     '（掃了 %d 級；第 %d~%d 級連續 %d 次 < %.3f%%，' ...
                     '窗內偏差 %.4f~%.4f%%）\n'], ...
                    sen_tri, prod(sen_tri), q, i0, i0+SEN_KWIN-1, SEN_KWIN, ...
                    SEN_TOL*100, min(sen_dev(i0:i0+SEN_KWIN-1))*100, ...
                    max(sen_dev(i0:i0+SEN_KWIN-1))*100);
            sen_hit = true;   break
        end
    end
    assert(sen_hit, ['sensor：%d 級內沒有連續 %d 級 < %.3f%%' ...
                     '（最後一級 (%d,%d,%d) 偏差 %.4f%%）'], ...
           CONV_NDMAX, SEN_KWIN, SEN_TOL*100, si.ladder(end,:), sen_dev(end)*100);
end

%% ---- ②③④ 產點+取場 -> 校正 -> 判收斂（迴圈由 main 驅動）--------------------
%  [ADDED 2026-08-23 使用者拍板] conv_design_ws **只做「決定內插點位置 + 三線性取場」**；
%  校正（fitting / solve_*）與判斷都在這裡。
wsopt = struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, 'frame','actuator');
AUTO  = ischar(GRID_NRPT) || isstring(GRID_NRPT);
if AUTO
    assert(strcmpi(GRID_NRPT,'auto'), 'GRID_NRPT 必為 ''auto'' | [] | [N_r N_phi N_theta]');
    [~,~,~,wi] = conv_design_ws(CONV_SEED(1), CONV_SEED(2), CONV_SEED(3), R_select, ...
                                struct('ladder',CONV_NDMAX));
    LADW = wi.ladder;
else
    LADW = [];                                              % 不迴圈，只跑一次
end

nq  = 1;   if AUTO, nq = size(LADW,1); end
SER = nan(nq,2);   OKV = false(1,nq);   tri_ws = [];   ws_hit = ~AUTO;
for q = 1:nq
    % ② 產點 + 三線性取場
    if AUTO
        tri_ws = LADW(q,:);
        [P, Bstack, ~, gi] = conv_design_ws(tri_ws(1), tri_ws(2), tri_ws(3), R_select, wsopt);
        npts = gi.npts_kept;
    elseif isempty(GRID_NRPT)
        [P, Bstack, npts] = cfg.select_ball(ad, R_select);  % 全格點基準
    else
        tri_ws = GRID_NRPT(:).';
        [P, Bstack, ~, gi] = conv_design_ws(tri_ws(1), tri_ws(2), tri_ws(3), R_select, wsopt);
        npts = gi.npts_kept;
        fprintf('[grid] 手動指定設計 (%d,%d,%d) → %d 點\n', tri_ws, npts);
    end

    % ③ 校正（在 main）
    [e, l_hat, J] = fitting(P, Bstack, Pc_base, l0, USE_BIAS);
    [KI_bar, gI_hat, G, rm] = solve_current(l_hat, e, Pc_base, P, Bstack, F, K22_SET, P_ev, B_ev);
    if ~AUTO, break; end

    % ④ 判收斂：l_hat 與 g_I 都穩 ∧ K_I_bar 合物理
    SER(q,:) = [l_hat*1e6, gI_hat];
    od = KI_bar(~eye(6));   [~,am] = max(abs(KI_bar),[],2);
    OKV(q) = all(diag(KI_bar) > 0) && isequal(am(:).',1:6) && (~CONV_KI_GATE || all(od < 0));
    is = fstab_all(SER(1:q,:), CONV_TOL, CONV_KWIN);
    ik = 1;   if CONV_KI_REQ, ik = ftrue(OKV(1:q), CONV_KWIN); end
    if ~isnan(is) && ~isnan(ik)
        tri_ws = LADW(max([is ik]), :);                     % 收斂點是**判準起點**那一級
        fprintf('[ws] 收斂設計 (%d,%d,%d)、%d 點（掃了 %d 級）\n', ...
                tri_ws, prod(tri_ws), q);
        [P, Bstack, ~, gi] = conv_design_ws(tri_ws(1), tri_ws(2), tri_ws(3), R_select, wsopt);
        npts = gi.npts_kept;
        [e, l_hat, J] = fitting(P, Bstack, Pc_base, l0, USE_BIAS);
        [KI_bar, gI_hat, G, rm] = solve_current(l_hat, e, Pc_base, P, Bstack, F, K22_SET, P_ev, B_ev);
        ws_hit = true;   break
    end
end
% [FIXED 2026-08-24] 原本寫成 assert(ws_hit, ..., LADW(end,:))：MATLAB 會**先求值全部
%   引數**再呼叫 assert，所以 GRID_NRPT=[] / 手動三元組（不迴圈、LADW 為空）時
%   LADW(end,:) 必定丟「Index in position 1 is invalid」。等於 2026-08-23 改寫後
%   這兩條非 auto 的路從來沒被跑通過。改成只在真的沒收斂時才去索引 LADW。
if ~ws_hit
    error('main:wsNotConverged', ...
          'workspace：%d 級內未達判準（最後一級 (%d,%d,%d)）', CONV_NDMAX, LADW(end,:));
end

%% ---- 共用結果紀錄（rec = 結果 + 設定，自描述）------------------------------
rec = struct('base',BASE, 'l_hat',l_hat, 'e',e, 'J',J, 'npts',npts, ...
             'model',MODEL, 'GEOM',GEOM, 'VARIANT',VARIANT, 'DATASET',DATASET, 'USE_BIAS',USE_BIAS, ...
             'INTERP_TO',INTERP_TO, 'V_METHOD',V_METHOD, ...
             'R_select',R_select, 'l0',l0, 'I_actual',I_actual, ...
             'apdl_to_paper_idx',cfg.apdl_to_paper_idx, 's_source',cfg.s_source);

%% ---- 後段（base 分兩路）-----------------------------------------------------
switch BASE
    case 'current'
        % KI_bar / gI_hat / G / rm 已在 ③ 算出（判準與最終結果共用同一次 solve_current）
        rec.KI_bar = KI_bar;  rec.gI_hat = gI_hat;  rec.G = G;  rec.Fmap = F;
    case 'voltage'
        % V 已在 ① 由收斂迴圈定出（conv_design_sensor 產點 -> build_V_matrix 組 V）
        rec.SOFF_upper = SOFF_upper;  rec.SOFF_lower = SOFF_lower;   % emit_results 據此加檔名後綴
        % sensor 取樣參數一併記進 .mat（自描述）
        rec.sen_tri = sen_tri;  rec.sensor_r = sensor_r;  rec.axial_tol = axial_tol;
        rec.SEN_REF = SEN_REF;  rec.SEN_TOL = SEN_TOL;  rec.SEN_KWIN = SEN_KWIN;
        [D_bar, gV_hat, G, rm] = solve_voltage(l_hat, e, Pc_base, P, Bstack, V, K22_SET, P_ev, B_ev);
        rec.D_bar = D_bar;  rec.gV_hat = gV_hat;  rec.G = G;  rec.V = V;
    otherwise
        error('BASE 必為 ''current'' | ''voltage''');
end
rec.C_mean = rm.C_mean;  rec.kappa_mean = rm.kappa_mean;  rec.C_min = rm.C_min;  rec.kappa_worst = rm.kappa_worst;
if isfield(rm,'RMSPE'), rec.RMSPE = rm.RMSPE; end
if isfield(rm,'NMAE'),  rec.NMAE  = rm.NMAE;  end     % [ADDED 2026-08-15] 擬合 NMAE [%]（PDF 印這個）
% [ADDED 2026-08-26] NMAE 現在是 out-of-sample（在 P_ev 全格點上算）。in-sample 版與
%   評估點雲的描述一併存進 .mat，讓每顆結果都能自證是哪種算法。
if isfield(rm,'NMAE_in'), rec.NMAE_in = rm.NMAE_in; end
if isfield(rm,'NMAE_on'), rec.NMAE_on = rm.NMAE_on; end
if isfield(rm,'Np_eval'), rec.Np_eval = rm.Np_eval; end

% [ADDED 2026-08-10] K22 約束 → 輸出變體名加 tag（**資料載入仍用原 VARIANT**，只有輸出改名，
%   故自由擬合版 model_results_*_maxwell.pdf 不會被覆蓋）。
if ~isempty(K22_SET)
    rec.K22_set = rm.K22_set;  rec.K22_free = rm.K22_free;
    VAR_OUT     = sprintf('%s_k22_%s', VARIANT, strrep(sprintf('%.4f', K22_SET), '.', 'p'));
    rec.VARIANT = VAR_OUT;                             % emit_results 依此決定 PDF 檔名
else
    VAR_OUT = VARIANT;
end
% [ADDED 2026-08-13] 等測度網格取樣 → 檔名加 _convN<點數>，與全格點版並存
%   [MODIFIED 2026-08-23] tri_ws 空 = 全格點模式（不加 tag）；'auto' 與手動指定都加。
if ~isempty(tri_ws)
    VAR_OUT       = sprintf('%s_convN%d', VAR_OUT, npts);
    rec.VARIANT   = VAR_OUT;
    rec.GRID_NRPT = tri_ws;
    rec.conv_auto = AUTO;
end

%% ---- 存 .mat（自描述）------------------------------------------------------
tag = 'single'; if USE_BIAS, tag = 'eighteen'; end
matdir = fullfile(CAL, 'data', MODEL, '.mat');
if ~exist(matdir, 'dir'), mkdir(matdir); end
% voltage 且 sensor 距非定案 4.572mm → .mat 檔名加後綴，不覆蓋既有結果
%   ⚠ [2026-08-15] 不要把 soff tag 併進 VAR_OUT —— emit_results 會**自己**再加一次
%   （它由 rec.SOFF_upper 判斷），併進去會得到 `..._soff3mm_soff3mm.pdf`。踩過。
soffsfx = '';
if strcmp(BASE,'voltage') && (abs(SOFF_upper-4.572e-3) > 1e-9 || abs(SOFF_lower-4.572e-3) > 1e-9)
    soffsfx = sprintf('_soff%gmm', SOFF_upper*1e3);
end
matfile = fullfile(matdir, sprintf('calib_%s_%s%s_R%03d_%s.mat', BASE, VAR_OUT, soffsfx, round(R_select*1e6), tag));
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

function i0 = fstab_all(SER, tol, K)
% 第一個 i，使其後連續 K 步「每一個序列」相對前一步的變化率都 < tol。
    i0 = NaN;   if size(SER,1) < K+1, return; end
    rel = abs(diff(SER,1,1)) ./ abs(SER(1:end-1,:));       % (n-1) x nser
    mx  = max(rel, [], 2);                                  % 每步取最嚴的那個序列
    for i = 1:numel(mx)-K+1
        if all(mx(i:i+K-1) < tol), i0 = i;  return; end
    end
end

function i0 = ftrue(v, K)
% 第一個 i，使 v(i..i+K-1) 全為 true。
    i0 = NaN;
    for i = 1:numel(v)-K+1
        if all(v(i:i+K-1)), i0 = i;  return; end
    end
end
