function plot_gain_iso_hist(USE_BIAS, R_FIT, R_EVAL, force, BINW, CMP, NTK, ZBASE)
% plot_gain_iso_hist -- 兩個六極設計的控制指標分布疊圖（Long Fei vs Zhi-Peng）
% =========================================================================
%   [MODIFIED 2026-08-21 使用者要求] **校正半徑與評估半徑分離**：
%     校正 = R <= R_FIT（預設 150 µm）的 **N_c 收斂點設計**（等測度網格降取樣）
%     評估 = R <= R_EVAL（預設 500 µm）的 **全部真實 .fld 格點**
%   （舊版兩者都用同一個 R_um。）兩個模型用同一組 (R_FIT, R_EVAL) → 可直接互比。
%
%     圖 1  Actuation volume 的立方根  C^(1/3)  [mT/A]
%           C = ∏ σ_k（σ = svd(S(p)·ᴮĤ_I) 的奇異值）→ 開三次方後回到 mT/A，
%           與 ĝ_I 同量綱、可直接比較兩個設計的「等效致動強度」。
%     圖 2  Isotropy  κ = σ₃/σ₁  [無因次]
%
%   ⚠ C 與 κ 是**點電荷模型導出的解析函數**，不是 raw FEM 場 —— 這裡的「節點」只用來
%     決定「在哪些位置求值」，並未讀取該點的 FEM 場值。
%   ⚠ **R_EVAL > R_FIT 時外圈是外推**：模型只在 R <= R_FIT 內被資料約束過。
%
%   ⚠ 兩個模型現在**工作區尺度相同**：Long Fei R_norm = 500 µm、
%     Zhi-Peng **R500** 變體 R_norm = 500 µm → R_EVAL = 500 對兩者都正好是工作球。
%     （舊版用 zhi_peng R594，兩者尺度不同、需另行換算，已不適用。）
%
%   ⚠ Zhi-Peng 的 'maxwell_split' 場（2026-08-20 改導磁係數後重解）六極已趨於等強，
%     K̄_I 非對角 30/30 全負、對角極差 12% —— 舊版檔頭寫的「±x 兩根強 2.3 倍、
%     8/30 非對角為正」是**舊場**（'maxwell'）的特性，已不適用。
%
%   資料量：兩個模型的六個 .fld 各 ~1.15 GB（201³ 格點），讀一輪約 1.6 分 →
%   結果存快取 data/gain_iso_hist_<model>_<variant>_fit<R>_eval<R>_<tag>.mat，第二次起秒級。
%
%   風格：①粗體框圖 + 疊圖直方圖慣例（nb=180、**共用 edges → 兩組 bin 寬相同**、
%   百分比縱軸）。[MODIFIED 2026-08-21] 照 err_hist 家族的定案：**不畫 mean 虛線、
%   長條不描黑邊、圖例只列系列名**；mean / CV / min / max 一律印在 console。
%   輸出 → figures/paper_fig/Section4_C/{gain_cbrt,iso}_hist_maxwell_<tag>.png
%
%   [ADDED 2026-08-24] 第 6 個引數 CMP 選「比哪兩組」（其餘流程完全共用）：
%     'design'（預設）Design A = long2016 / Design B = zhi_peng R500 maxwell_split
%                     -> 輸出檔名與行為與先前**逐字相同**
%     'gap'           No gap   = zhi_peng R500 maxwell_split（無氣隙）
%                     Have gap = zhi_peng R500 maxwell_gap（100 um 氣隙重解）
%                     -> {gain_cbrt,iso}_hist_gap_<tag>_R<eval>.png
%     顏色（2026-08-24 使用者要求，由原本的相反配色互換而來）：
%       No gap   -> 深藍 [0.05 0.10 0.95]（值高、落右）
%       Have gap -> 紅   [0.85 0.10 0.10]（值低、落左）
%     ⚠ 連帶：'design' 圖裡 zhi_peng maxwell_split = Design B 是**紅**，本圖它是
%       **藍** -> 同一份場在兩張圖不同色，跨圖對色時要注意。'ab_gap' 的 Have gap
%       則與本圖一致（都紅）。
% =========================================================================
    clc;
    if nargin < 1 || isempty(USE_BIAS), USE_BIAS = false; end
    if nargin < 2 || isempty(R_FIT),    R_FIT    = 150;   end   % 校正半徑 [µm]
    if nargin < 3 || isempty(R_EVAL),   R_EVAL   = 500;   end   % 評估半徑 [µm]
    if nargin < 4 || isempty(force),    force    = false; end
    % [ADDED 2026-08-21] BINW = 直接指定 C^(1/3) 圖的 bin 寬 [mT/A]；[] = 沿用 nb=180。
    %   起因：兩個設計相距 23.1 mT/A（是 Design A 自身寬度的 39 倍），共用 edges 之下
    %   180 根裡有 157 根落在中間的空白 -> bin 寬被空白決定，A 的整個分布只佔 4 根、
    %   看起來像一根尖刺。給 BINW 可把格子切細（例 0.05 -> A 約 12 根）。
    %   ⚠ 只套用在 C^(1/3) 那張；κ 的值域是 [0,1]、尺度完全不同，仍用 nb=180。
    if nargin < 5, BINW = []; end
    % [MODIFIED 2026-08-21] BINW 可給 1x2 = [gain_bin, kappa_bin]，分別套到兩張圖；
    %   給純量 = 只套 C^(1/3)（舊行為）；[] = 兩張都用 nb=180。
    %   預設值由 Freedman-Diaconis 決定（h = 2*IQR*n^(-1/3)）：四組資料一致指向
    %   「每個分布切 16~20 根」。C 的 A/B 寬度差 4.6 倍、FD 各要 0.033/0.165，
    %   共用 edges 只能取一個 -> 折衷 0.05（A 12 根、B 雜訊 17%%）。
    %   kappa 兩組寬度幾乎相同、FD 給 0.0106/0.0109 -> 取 0.011（各 19 根、雜訊 10%%）。
    %   ⚠ 舊值 nb=180 對 kappa 是**過度解析 3.5 倍**（每根僅 26 點、雜訊 20%%），
    %     圖上的鋸齒是統計噪音而非結構。
    if isempty(BINW),      BW_G = 0.05;   BW_K = 0.011;
    elseif isscalar(BINW), BW_G = BINW;   BW_K = 0.011;
    else,                  BW_G = BINW(1); BW_K = BINW(2);
    end

    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section4_C');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    % [MODIFIED 2026-08-28] Tree moved under matlab/Flux/; the old hardcoded path is gone,
    %   so addpath silently failed and an APDL copy of model_config/solve_* could shadow it.
    MAINR = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main';
    CAL   = fullfile(MAINR,'matlab','Flux','Maxwell');
    addpath(fullfile(CAL,'function'), fullfile(CAL,'utils'), fullfile(CAL,'common_path'));
    % [MODIFIED 2026-08-21] 移除 addpath(utils/long2016_hexapole_halfcut)：utils/ 已扁平化，
    %   conv_design_ws 在 function/（上一行已涵蓋），舊路徑不存在。

    tag = 'single';  if USE_BIAS, tag = 'eighteen'; end
    % {model, geom, variant, 圖例名, 顏色}
    %   [MODIFIED 2026-08-21] ① zhi_peng 改用 **R500** 幾何 + **maxwell_split** 場
    %   （2026-08-20 改導磁係數後重解的版本）；② 加 variant 欄位 —— 不明給的話
    %   conv_design_ws 會退回 cfg.default_variant，**靜默**用舊場。
    %   顏色照「顏色 = 模型」慣例：Long Fei 深藍、Zhi-Peng 紅。
    %   [MODIFIED 2026-08-21 使用者拍板] 圖例名改成 **Design A / Design B**（原 Long Fei /
    %   Zhi-Peng）；console 仍印 model 名，追溯不受影響。
    if nargin < 6 || isempty(CMP), CMP = 'design'; end
    % [ADDED 2026-09-02 使用者指定] NTK = [gain 圖根數, iso 圖根數]，兩軸同值。
    %   原本 axis_odd / ylim_odd 都寫死 3。使用者對不同圖指定不同根數
    %   （gap 的 gain = 4、iso = 3），故參數化。預設 [3 3] -> 既有圖重跑輸出不變。
    if nargin < 7 || isempty(NTK), NTK = [3 3]; end
    if isscalar(NTK), NTK = [NTK NTK]; end
    % [ADDED 2026-09-02 使用者指定] ZBASE = 另外出一版**水平軸從 0 起**的 Design A：
    %     gain  0~15 mT/A（只有 A 的 C^(1/3)~14.6 放得進；B 39.5 / Bgap 24.0 會被裁光）
    %     iso   0~1  （kappa 的物理滿量程）
    %   目的是「以絕對尺度看分布在全量程中的位置」，與現行的放大版並存，
    %   故檔名加後綴 **_z0**（zero-based），不覆蓋原檔。僅對 CMP='solo' 生效、只出 Design A。
    if nargin < 8 || isempty(ZBASE), ZBASE = false; end
    switch lower(CMP)
        case 'design'
            MD = { 'long2016_hexapole_halfcut', 'tip40um', '',              'Design A',  [0.05 0.10 0.95];
                   'zhi_peng',                  'R500',    'maxwell_split', 'Design B',  [0.85 0.10 0.10] };
            gstem = sprintf('gain_cbrt_hist_maxwell_%s_R%d', tag, R_EVAL);
            kstem = sprintf('iso_hist_maxwell_%s_R%d',       tag, R_EVAL);
        case 'gap'
            %   [MODIFIED 2026-08-24 使用者要求] 兩組顏色**互換**：No gap 深藍、Have gap 紅。
            MD = { 'zhi_peng', 'R500', 'maxwell_split', 'No gap',   [0.05 0.10 0.95];
                   'zhi_peng', 'R500', 'maxwell_gap',   'Have gap', [0.85 0.10 0.10] };
            gstem = sprintf('gain_cbrt_hist_gap_%s_R%d', tag, R_EVAL);
            kstem = sprintf('iso_hist_gap_%s_R%d',       tag, R_EVAL);
        case 'ab_gap'
            % [ADDED 2026-08-24] 同 'design'，但 Design B 換成**帶 100 um 氣隙**的重解
            %   （zhi_peng R500 maxwell_gap）。圖例仍用 Design A / Design B —— 這是
            %   同一組跨設計比較，只是 B 用了氣隙版；顏色也沿用 design 的藍/紅。
            MD = { 'long2016_hexapole_halfcut', 'tip40um', '',            'Design A',  [0.05 0.10 0.95];
                   'zhi_peng',                  'R500',    'maxwell_gap', 'Design B (Have gap)',  [0.85 0.10 0.10] };
            gstem = sprintf('gain_cbrt_hist_abgap_%s_R%d', tag, R_EVAL);
            kstem = sprintf('iso_hist_abgap_%s_R%d',       tag, R_EVAL);
        case 'v2'
            % [ADDED 2026-08-29] 志鵬自身的改版比較：Design B = R500（maxwell_split）vs Design C = R500_V2（maxwell_v2）。
            %   兩者極尖幾何相同（CAD 實測），差別在外圍（bbox 26.5->20 mm、薄舌片 15->10 mm）。
            %   配色沿用 'gap' 的慣例：值高者深藍、值低者紅 -> Design C 的 kappa 高 -> Design C 藍。
            MD = { 'zhi_peng', 'R500', 'maxwell_split', 'Design B', [0.85 0.10 0.10];
                   'zhi_peng', 'R500', 'maxwell_v2',    'Design C', [0.05 0.10 0.95] };
            gstem = sprintf('gain_cbrt_hist_v2_%s_R%d', tag, R_EVAL);
            kstem = sprintf('iso_hist_v2_%s_R%d',       tag, R_EVAL);
        case 'solo'
            % [ADDED 2026-08-31 使用者要求] **每個模型各自出圖**（不疊圖）：
            %   3 個模型 x {C^(1/3), kappa} = 6 張。顏色依**物理量**固定，不再依模型分色：
            %     capacity C^(1/3) -> 深藍 [0.05 0.10 0.95]（figure-style「直方圖用深藍」）
            %     isotropy  kappa  -> 紅   [0.85 0.10 0.10]
            %   單一系列不受「共用 bin 寬」硬條件約束（沒有比較對象）-> 每張各自用
            %   Freedman-Diaconis 選 bin 寬（見 fd_bin）；nb=180 對 ~1.8k 點是過度解析。
            MD = { 'long2016_hexapole_halfcut', 'tip40um', '',              'Design A',            [0.05 0.10 0.95];
                   'zhi_peng',                  'R500',    'maxwell_split', 'Design B',            [0.05 0.10 0.95];
                   'zhi_peng',                  'R500',    'maxwell_gap',   'Design B (Have gap)', [0.05 0.10 0.95] };
            % [MODIFIED 2026-09-01 使用者要求] C^(1/3) 三張的 bin 改細：FD 給的
            %   0.025/0.2/0.1（24/14/17 根）太粗，尤其 Design B 只有 14 根。改成
            %   **指定目標根數、各自反解 nice 寬度**（三張值域差 0.59/2.68/1.69，
            %   不能共用同一個寬度）。30 根 -> 每根仍有 ~50-66 點，形狀不會抖。
            %   kappa 三張維持 fd_bin 自動（使用者未要求更動）。
            NBAR_G = 30;
            % [ADDED 2026-09-02 使用者要求] 釘死每張的 bin 寬（[] = 走 fd_bin 自動）。
            %   Design A 釘 0.02 是為了與 sigma_hist_R150_maxwell **共用同一把尺**
            %   （兩張同為 mT/A、同在 R=150 評估，長條高度才能互比）。
            %   0.02 恰好也是 fd_bin(nbar=30) 當前算出的值 -> 本張外觀不變；
            %   釘死是為了防 fd_bin 隨資料更新而漂掉。
            %   ⚠ 要改 BW_SOLO_G{1} 必須同步改 plot_sigma_hist.m 的 BINW 預設。
            BW_SOLO_G = {0.02, [], []};
            % [ADDED 2026-09-02 使用者指定] kappa 三張也釘 **0.02**（原本 fd_bin 自動給 0.01）。
            %   三個設計是要互相比較的 -> bin 寬不同就不能比高度（面積 = 100%% x bin 寬），
            %   故三張一起改，不只改 Design A。kappa 無因次 -> 0.02 也是無因次。
            BW_SOLO_K = {0.02, 0.02, 0.02};            STAG  = {'A','B','Bgap'};
            gstem = '';   kstem = '';
        otherwise
            error('plot_gain_iso_hist:CMP', 'CMP 必為 design | gap | ab_gap | v2 | solo（給了 %s）', CMP);
    end
    SOLO = strcmpi(CMP,'solo');

    nM = size(MD,1);   D = cell(1,nM);
    for a = 1:nM
        vt = MD{a,3};  if isempty(vt), vt = 'default'; end
        cf = fullfile(here, 'data', sprintf('gain_iso_hist_%s_%s_fit%03d_eval%03d_%s.mat', ...
                                            MD{a,1}, vt, R_FIT, R_EVAL, tag));
        if exist(cf,'file') && ~force
            D{a} = load(cf);   fprintf('由快取載入 %s\n', cf);
        else
            D{a} = compute_one(MD{a,1}, MD{a,2}, MD{a,3}, R_FIT, R_EVAL, USE_BIAS, here, CAL);
            S = D{a};   save(cf, '-struct', 'S');
            fprintf('已存 %s\n', cf);
        end
    end

    % ---- console 統計（使用者要求：mean 印在 console，不進圖）----
    fprintf('%s\n', repmat('=',1,78));
    fprintf('校正 R <= %d um（N_c 降取樣）｜評估 R <= %d um（真實 .fld 格點）｜%s\n', ...
            R_FIT, R_EVAL, tag);
    fprintf('%s\n', repmat('-',1,78));
    for a = 1:nM
        S = D{a};
        fprintf('[%-9s] %s / %s：N_c=%d (%d,%d,%d)  評估 %d 點  l_hat=%.1f um  g_I=%.4f mT/A  NMAE=%.2f%%\n', ...
                MD{a,4}, S.model, S.variant, S.Nc, S.tri, S.npts, S.l_hat*1e6, S.gI, S.NMAE);
        fprintf('    C^(1/3)  mean = %8.4f mT/A   CV = %5.2f%%   min = %8.4f   max = %8.4f\n', ...
                mean(S.Ccbrt), std(S.Ccbrt)/mean(S.Ccbrt)*100, min(S.Ccbrt), max(S.Ccbrt));
        fprintf('    kappa    mean = %8.4f         CV = %5.2f%%   min = %8.4f   max = %8.4f\n', ...
                mean(S.kap), std(S.kap)/mean(S.kap)*100, min(S.kap), max(S.kap));
    end
    fprintf('%s\n', repmat('-',1,78));
    if nM == 2                                   % [MODIFIED 2026-08-31] 比值只在成對模式有意義
        fprintf('mean 比值（%s / %s）：C^(1/3) %.3f 倍｜kappa %.3f 倍\n', MD{2,4}, MD{1,4}, ...
                mean(D{2}.Ccbrt)/mean(D{1}.Ccbrt), mean(D{2}.kap)/mean(D{1}.kap));
    end
    fprintf('%s\n', repmat('=',1,78));

    % ---- solo：每個模型各自一張 C^(1/3)（藍）+ 一張 kappa（紅）----
    if SOLO
        COL_C = [0.05 0.10 0.95];                % capacity 深藍
        COL_K = [0.85 0.10 0.10];                % iso 紅
        % [ADDED 2026-08-31] 明給 x 範圍與刻度。xlim_pick 對這六組資料湊不出
        %   「乾淨端點 + 整數/0.N 刻度」—— kappa 三張更是直接掉進保底分支
        %   （端點 0.359966 / 0.584401、刻度 0.4098 這種四位小數，違反 figure-style）。
        %   下列六組都滿足：端點乾淨、刻度整數或 0.N、等距、**兩端留白 = 刻度間距**。
        %   kappa 為湊「0.N 刻度」各取 2 根（偶數）—— 規則明訂為了乾淨數字
        %   犧牲奇數是可接受的（優先序：乾淨數 > 奇數 > 等距 > 留白）。
        %   ⚠ 針對 **R_EVAL=150** 的資料範圍挑的；換評估半徑請重挑，否則走自動。
        % [MODIFIED 2026-09-02 使用者回報「kappa 圖偏一邊」] 實量發現實作已足夠：
        %   axis_odd 自己就能挑出這三組 gain 視野（與舊的寫死值逐位相同），
        %   且在 kappa 上能選到更平衡的落點 -> 改全部走自動，不再寫死。
        [XRG{1:3}] = deal([]);   [XTG{1:3}] = deal([]);
        % [MODIFIED 2026-09-02 使用者逐張指定端點] kappa 三張的水平軸全部寫死：
        %     A    起 0.7、訖 1.0（kappa = sigma3/sigma1 不可能 > 1）
        %     B    起 0.4、訖 0.7
        %     Bgap 起 0.3、訖 0.6
        %   三張 span 都是 **0.3** -> 水平尺度一致、形狀寬度可直接互比。
        %   ⚠ span 0.3 搭「一位小數刻度」只能給 **兩根**（間距 0.1、兩端留白也是 0.1）；
        %     要三根就得 s = 0.3/4 = 0.075 -> 0.775/0.85/0.925 三位小數，違反規則 9。
        %     使用者明確指定了端點 -> 以端點為優先，接受偶數根（同 axsh_gain 的既有例外）。
        %   其餘（gain 三張）仍走 axis_odd 自動。⚠ 刻度必須是 0.1 的倍數（規則 9）-> 間距至少 0.1
        %   -> 視野至少 0.4 寬，而 kappa 資料寬只有 ~0.21 -> **資料只能塔滿約一半視野**，
        %   這是規則 5+9 的硬限制，不是選錯落點；axis_odd 能做的是把它**擺中**。
        XRK = {[0.7 1.0], [0.4 0.7], [0.3 0.6]};
        XTK = {[0.8 0.9],  [0.5 0.6],  [0.4 0.5]};
        if R_EVAL ~= 150 || nM ~= 3                  % 非定案情境 -> 全部回到自動
            [XRG{:}] = deal([]);  [XTG{:}] = deal([]);
            [XRK{:}] = deal([]);  [XTK{:}] = deal([]);
        end
        % [ADDED 2026-09-02] ZBASE：另一版「從 0 起」的 Design A。兩軸五根。
        %   gain [0,15]：n=5 -> s = 15/6 = 2.5 -> 2.5/5/7.5/10/12.5，
        %     兩端留白 2.5 = 刻度間距、一位小數 -> 規則 5+9 全合。
        %   iso  [0,1] ：奇數五根且兩端留白 = 間距需 s = 1/6 = 0.1667（四位小數、違規則 9），
        %     故取 0.1/0.3/0.5/0.7/0.9（五根、等距 0.2、兩端留白 0.1 彼此相等但是間距的一半）。
        %     這是 [0,1] + 五根 + 一位小數三條下的唯一可行解，已向使用者說明。
        %   [MODIFIED 2026-09-02 使用者要求] **縱軸的根數與間距要跟非-z0 版相同** ->
        %     不再寫死 5，改用 NTK（預設 [3 3]）。水平軸因為明給 xt，不受 NTK 影響，
        %     仍是五根。同一筆資料 + 同一 NTK -> ylim_odd 給出與非-z0 版逐位相同的縱軸。
        if ZBASE && ~isempty(BINW)
            % [ADDED 2026-09-02 使用者要求] ZBASE + 明給 BINW = **只重出 gain 這張**、
            %   用指定的 bin 寬，水平軸維持 [0,15]。檔名再加 _bw<值>（小數點寫 p），
            %   不覆蓋前一版的 _z0（bin 0.02），兩版可直接對照。
            bwg_ = BINW(1);
            tagb_ = ['_bw' strrep(num2str(bwg_), '.', 'p')];
            render_overlay({D{1}.Ccbrt}, MD(1,4), {COL_C}, ...
                '$\mathbf{\mathcal{C}^{1/3}\;(mT/A)}$', ...
                fullfile(figdir, sprintf('gain_cbrt_hist_%s_%s_R%d_z0%s.png', STAG{1}, tag, R_EVAL, tagb_)), ...
                [0 15], bwg_, [2.5 5 7.5 10 12.5], NTK(1));
            return
        end
        if ZBASE
            render_overlay({D{1}.Ccbrt}, MD(1,4), {COL_C}, ...
                '$\mathbf{\mathcal{C}^{1/3}\;(mT/A)}$', ...
                fullfile(figdir, sprintf('gain_cbrt_hist_%s_%s_R%d_z0.png', STAG{1}, tag, R_EVAL)), ...
                [0 15], pick_bw(BW_SOLO_G{1}, D{1}.Ccbrt, NBAR_G), [2.5 5 7.5 10 12.5], NTK(1));
            render_overlay({D{1}.kap}, MD(1,4), {COL_K}, ...
                '$\mathbf{\kappa}$', ...
                fullfile(figdir, sprintf('iso_hist_%s_%s_R%d_z0.png', STAG{1}, tag, R_EVAL)), ...
                [0 1], pick_bw(BW_SOLO_K{1}, D{1}.kap), [0.1 0.3 0.5 0.7 0.9], NTK(2));
            return
        end
        for a = 1:nM
            render_overlay({D{a}.Ccbrt}, MD(a,4), {COL_C}, ...
                '$\mathbf{\mathcal{C}^{1/3}\;(mT/A)}$', ...
                fullfile(figdir, sprintf('gain_cbrt_hist_%s_%s_R%d.png', STAG{a}, tag, R_EVAL)), ...
                XRG{a}, pick_bw(BW_SOLO_G{a}, D{a}.Ccbrt, NBAR_G), XTG{a}, NTK(1));
            render_overlay({D{a}.kap}, MD(a,4), {COL_K}, ...
                '$\mathbf{\kappa}$', ...
                fullfile(figdir, sprintf('iso_hist_%s_%s_R%d.png', STAG{a}, tag, R_EVAL)), ...
                XRK{a}, pick_bw(BW_SOLO_K{a}, D{a}.kap), XTK{a}, NTK(2));
        end
        return
    end

    % ---- 兩張疊圖 ----
    % [ADDED 2026-08-21 使用者拍板] 水平軸範圍可明給（[] = 由 xlim_pick 自動選）。
    %   C^(1/3) 定成 **[10, 50]**：自動選出來的 [5,55] 兩端各留 s/2 = 5，右邊那 4.1 個
    %   單位是空的、看起來「右邊留太多」。改成 10~50 後端點就是刻度、完全不留白。
    %   ⚠ 代價：Design B 有一小段尾巴 > 50 會落在視野外（比例印在 console）。
    %   ⚠ [10,50] 是**針對 R_EVAL=500 的資料**挑的；換評估半徑資料範圍就不同（R<=150
    %     的 C^(1/3) 只落在 12.2~16.6 / 26.9~35.6），沿用會空掉一大半 → 其餘半徑走自動。
    XR_GAIN = [];   if R_EVAL == 500, XR_GAIN = [10 50]; end
    % [ADDED 2026-08-24 使用者拍板] gap 比較（資料 23.13~40.77）明給 [20,45] + 刻度
    %   25/30/35/40：整數、等距 5、且**兩端留白正好 = 刻度間距**（20->25、40->45 皆 5），
    %   正是 figure-style 對 2D 軸的理想。代價是 4 根（偶數）—— 規則明訂
    %   「為了湊整數而讓刻度變成 4 個是可以接受的」（優先序：整數 > 奇數 > 等距 > 留白）。
    %   為何要明給：自動選刻度的 xticks_in **只收奇數根**，在 [20,45] 內唯一的 5 根解是
    %   22.5/27.5/32.5/37.5/42.5（小數，違反「刻度一律整數」）；而 xlim_pick 若走自動
    %   會挑總留白最小的 [23,41]，端點 23 與刻度 24（40 與 41 同理）只差一格的 5.6%、
    %   數字擠成一團。⚠ 兩者都是**針對這組資料**挑的；換 R_EVAL 或換資料請拿掉走自動。
    XT_GAIN = [];
    if strcmpi(CMP,'gap') && R_EVAL == 150, XR_GAIN = [20 45];  XT_GAIN = [25 30 35 40]; end
    XR_ISO  = [];   XT_ISO = [];        % κ 的自動結果 [0,1] 已是最緊、不必 override
    % [ADDED 2026-08-29] 'v2' 的 κ 只落在 0.417~0.641，自動選會把端點訂在資料邊界
    %   (0.395007 / 0.663597) 且刻度變成四位小數 -> 違反 figure-style「端點必須是乾淨的數、
    %   刻度整數或 0.N」。明給 [0.3,0.7] + 0.4/0.5/0.6：3 根（奇數）、等距 0.1、
    %   兩端留白 0.1 = 刻度間距，端點 0.3/0.7 皆為 nice 值。
    if strcmpi(CMP,'v2'), XR_ISO = [0.3 0.7];  XT_ISO = [0.4 0.5 0.6];  end
    % [ADDED 2026-09-02 使用者指定「bin 寬盡量抓整數」] 'v2' 兩張的 bin 寬釘成整齊值：
    %     kappa    0.011 -> **0.02**（與其他所有 kappa 圖同寬，合併值域 0.224 -> 12 根）
    %     C^(1/3)  0.05  -> **0.1** （與 Design B 自己的 sigma / C 圖同寬）
    %   同時明給 C^(1/3) 的視野：B（37.96~40.64）與 C（29.44~31.50）相距很遠，
    %   自動選只能取 s=5 -> [25,45]（資料只填 56%%）；[28,44] + 32/36/40 同樣是
    %   三根等距、兩端留白 = 刻度間距（都是 4）、刻度皆整數，但填充率 70%%。
    if strcmpi(CMP,'v2') && R_EVAL == 150
        if isempty(BINW), BW_G = 0.1;   BW_K = 0.01;   end   % [2026-09-02] kappa 使用者改指定 0.01
        XR_GAIN = [28 44];   XT_GAIN = [32 36 40];
    end
    % [MODIFIED 2026-08-21] 檔名帶評估半徑 _R<eval>，讓不同 R_EVAL 的圖並存
    %   （原本沒帶，R150 版會直接蓋掉 R500 版）。
    render_overlay({D{1}.Ccbrt, D{2}.Ccbrt}, MD(:,4), MD(:,5), ...
        '$\mathbf{\mathcal{C}^{1/3}\;(mT/A)}$', ...
        fullfile(figdir, [gstem '.png']), XR_GAIN, BW_G, XT_GAIN, NTK(1));
    render_overlay({D{1}.kap,   D{2}.kap},   MD(:,4), MD(:,5), ...
        '$\mathbf{\kappa}$', ...
        fullfile(figdir, [kstem '.png']), XR_ISO, BW_K, XT_ISO, NTK(2));
end

% ============================================================================
function S = compute_one(model, geom, variant, R_FIT, R_EVAL, USE_BIAS, here, CAL)
% 校正用 **R<=R_FIT 的收斂點設計 N_c**（降取樣，等測度網格）；
% 評估（畫直方圖）用 **R<=R_EVAL 的全部真實 .fld 格點**。
    fprintf('--- %s：載入 + 校正（N_c 降取樣 @R<=%d um）---\n', model, R_FIT);
    cfg = model_config(model, geom);
    if isempty(variant), variant = cfg.default_variant; end
    raw = extract_maxwell_data(cfg, 'all', variant);
    ad  = build_actuator_data(raw, cfg);
    [P, ~, np] = cfg.select_ball(ad, R_EVAL*1e-6);    % 評估集（真實格點）

    F = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end

    % ---- 收斂點設計 N_c（per model）----
    %   zhi_peng 的六極不等強 → K̄_I 非對角恆有正值，故關掉那道閘（ki_gate=false），
    %   只保留「ℓ̂ 穩定 ∧ ĝ_I 穩定 ∧ 對角全正 ∧ 對角占優」。
    %   與 plot_svd_polar 同一把尺（使用者拍板 2026-08-15）：非 long2016 時 K̄_I 的
    %   物理結構條件不適用（六極不等強）→ **完全不納入判準**（ki_req=false），
    %   收斂點只由 ℓ̂ 與 ĝ_I 決定。
    is_l2016 = strcmp(model,'long2016_hexapole_halfcut');
    sg = struct('model',model, 'geom',geom, 'variant',variant, ...
                'ki_gate',is_l2016, 'ki_req',is_l2016);
    % ---- 讀收斂設計並就地校正（axes-shells 取樣器）------------------------
    %   [MODIFIED 2026-08-28] main.m 的 calib_*_convN*.mat 出自**舊的等測度取樣器**，
    %   與現行判準不一致（sample_axes_shells：六根致動軸 x Nr 層等距球殼 + 中心、
    %   N=6*Nr+1；穩態 = 連續 20 級步進 < 0.01%，收斂 = 連續 10 級落在穩態值 +/-0.2%）。
    %   改成讀 axsh 掃描的 N_c、在此重建設計就地擬合 —— 與 plot_svd_polar /
    %   plot_sigma_hist 同一把尺。快取優先序（variant 專屬 > R_FIT 專屬 > R-sweep）：
    %     ..._axsh<zsuf>_R<R>_<vtag>.mat  >  ..._axsh<zsuf>_R<R>.mat  >  ..._axsh<zsuf>.mat
    %   ⚠ split 與 gap 的 g_I 差約 34%，收斂點不保證相同 -> 不可共用一顆快取。
    zsuf_ = '';   if strcmp(model,'zhi_peng'), zsuf_ = '_zhi'; end
    vtag_ = regexprep(variant,'^maxwell_?','');
    dfv_ = fullfile(here,'data',sprintf('full_vs_conv_vs_R_maxwell_axsh%s_R%d_%s.mat',zsuf_,R_FIT,vtag_));
    dfx_ = fullfile(here,'data',sprintf('full_vs_conv_vs_R_maxwell_axsh%s_R%d.mat',zsuf_,R_FIT));
    df_  = fullfile(here,'data',['full_vs_conv_vs_R_maxwell_axsh' zsuf_ '.mat']);
    if     exist(dfv_,'file')==2, df_ = dfv_;
    elseif exist(dfx_,'file')==2, df_ = dfx_;   end
    assert(exist(df_,'file')==2, 'missing %s', df_);
    D_  = load(df_);
    NCv = D_.n_c1;   if USE_BIAS, NCv = D_.n_c2; end
    FBv = NCv(:).' == D_.npts_f(:).';        % 退回全格點的 R（校正集 = 評估集）要跳過
    Rv_ = D_.R_um(:).';   ok_ = isfinite(NCv(:).') & ~FBv;
    if ~any(ok_), ok_ = isfinite(NCv(:).'); end
    cd_ = find(ok_);   [~,jj_] = min(abs(Rv_(cd_) - R_FIT));   ix_ = cd_(jj_);
    Nc  = NCv(ix_);
    o_  = struct('model',model,'geom',geom,'variant',variant,'frame','actuator','quiet',true);
    Pd_ = [];   Bd_ = [];
    for Nr_ = ceil((Nc-1)/6) + (0:2)         % 濾鐵可能吃掉幾點 -> N 不一定剛好是 6*Nr+1
        Pq_ = conv_design_ws(Nr_, R_FIT*1e-6, struct('points_only',true,'R_act',cfg.R_act,'quiet',true));
        evalc('[Pk_,Bk_] = conv_design_ws([], R_FIT*1e-6, setfield(o_,''query'',Pq_));');
        if size(Pk_,1) == Nc,  Pd_ = Pk_;  Bd_ = Bk_;  break;  end
    end
    assert(~isempty(Pd_), 'no axes-shells Nr gives N=%d at R=%g um', Nc, R_FIT);
    [Pf_,Bf_]       = cfg.select_ball(ad, R_FIT*1e-6);   % out-of-sample NMAE 的評估集
    [e, l_hat]      = fitting(Pd_, Bd_, cfg.Pc_base, 0.5e-3, USE_BIAS);
    [KI, gI, ~, rm] = solve_current(l_hat, e, cfg.Pc_base, Pd_, Bd_, F, [], Pf_, Bf_);
    tri = [0 0 0];                           % 階梯是單一索引，不再是 (Nr,Nphi,Ntheta) 三元組
    fprintf(['  校正設計：axes-shells N_c=%d (Nr=%g) @R<=%d um（判準取自 R=%d）；' ...
             '評估集 = %d 個真實格點 @R<=%d um' newline], ...
            Nc, (Nc-1)/6, R_FIT, Rv_(ix_), np, R_EVAL);


    Pc   = make_Pc(e, cfg.Pc_base);
    Hhat = gI * KI;                          % ᴮĤ_I [mT/A]
    C    = zeros(np,1);   kap = zeros(np,1);
    for i = 1:np
        d  = P(i,:)/l_hat - Pc.';            % 6×3
        Sm = (d ./ (vecnorm(d,2,2).^3)).';   % 3×6
        sv = svd(Sm * Hhat);
        C(i) = prod(sv);   kap(i) = sv(3)/sv(1);
    end
    S = struct('model',model, 'geom',geom, 'variant',variant, ...
               'R_FIT',R_FIT, 'R_EVAL',R_EVAL, 'USE_BIAS',USE_BIAS, ...
               'npts',np, 'Nc',Nc, 'tri',tri, ...
               'l_hat',l_hat, 'gI',gI, 'KI',KI, 'NMAE',rm.NMAE, ...
               'C',C, 'Ccbrt',C.^(1/3), 'kap',kap);
end

% ============================================================================
function render_overlay(vals, names, cols, xlab, out, xrset, BINW, xtset, NTK)
% N 組資料的疊圖直方圖（各自正規化成百分比 → 比較的是分布形狀）
% [MODIFIED 2026-08-31] 由寫死 2 組泛化成 numel(vals) 組；N=2 的輸出逐位不變。
%   N=1（solo）時「共用 bin 寬」的硬條件自動失效（沒有比較對象）。
% [MODIFIED 2026-08-21 使用者拍板，比照 err_hist 家族]：
%   ① 長條**不描黑邊**（EdgeColor 'none'）——細長條逐根描邊會糊成一團黑
%   ② **不畫 mean 虛線**；③ 圖例**只列系列名**（mean / CV 改印 console）
%   ④ 縱軸換成 figure-style 2026-08-20 定案：內部刻度標數字、起點與終點都不標、
%      [0,T] 平分 N+1 段、步長為 0.1 的倍數且數字為整數或 0.N
%   ⑤ 兩組**共用同一組 edges → bin 寬必然相同**（figure-style 硬條件：
%      面積 = 100% × bin 寬，寬度不同就不能比高度）
    if nargin < 6, xrset = []; end
    if nargin < 7, BINW  = []; end
    % [ADDED 2026-08-24] xtset = 明給刻度（[] = 由 xticks_in 自動選）。需要它是因為
    %   xticks_in 只收奇數根，湊不出「整數 + 偶數根」這種合法組合（見呼叫端註解）。
    if nargin < 8, xtset = []; end
    % [MODIFIED 2026-09-02 使用者指定] 套上 figure-style 的九條規則：規則1 刻度 60、
    %   規則2 圖例 45、規則3 框線加粗、規則6 畫布等邊、規則7 圖例框與座標框同粗。
    %   原本 FS=28 同時餵刻度與軸標題、圖例寫死 24 → 拆成 FS / FSLAB / FSLEG。
    %   軸標題 FSLAB 規則未指定值 → 沿用本腳本原值 36。
    if nargin < 9 || isempty(NTK), NTK = 3; end
    FS = 60;   FSLAB = 36;   FSLEG = 45;   LWBOX = 5.0;   CANV = 14.5;
    ALPH = 0.60;   nb = 180;
    nS   = numel(vals);
    allv = cell2mat(cellfun(@(v) v(:), vals(:), 'UniformOutput', false));
    if isempty(BINW)
        edg = linspace(min(allv), max(allv), nb+1);
    else
        % [ADDED 2026-08-21] 明給 bin 寬：自 min 起等寬鋪到蓋過 max。仍是兩組共用
        %   同一組 edges -> 面積 = 100%% x BINW 的可比性不變。
        edg = min(allv) : BINW : (min(allv) + ceil((max(allv)-min(allv))/BINW)*BINW);
        nb  = numel(edg) - 1;
    end
    % **兩組共用同一組 edges** → 兩組 bin 寬必然相同（面積可比的硬條件）
    ctr  = (edg(1:end-1) + edg(2:end))/2;
    fprintf('  bin：共用 %d 個、寬 %.4g（兩組同寬；由兩組合併後的 %.4g ~ %.4g 均分）\n', ...
            nb, edg(2)-edg(1), min(allv), max(allv));
    for q = 1:nS                                  % [ADDED] 各組實際佔幾根 / 峰值多高
        wv = vals{q}(:);   hc = histcounts(wv, edg);
        fprintf(['    %-10s 全寬 %.4g -> %.1f 根（非空 %d），峰值 %.1f%%' char(10)], ...
                names{q}, max(wv)-min(wv), (max(wv)-min(wv))/(edg(2)-edg(1)), ...
                sum(hc>0), max(hc)/numel(wv)*100);
    end

    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig);   hold(ax,'on');
    h = gobjects(1,nS);   pk = [];
    for a = 1:nS
        p = histcounts(vals{a}, edg) / numel(vals{a}) * 100;
        h(a) = bar(ax, ctr, p, 1, 'FaceColor',cols{a}, 'FaceAlpha',ALPH, 'EdgeColor','none');
        pk = [pk p]; %#ok<AGROW>
    end

    box(ax,'on');  grid(ax,'off');
    % [MODIFIED 2026-08-21 使用者拍板] tick 朝外（TickDir 'out'，同 plot_sigma_hist）
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out');
    if isempty(xrset)
        % [MODIFIED 2026-09-02] 自動選軸改用 axis_odd（規則 5+9：三根等距、兩端留白
        %   = 刻度間距、刻度最多一位小數），且在同分的候選裡挑**左右留白最平衡**者。
        %   舊的 xlim_pick 只看根數 3~7、不管平衡 -> kappa 圖會明顯偏一邊。
        [xr, xt] = axis_odd(min(allv), max(allv), NTK);
    else
        xr = xrset;                                  % 使用者指定範圍
        if isempty(xtset), xt = xticks_in(xr); else, xt = xtset; end
        for a = 1:nS                                 % 落在視野外的樣本比例（誠實回報）
            f = 100*mean(vals{a}(:) < xr(1) | vals{a}(:) > xr(2));
            if f > 0, fprintf('  ⚠ %s 有 %.2f%% 的點落在視野 [%g, %g] 之外\n', ...
                              names{a}, f, xr(1), xr(2)); end
        end
    end
    xlim(ax, xr);   set(ax,'XTick',xt);
    [yr, yt] = ylim_odd(max(pk), NTK);       % [MODIFIED 2026-09-02] 規則 4+5：等距、兩端間距相等
    ylim(ax, yr);   set(ax,'YTick', yt);   ytop = yr(2);

    % x 起訖：只標數字、不畫 tick mark（figure-style）
    for xv = xr
        text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    xlabel(ax, xlab, 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);

    % 圖例：框內右上角、一行一系列、只列系列名
    % [MODIFIED 2026-09-02 使用者指示] **單一系列（solo）不畫圖例** —— 只有一組資料時
    %   「Design A」那顆圖例是冗餘的，標在圖說即可。多系列（疊圖比較）仍須保留，
    %   否則分不出哪組是哪個設計。
    if nS > 1
    lg = legend(ax, h, names, 'Interpreter','tex', 'Location','northeast', 'NumColumns',1);
    lg.FontSize = FSLEG;  lg.FontWeight = 'bold';
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
    % [規則2] ItemTokenSize 預設 30 pt 不隨 FontSize 放大，45 pt 字級下色塊會貼到框線。
    %   與 axsh_gain / axsh_kfro 用同一組值，保持跨圖一致。
    lg.ItemTokenSize = [55 25];
    end
    ax.Toolbar.Visible = 'off';   hold(ax,'off');

    % [MODIFIED 2026-09-02] 規則 6：exportgraphics 會裁掉畫布白邊，正方形畫布匯出後
    %   就不等邊 → 改用 print + 明確的 PaperPosition。
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');                % 14.5 in x 200 dpi = 2900 px 見方
    fprintf('wrote %s\n', out);
    close(fig);
end

% ============================================================================
function [xr, xt] = axis_odd(lo, hi, NTK)
% [ADDED 2026-09-02] 水平軸：規則 5（等距 tick、與起點/終點的間距也相等）+ 規則 9
%   （刻度最多一位小數）+ 使用者指定的「三根」。與 plot_sigma_hist.m 內同名函式相同。
%   形式：xr = [x0, x0+4*s]、刻度 = x0 + (1:3)*s。
%   排序鍵：① 視野最小 ② 刻度是整數優先 ③ 根數接近 3 ④ **左右留白對稱**。
%   第④ 鍵就是 2026-09-02 修「kappa 圖偏一邊」的那一項。
    best = [inf inf inf inf];   xr = [floor(lo) ceil(hi)];   xt = [];
    for k = -3:4
        for c = [1 1.5 2 2.5 5]
            s = c*10^k;
            if abs(s*10 - round(s*10)) > 1e-9 && s < 1, continue; end   % s 也要一位小數
            for n = NTK               % [MODIFIED 2026-09-02] 根數由呼叫端給（預設 3）
                span = (n+1)*s;
                if span < (hi-lo) - 1e-9, continue; end
                for m = (floor(hi/s) - n) + (-3:1)
                    x0 = m*s;   x1 = x0 + span;
                    if x0 > lo + 1e-9 || x1 < hi - 1e-9, continue; end
                    t = x0 + (1:n)*s;
                    if any(abs(t*10 - round(t*10)) > 1e-9), continue; end
                    sc = [span, double(any(abs(t-round(t)) > 1e-9)), abs(n-NTK), ...
                          abs((lo-x0) - (x1-hi))];
                    if lexlt(sc, best), best = sc;   xr = [x0 x1];   xt = t; end
                end
            end
        end
    end
end

% ============================================================================
function [lim, tk] = ylim_odd(maxv, NTK)
% [ADDED 2026-09-02] 縱軸：規則 4（起點與終點都不標）+ 規則 5（等距、且與兩端間距也
%   相等）+ 使用者指定的「三根 tick」。上緣 = 4*s、刻度 = (1:3)*s，0 與上緣不進 tk。
%
% [MODIFIED 2026-09-02 使用者回報「上面留白太多」]
%   舊版的步長只能取 [1 1.5 2 2.5 3 4 5 ...] 這組 nice 值。三根下上緣被鎖成 4*s，
%   步長一跳上緣就跳 → sigma 圖峰值 4.2%% 需 s>=1.05，舊清單只能給 1.5 → 上緣 6、
%   填充率僅 70%%，上方空一大塊。
%   現在允許**兩位有效數字**的步長（1.1 / 1.7 / 2.3 …），取滿足上緣的最小者；
%   但若 nice 步長只多浪費 <= 20%%，就還是用 nice 的（讓 2/4/6、2.5/5/7.5 這種
%   好讀的刻度不會被 1.7/3.4/5.1 取代）。
%   實測：sigma 5.5%% → 1.5/3/4.5 上緣 6 改成 1.1/2.2/3.3 上緣 4.4（填充 70%%→95%%）；
%          gain A 6.3%% 仍是 2/4/6（tight=1.7、nice=2，2/1.7=1.18 <= 1.20）。
%   縱軸不受規則 9（一位小數）約束 —— 那條只綁水平軸。
    n = NTK;                                      % [MODIFIED 2026-09-02] 由呼叫端給
    smin = 1.02*maxv/(n+1);                       % 上緣至少 1.02*maxv（長條不貼框）
    k     = floor(log10(max(smin, realmin)));
    tight = ceil(smin/10^(k-1) - 1e-9) * 10^(k-1);  % 兩位有效數字的最小合格步長
    nice  = inf;
    for kk = k-1:k+1
        for c = [1 1.5 2 2.5 3 4 5 6 8]
            s = c*10^kk;
            if s >= smin - 1e-12 && s < nice, nice = s; end
        end
    end
    s = tight;
    if nice <= 1.20*tight, s = nice; end          % nice 只貴 <=20%% 就用 nice
    tk  = (1:n)*s;
    lim = [0, (n+1)*s];
end

% ============================================================================
function w = pick_bw(wfix, v, nbar)
% 釘死的 bin 寬優先（跨圖共用同一把尺用）；[] 才退回 fd_bin 自動選。
    if ~isempty(wfix), w = wfix;  else, w = fd_bin(v, nbar);  end
end

% ============================================================================
function w = fd_bin(v, nbar)
% Freedman-Diaconis bin 寬（h = 2*IQR*n^(-1/3)），四捨五入到 nice 值 [1 2 2.5 5]x10^k。
%   只給**單一系列**的圖用 —— 多系列必須共用 bin 寬（figure-style 硬條件）。
%   分位數自己算，不依賴 Statistics Toolbox 的 iqr/prctile。
%   [ADDED 2026-09-01] 給了 nbar 就改成**目標根數**模式：h = range/nbar，一樣四捨五入
%   到 nice 值。用於各圖值域差很大、但希望根數相近的場合（FD 是統計最適、不保證根數）。
    if nargin < 2, nbar = []; end
    v = v(:);   n = numel(v);   sv = sort(v);
    if ~isempty(nbar)
        h = (max(v) - min(v)) / nbar;
        k = floor(log10(h));   c = [1 2 2.5 5 10];
        [~,j] = min(abs(c*10^k - h));   w = c(j)*10^k;   return
    end
    q  = interp1(linspace(0,1,n).', sv, [0.25; 0.75]);
    h  = 2*(q(2)-q(1))*n^(-1/3);
    if ~isfinite(h) || h <= 0, w = (max(v)-min(v))/20;   return;   end
    k  = floor(log10(h));   c = [1 2 2.5 5 10];
    [~,j] = min(abs(c*10^k - h));   w = c(j)*10^k;
end

% ============================================================================
function [xr, xt] = xlim_pick(lo, hi)
% 橫軸（figure-style 2026-08-19/20/21）。硬條件，優先序由高到低：
%   ① 內部刻度**數量奇數**（3 或 5）  ② 刻度**等距**
%   ③ 每個刻度與兩個端點都是**乾淨的數**（整數或一位小數 0.N）
%   ④ 兩端留白盡量小（端點另以 text 標數字、不畫 tick）
% ⚠ 別走回頭路：舊版「先由資料定 x0=floor(lo/s)*s、x1=ceil(hi/s)*s，再取內部刻度」
%   常常無解 —— 它要求 (x1−x0)/s 剛好是 4 或 6，而 x0/x1 又被資料釘死。2026-08-21 實測
%   C^(1/3)（12.2~50.9）在那個寫法下完全無解，掉進 fallback 後端點變成 8.3426 / 54.7651、
%   刻度互相疊字。
% 現行做法 = 「**先定刻度、再由刻度反推 xlim**」：
%   刻度取 o + m*s（o ∈ {0, s/2} 兩種相位），xlim = [t_first − s/2, t_last + s/2]
%   → 兩端留白恆為 s/2（對稱、且是最小可能值）。挑「總留白最小」者，同分再挑
%   「單邊留白最大值最小」者（避免一邊貼很緊、另一邊空一大塊 —— 使用者反映的
%   「右邊留太多」正是這種：舊版 [10,60] 右邊空 9.1，新版 [5,55] 只空 4.1）。
%   ① 端點取 **s/2 或 s/4 的整數倍**（都保證乾淨），由資料 floor/ceil 而來 → 貼著資料
%   ② 刻度取 o + k*s（o ∈ {0, s/2} 兩種相位）落在 (x0,x1) **嚴格內部**者
%   ③ 端點各再試「往外推一格」一次 —— 只用 ①② 常常湊不出奇數根（實測 R500 的
%      12.2~50.9 在基礎版全是偶數 4 或 8 根），往外推可以翻轉奇偶
%   排序：**整數刻度優先** → 留白小 → 根數接近 5
% [MODIFIED 2026-08-21 使用者反映「兩邊留太寬」] 端點格加入 **s/4**：
%   只用 s/2 時，「整數刻度」會把端點鎖在 s/2 的倍數上 —— C^(1/3) R150（14.36~40.77）
%   因此只能取 [10,45]（兩邊各空 4.4 / 4.2）。加了 s/4 之後可取 [12.5, 42.5]，
%   刻度仍是整數 20/30/40，但留白縮到 1.86 / 1.73（少 58%）。
    rng_ = max(hi-lo, realmin);
    % [MODIFIED 2026-08-24] 步長候選補上 **3 與 4**，與 figure-style.md 的 nice 清單
    %   [1 2 2.5 3 4 5 10] 對齊。少了 4 的時候，23.13~40.77 這種資料湊不出
    %   「整數刻度 + 奇數根 + 貼緊資料」的解 -> 只能退到 [15,45] 配 20/30/40，
    %   左邊留白 8.13（是右邊的 1.9 倍，使用者反映「左邊留太多」）。
    %   有了 s=4 就有 [23,41] 配 24/28/32/36/40：5 根、整數、等距、留白 0.13/0.23。
    cand = [1 2 2.5 3 4 5 10];
    xr = [];   xt = [];   best = [inf inf inf];
    for k = (floor(log10(rng_))-2) : (floor(log10(rng_))+1)
        for c = cand
            s = c*10^k;
            for hq = [s/2, s/4]                            % 端點格（粗 / 細兩種都試）
                b0 = floor(lo/hq)*hq;    b1 = ceil(hi/hq)*hq;
                for o = [0, s/2]                           % 刻度相位
                    for ext = 1:4
                        x0 = b0;   x1 = b1;
                        % [MODIFIED 2026-08-24] 加 ext==4「兩端同時往外推」。原本只有單邊
                        %   (ext 2 或 3)，kappa 的 0.381~0.6439 因此無解 -> 掉進保底分支，
                        %   端點變成 0.354683 / 0.670153、刻度 0.4248 這種四位小數還疊字。
                        %   兩端同推才到得了唯一乾淨解 [0.3, 0.7]（刻度 0.4/0.5/0.6）。
                        %   ⚠ 兩端同推留白必然較大 -> 排序上本來就墊底，只有在其他候選
                        %     全被 isclean 刷掉時才會出線；既有可解的圖不受影響。
                        if     ext == 2, x0 = b0 - hq;
                        elseif ext == 3, x1 = b1 + hq;
                        elseif ext == 4, x0 = b0 - hq;   x1 = b1 + hq;
                        end
                        t = (ceil((x0-o)/s + 1e-9) : floor((x1-o)/s - 1e-9))*s + o;
                        t = t(t > x0+1e-9 & t < x1-1e-9);
                        n = numel(t);
                        if mod(n,2)~=1 || n<3 || n>5, continue; end
                        if ~(isclean(x0) && isclean(x1) && isclean(t)), continue; end
                        isInt = double(any(abs(t - round(t)) > 1e-9));   % 0 = 全整數（優先）
                        sc = [isInt, (lo-x0)+(x1-hi), abs(n-5)];
                        if lexlt(sc, best), best = sc;   xr = [x0 x1];   xt = t; end
                    end
                end
            end
        end
    end
    if ~isempty(xr), return; end
    s  = rng_/6;   xr = [lo-0.1*rng_, hi+0.1*rng_];   xt = lo + (1:5)*s;   % 保底
end

% ============================================================================
function tf = lexlt(a, b)      % 字典序比較（容忍浮點）
    for i = 1:numel(a)
        if a(i) < b(i) - 1e-9, tf = true;  return; end
        if a(i) > b(i) + 1e-9, tf = false; return; end
    end
    tf = false;
end

% ============================================================================
function xt = xticks_in(xr)
% 給定的 [x0 x1] 內挑內部刻度：等距、**奇數個**（3 或 5）、每個都是乾淨的數。
%   兩種相位都試（s 的整數倍 / 再偏半格），取最接近 5 根的。端點本身不放 tick
%   —— 端點值由呼叫端以 text 標，照 figure-style「端點只標數字、不畫 tick mark」。
    span = xr(2) - xr(1);
    % [MODIFIED 2026-08-24] 步長候選補上 **3 與 4**，與 figure-style.md 的 nice 清單
    %   [1 2 2.5 3 4 5 10] 對齊。少了 4 的時候，23.13~40.77 這種資料湊不出
    %   「整數刻度 + 奇數根 + 貼緊資料」的解 -> 只能退到 [15,45] 配 20/30/40，
    %   左邊留白 8.13（是右邊的 1.9 倍，使用者反映「左邊留太多」）。
    %   有了 s=4 就有 [23,41] 配 24/28/32/36/40：5 根、整數、等距、留白 0.13/0.23。
    cand = [1 2 2.5 3 4 5 10];
    xt = [];   bestd = inf;
    for k = (floor(log10(span))-2) : (floor(log10(span))+1)
        for c = cand
            s = c*10^k;
            for o = [0, s/2]
                t = (ceil((xr(1)-o)/s + 1e-9) : floor((xr(2)-o)/s - 1e-9))*s + o;
                t = t(t > xr(1)+1e-9 & t < xr(2)-1e-9);
                n = numel(t);
                if mod(n,2)~=1 || n<3 || n>5 || ~isclean(t), continue; end
                if abs(n-5) < bestd, bestd = abs(n-5);  xt = t; end
            end
        end
    end
    if isempty(xt), xt = xr(1) + (1:3)*(span/4); end          % 保底
end

% ============================================================================
function tf = isclean(v)
% 「乾淨的數」＝整數或一位小數 0.N（figure-style 2026-08-20：端點不可是 0.123 這種尾數）
    tf = all(abs(v*10 - round(v*10)) < 1e-9);
end

% ============================================================================
function [lim, tk] = ylim_from_zero(maxv, N)
% 縱軸（figure-style 2026-08-20 定案，與 plot_err_hist_shell 同一支）：
%   • 自 0 起；**終點（軸上緣）與起點 0 都不標數字、也不畫 tick**（端點值由水平軸負責）
%   • **內部刻度要標數字**，每個必須是整數或一位小數 0.N
%   • 「N 根 tick」= 不含端點的內部刻度數；[0,T] 平分 N+1 段 → tk=(1:N)*s、T=(N+1)*s
%   • s 必為 0.1 的倍數且 (N+1)*s >= 1.08*maxv；在 [smin, 1.15*smin] 內挑最漂亮的
%     （整數 > 0.5 的倍數 > 0.2 的倍數 > 其餘；同分取最小 s，填充率最高）
    if nargin < 2 || isempty(N), N = 4; end
    smin = 1.08*maxv/(N+1);
    k0 = max(1, ceil(smin/0.1 - 1e-9));   k1 = max(k0, ceil(1.15*smin/0.1));
    best = k0;   bs = -1;
    for k = k0:k1
        if     mod(k,10) == 0, sc = 3;               % 整數
        elseif mod(k,5)  == 0, sc = 2;               % 0.5 的倍數
        elseif mod(k,2)  == 0, sc = 1;               % 0.2 的倍數
        else,                  sc = 0;
        end
        if sc > bs, bs = sc;  best = k; end
    end
    s   = best*0.1;
    tk  = round((1:N)*s*10)/10;
    lim = [0 (N+1)*s];
end

% ============================================================================
function Pc = make_Pc(e17, Pc_base)
    if isempty(e17) || all(e17(:) == 0), Pc = Pc_base;  return; end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
