function [P, Bstack, info] = conv_design_ws(Nr, R, opt)
% conv_design_ws -- 工作空間（球）：決定取樣點位置 + 三線性內插取場
% =========================================================================
%   [SPLIT  2026-08-23 使用者拍板] 由 conv_design.m 拆出（另一支 = conv_design_sensor）。
%   [MERGE  2026-08-23 使用者拍板] sphere_grid_sample.m 整支併進來（該檔已刪除）。
%   [REPLACE 2026-09-02 使用者拍板] **等測度球格（配比 1 : 3 : 3*pi）整套廢除**，
%     產點方法換成**六軸殼層**（原 sample_axes_shells.m，該檔整支併進來並刪除）。
%     使用者指示：「所有產點的功能請都在這個函式中」「介面改成通用一點的模式」。
%     取場引擎（濾鐵 + 規則格三線性內插 + frame 轉換）**逐字不動**。
%
%   ⚠ **本函式只做兩件事：①決定取樣點位置 ②三線性內插取場。**
%     校正（fitting / solve_*）與收斂判斷**一律在 main.m**，迴圈也由 main.m 驅動：
%
%       [~,~,wi] = conv_design_ws([], R, struct('ladder',30));   % 先取階梯表
%       for q = 1:numel(wi.ladder)
%           [P,Bs] = conv_design_ws(wi.ladder(q), R, o);         % 產點 + 取場
%           ... main 自己 fitting + solve_* ...                  % 校正
%           ... main 自己累積序列、判收斂，過了就 break ...        % 判斷
%       end
%
%   ⚠ 這是**引擎**（無對應圖），被 main.m 與多支繪圖腳本呼叫，勿當孤兒刪除。
%
%   ── 六步流程 ────────────────────────────────────────────────
%     1. 決定設計       Nr（殼層數）／直接給點／由 N_target 反解
%     2. 建查詢點       六軸 x Nr 層 + 中心 -> 轉 measure frame -> z + SPH_OFST
%     3. 濾鐵           filter_iron_nodes(x,y,z,cfg)，吃全域座標、零轉換
%     4. 定位           fx=(x-x0)/h -> i=floor(fx), tx=fx-i；取 8 個角點
%     5. 三線性內插     7 次 lerp（4 沿 x -> 2 沿 y -> 1 沿 z）-> B [T] -> x1e3 -> mT
%     6. 最後轉換       平移 z - SPH_OFST，再旋轉 cfg.R_act -> actuator frame
%
%   ── Step 1+2：六軸殼層（決定「點擺在哪」）────────────────────
%     六個致動軸方向（ACTUATOR frame，順序同 Pc_base 的 P1..P6）：
%         (+1,0,0) (-1,0,0)   <- u 軸，P1 / P2
%         (0,+1,0) (0,-1,0)   <- v 軸，P3 / P4
%         (0,0,+1) (0,0,-1)   <- w 軸，P5 / P6
%     殼層半徑   r_k = (R/Nr)*k ,  k = 1 ... Nr      （最外層正好落在 r = R）
%     點數       N = 6*Nr + 1                        （+1 是中心點）
%         Nr = 1 ->  7     Nr = 3 -> 19     Nr = 10 ->  61
%         Nr = 2 -> 13     Nr = 4 -> 25     Nr = 166 -> 997
%
%   ── 為什麼是這個佈局 ────────────────────────────────────────
%     **每顆極在每一層都有自己的點，且 +/- 兩側對等。** 這是 4 點設計做不到的：
%     只取兩根軸的 +/- 成對在數學上必然奇異（rank 5，未取樣那根軸的電荷對相對
%     取樣平面鏡像對稱）；只取三根軸的 + 側雖然滿秩，卻永遠不取樣 P2/P4/P6 ->
%     它們的電荷強度不可辨識、K_I_bar 失去物理結構。
%
%     徑向鋪點也讓 l_hat 的槓桿最大：沿極軸移動正是「到該電荷的距離」變化最快的
%     方向，而場的衰減是唯一攜帶長度尺度的東西，等於正面取樣它。（相對地，在中心
%     點 S 完全與 l_hat 無關。）
%
%   ── 已知弱點 ────────────────────────────────────────────────
%     所有點落在三條直線上，沒有離軸的角向覆蓋。電荷強度與 l_hat 會很乾淨，但 e
%     裡**垂直於各極軸**的分量槓桿很小（把電荷側向移動對它自己軸上的場是二階
%     效應）。預期 single 表現好，eighteen 有些分量會偏軟。
%
%   ── 用法 ────────────────────────────────────────────────────
%     [~,~,i]  = conv_design_ws([], R, struct('ladder',30))    % 只出階梯表
%     [P,Bs]   = conv_design_ws(3, 150e-6)                     % Nr=3 -> 19 點
%     [P,Bs]   = conv_design_ws(Pq, 150e-6, o)                 % 直接給點（任意點取場）
%     [P,Bs]   = conv_design_ws([], R, struct('N_target',61))  % 由目標點數反解 Nr
%
%   輸入
%     Nr     : **通用設計引數**，三種形式擇一
%                純量正整數  Nr（殼層數）-> 六軸殼層，N = 6*Nr + 1
%                Np x 3 矩陣 直接給查詢點（**measure frame** [m]），跳過 Step 1-2
%                []          由 opt 決定（opt.ladder 只出表；opt.N_target 反解 Nr）
%     R      : 取樣球半徑 [m]（直接給點時可省）
%     opt    : 選項 struct（皆可省略）
%       .ladder   M   **只**回 info.ladder（Mx1 的 Nr 階梯 = (1:M)'），不產點
%       .N_target N   目標點數 -> Nr = max(1, round((N-1)/6))
%       .center   T/F 是否含中心點（預設 true）
%       .points_only T/F **只產點、不取場**（預設 false）。跳過 Step 3-6，Bstack 回空。
%                 P 回 **measure frame**、info.P_act 回 actuator frame，與已退役的
%                 sample_axes_shells 契約完全相同（那支的 16 個呼叫端靠這個模式遷移）。
%                 ⚠ 同時給 opt.R_act 時**連 model_config 都不載入**（不需要 .fld）——
%                   Force 那條線只要殼層幾何、不要場，讀 2 GB .fld 是純浪費。
%       .R_act    3x3 **產點用**的致動軸基底（預設 cfg.R_act）
%                 ⚠ 只影響「軸往哪擺」，**不影響 Step 6 的輸出 frame**（那一律用
%                   cfg.R_act）。axsh_vs_R.m 的旋轉不變性測試就是餵一個轉過的
%                   R_act 給這裡、模型端維持 cfg.R_act，兩者必須分開。
%       .query    Np x 3  **相容別名**，等同把點放進 Nr（舊呼叫端用）
%       .frame    'actuator'（預設）| 'measure'   輸出座標與 B 的 frame
%       .drop_iron true/false  是否剔除落在鐵件內的點（預設 true）
%       .model .geom .variant  資料路由（預設 long2016_hexapole_halfcut/tip40um）
%       .raw      直接給已載入的 raw（跳過 extract_maxwell_data）
%       .quiet    true 則不印報告（逐級呼叫時用；預設 true）
%
%   輸出
%     P       Np x 3   取樣點座標 [m]（frame 依 opt.frame；'actuator' 時球心為原點）
%     Bstack  3Np x 6  磁場 [mT]，逐點 [bx;by;bz] 堆疊、第 2 維 = 6 個激發（paper 序）
%                      （已是 fitting / solve_* 要的形狀，main 拿了直接用）
%     info    struct   .Nr .radii .center .npts_design .npts_kept .ladder .R
%                      .n_iron .n_outbox .keep .N .frame .variant .B
%                      （.B = Np x 3 x N_I 的原始形狀，需要時可直接取用）
%
%   ⚠ 回傳點數可能少於設計點數：落在鐵件內、或落在 .fld 格盒外的點會被剔除
%     （用 info.keep 對回原本的產點順序）。R > 401 um 時磁極會伸進取樣球。
%   ⚠ 磁場是**內插**值（匯出格距 20 um），非 FEM 節點原值。內插不增加資訊。
% =========================================================================
    if nargin < 1, Nr = []; end
    if nargin < 2, R      = []; end
    if nargin < 3 || isempty(opt), opt = struct(); end
    gv = @(f,d) getdef_(opt, f, d);

    LADM  = gv('ladder',   []);
    PONLY = logical(gv('points_only', false));
    NTGT  = gv('N_target', []);
    CENTER    = logical(gv('center', true));
    FRAME     = lower(gv('frame', 'actuator'));
    DROP_IRON = gv('drop_iron', true);
    QUIET     = logical(gv('quiet', true));
    assert(any(strcmp(FRAME,{'actuator','measure'})), 'opt.frame 必為 actuator | measure');

    % ---- 拆解通用的 Nr 引數（純量 = 殼層數；矩陣 = 直接給點）----------------
    QRY = gv('query', []);                     % 相容別名（舊呼叫端）
    if ~isempty(Nr) && ~isscalar(Nr)
        assert(size(Nr,2) == 3, 'conv_design_ws:Nr', ...
               'Nr 必為純量殼層數，或 Np x 3 的查詢點（measure frame）');
        QRY = Nr;   Nr = [];
    elseif ~isempty(Nr)
        validateattributes(Nr, {'numeric'}, {'scalar','positive','integer'}, ...
                           'conv_design_ws', 'Nr');
    end

    % ---- ladder 模式：只出階梯表、不產點 -----------------------------------
    %   六軸殼層的階梯就是 Nr = 1,2,3,...（每個 Nr 都是一級，點數 6*Nr+1）。
    if ~isempty(LADM)
        validateattributes(LADM, {'numeric'}, {'scalar','positive','integer'});
        P = [];   Bstack = [];
        info = struct('Nr',[], 'radii',[], 'center',CENTER, 'npts_design',[], ...
                      'npts_kept',[], 'ladder',(1:LADM).', 'R',R, 'frame',FRAME);
        return
    end

    % ---- points_only + 明給 R_act：完全不需要模型，直接產點後回傳 ----------
    if PONLY && ~isempty(gv('R_act', []))
        [P, info] = axes_shells_(Nr, R, gv('R_act',[]), CENTER, NTGT);
        Bstack = [];   info.frame = 'measure';   return
    end

    solver_path();
    MODEL = gv('model', 'long2016_hexapole_halfcut');
    GEOM  = gv('geom',  '');
    if isempty(GEOM) && strcmp(MODEL,'long2016_hexapole_halfcut'), GEOM = 'tip40um'; end
    cfg     = model_config(MODEL, GEOM);
    VARIANT = gv('variant', cfg.default_variant);
    RAWIN   = gv('raw', []);
    RACT    = gv('R_act', cfg.R_act);           % **產點用**的基底（見檔頭警語）

    info = struct('Nr',[], 'radii',[], 'center',CENTER, 'npts_design',[], ...
                  'npts_kept',[], 'ladder',[], 'R',R, 'frame',FRAME, 'variant',VARIANT, ...
                  'n_iron',0, 'n_outbox',0, 'keep',[], 'N',[], 'B',[]);

    %% ==== Step 1+2：決定 Nr -> 建六軸殼層查詢點（直接建在全域座標）========
    if isempty(QRY)
        [Pm, gi_] = axes_shells_(Nr, R, RACT, CENTER, NTGT);
        if PONLY                                      % 只產點：跳過 Step 3-6
            P = Pm;   Bstack = [];   info = gi_;   info.frame = 'measure';   return
        end
        x = Pm(:,1);   y = Pm(:,2);   z = Pm(:,3) + cfg.SPH_OFST;   % **全域座標**
        info.Nr = gi_.Nr;   info.radii = gi_.radii;   info.npts_design = gi_.npts_design;
    else
        % 通用介面：直接給 measure frame 的查詢點
        assert(size(QRY,2) == 3, 'Nr / opt.query 必須是 Np x 3（measure frame [m]）');
        x = QRY(:,1);   y = QRY(:,2);   z = QRY(:,3) + cfg.SPH_OFST;
        info.npts_design = size(QRY,1);
    end
    info.N = numel(x);

    %% ==== Step 3：濾鐵（filter_iron_nodes 吃全域座標）====================
    keep = true(numel(x),1);
    if DROP_IRON
        keep = filter_iron_nodes(x, y, z, cfg);   keep = keep(:);
        info.n_iron = nnz(~keep);
    end

    %% ==== Step 4+5：定位 + 三線性內插 -> mT ==============================
    lat = lattice_cached(cfg, VARIANT, RAWIN);    % 規則格查表（persistent 快取）
    [Bt, inbox] = trilerp(lat, x(keep), y(keep), z(keep));       % [T]
    info.n_outbox = nnz(~inbox);
    idx = find(keep);   keep(idx(~inbox)) = false;               % 盒外一併剔除
    B = 1e3 * Bt;                                                % T -> mT

    %% ==== Step 6：平移 z - SPH_OFST，再旋轉 R_act ========================
    %   ⚠ 這裡一律用 **cfg.R_act**（模型的 frame），不是 opt.R_act（產點用的基底）。
    x = x(keep);   y = y(keep);   z = z(keep) - cfg.SPH_OFST;    % (1) 平移
    if strcmp(FRAME,'actuator')                                  % (2) 旋轉
        Pr = (cfg.R_act * [x, y, z].').';
        x = Pr(:,1);   y = Pr(:,2);   z = Pr(:,3);
        for m = 1:size(B,3)
            B(:,:,m) = (cfg.R_act * B(:,:,m).').';               % B 只旋轉
        end
    end

    P  = [x y z];   nkeep = size(P,1);
    Bstack = zeros(3*nkeep, size(B,3));                          % 逐點 [bx;by;bz]
    for j = 1:size(B,3), Bstack(:,j) = reshape(B(:,:,j).', [], 1); end

    info.keep = keep;   info.npts_kept = nkeep;   info.B = B;
    if ~QUIET, report(info); end
end

% ============================================================================
function [Pm, info] = axes_shells_(Nr, R, RACT, CENTER, NTGT)
% 六軸殼層產點 —— 實作**逐行對應** Algorithm 1「Generate Interpolation Points」。
%   Input  : N_a (= 本函式的 Nr，每根軸每個方向的等分數), R
%   Output : p_a  (N_tot x 3)，N_tot = 1 + 6*N_a
%
%   演算法原文（行號與下面的程式一一對應）：
%      1: d     <- R / N_a                       // division spacing
%      2: N_tot <- 1 + 6*N_a                     // total number of points
%      3: p_a   <- zeros(N_tot, 3)
%      4: E     <- eye(3)                        // unit vectors of three axes
%      5: cnt   <- 1;  p_a(cnt,:) <- [0,0,0]     // center point
%      6: for j = 1 to 3 do
%      7:     for s in {+1,-1} do                // +/- direction
%      8:         for i = 1 to N_a do
%      9:             cnt <- cnt + 1
%     10:             p_a(cnt,:) <- s*i*d*E(j,:) // i-th point on axis j
%     11:         end for
%     12:     end for
%     13: end for
%     14: return p_a
%
%   ⚠ **列的順序**是「軸 j -> 正負 s -> 等分 i」（+x 全部、-x 全部、+y…），
%     不是舊版的「逐層 -> 六個方向」。點集完全相同、擬合與順序無關，故數值不變
%     （已實測 main.m 的 l_hat / g_I 逐位相同）。
%
%   本函式回 p_a 轉到 **measure frame** 的 Pm；info.P_act 保留 actuator frame 的 p_a。
%   （原 sample_axes_shells.m，2026-09-02 整支併入本檔、該檔已刪除。）
    assert(~isempty(R) && isscalar(R) && R > 0, ...
           'conv_design_ws:noR', 'conv_design_ws 必須給取樣球半徑 R');
    if isempty(Nr)
        assert(~isempty(NTGT), 'conv_design_ws:noDesign', ...
               'Nr 為空時必須給 opt.N_target（或 opt.ladder）');
        Nr = max(1, round((NTGT - double(CENTER)) / 6));
    end

    Na   = Nr;                                      % 演算法的 N_a
    d    = R / Na;                                  % 1: division spacing
    Ntot = 1 + 6*Na;                                % 2: total number of points
    pa   = zeros(Ntot, 3);                          % 3
    E    = eye(3);                                  % 4: unit vectors of three axes
    cnt  = 1;   pa(cnt,:) = [0, 0, 0];              % 5: center point
    for j = 1:3                                     % 6
        for s = [+1, -1]                            % 7: +/- direction
            for i = 1:Na                            % 8
                cnt = cnt + 1;                      % 9
                pa(cnt,:) = s * i * d * E(j,:);     % 10: i-th point on axis j
            end                                     % 11
        end                                         % 12
    end                                             % 13

    % 選項（演算法沒有這一條）：CENTER=false 時把第 5 行的中心點拿掉。
    if ~CENTER, pa(1,:) = []; end

    % actuator -> measure：R_act 的**列**是 measure 座標下的致動軸單位向量，
    %   故 P_meas = R_act' * P_act；以列向量寫就是 pa * R_act。
    Pm = pa * RACT;                                 % 14: return p_a（轉到 measure frame）
    info = struct('Nr',Na, 'R',R, 'h',d, 'radii',(1:Na)*d, 'center',CENTER, ...
                  'npts_design',size(pa,1), 'npts',size(pa,1), 'P_act',pa, ...
                  'npts_kept',size(pa,1), 'ladder',[], 'n_iron',0, 'n_outbox',0, ...
                  'keep',true(size(pa,1),1), 'N',size(pa,1), 'B',[], 'frame','measure');
end

% ============================================================================
function [Bq, ok] = trilerp(lat, xq, yq, zq)
% 定位 + 三線性內插（7 次 lerp）。xq/yq/zq 為 Np x 1 **全域座標** [m]。
%   回 Bq（僅 in-box 的列，Np' x 3 x N_I，單位同 lat.V 即 Tesla）與 ok（是否在盒內）。
    nx = lat.n(1);   ny = lat.n(2);   nz = lat.n(3);

    % --- 定位：fx 的整數部分 = 第幾格（0 起算）、小數部分 = 格內正規化座標 ---
    fx = (xq - lat.o(1)) / lat.h(1);
    fy = (yq - lat.o(2)) / lat.h(2);
    fz = (zq - lat.o(3)) / lat.h(3);
    % ⚠ 容差不可省：落在**最外面那一面**的點（例如 .fld 自己的邊界格點）算出來會是
    %   60.0000000000001 這種值，無容差就會被誤判成盒外（實測 3000 點中誤殺 100 個）。
    %   TOL 以「格」為單位，1e-9 格 = 2e-14 m，物理上完全可忽略。
    TOL = 1e-9;
    ok = fx >= -TOL & fx <= nx-1+TOL & fy >= -TOL & fy <= ny-1+TOL ...
       & fz >= -TOL & fz <= nz-1+TOL;
    fx = fx(ok);   fy = fy(ok);   fz = fz(ok);

    % 夾在 [0, n-2]：落在最外面那一面時 i0 會等於 n-1，夾回去讓 t=1（數學等價）
    i0 = min(max(floor(fx), 0), nx-2);   tx = min(max(fx - i0, 0), 1);
    j0 = min(max(floor(fy), 0), ny-2);   ty = min(max(fy - j0, 0), 1);
    k0 = min(max(floor(fz), 0), nz-2);   tz = min(max(fz - k0, 0), 1);

    % --- 取 8 個角點（0 起算索引 -> MATLAB 要 +1）---
    a = i0 + 1;   b = j0 + 1;   c = k0 + 1;
    sz = [nx ny nz];
    L  = @(da,db,dc) sub2ind(sz, a+da, b+db, c+dc);
    V  = lat.V;                                   % (nx*ny*nz) x (3*N_I)
    V000 = V(L(0,0,0),:);   V100 = V(L(1,0,0),:);
    V010 = V(L(0,1,0),:);   V110 = V(L(1,1,0),:);
    V001 = V(L(0,0,1),:);   V101 = V(L(1,0,1),:);
    V011 = V(L(0,1,1),:);   V111 = V(L(1,1,1),:);

    % --- 三輪 lerp：8 -> 4 -> 2 -> 1（t 為 Np x 1，隱式擴張到 Np x 3*N_I）---
    lp = @(A,Bb,t) A + t .* (Bb - A);
    B00 = lp(V000, V100, tx);   B01 = lp(V001, V101, tx);     % 沿 x，4 次
    B10 = lp(V010, V110, tx);   B11 = lp(V011, V111, tx);
    C0  = lp(B00,  B10,  ty);   C1  = lp(B01,  B11,  ty);     % 沿 y，2 次
    Vq  = lp(C0,   C1,   tz);                                  % 沿 z，1 次

    Bq = reshape(Vq, [], 3, lat.N_I);            % 欄序 (c 快, m 慢) 與 lat.V 一致
end

% ============================================================================
function lat = lattice_cached(cfg, variant, rawin)
% 讀 .fld -> 組成規則格查表（persistent 快取，同 session 只讀一次）。
%   ⚠ 用座標反推 (i,j,k) 再散射，**與檔案列順序無關** —— 不假設 z 最快。
%   第 3 引數 rawin：呼叫者已載入的 raw（跳過 extract_maxwell_data）。
%   快取 key 改用「格點數 + 三個角落座標」當指紋，不同資料集不會互相污染。
    persistent LAT KEY
    if nargin < 3, rawin = []; end
    if isempty(rawin)
        lkey = [cfg.fld_dir '|' variant];         % key 須含 model
    else
        lkey = sprintf('raw|%d|%.12g|%.12g|%.12g|%.12g', numel(rawin.x), ...
                       rawin.x(1), rawin.y(1), rawin.z(1), rawin.z(end));
    end
    if ~isempty(KEY) && strcmp(KEY, lkey), lat = LAT;  return; end

    if isempty(rawin)
        raw = extract_maxwell_data(cfg, 'all', variant);      % [T]、Maxwell frame
    else
        raw = rawin;
    end
    N   = numel(raw.x);   N_I = size(raw.B,3);

    q  = @(vv) unique(round(vv*1e9))/1e9;                     % 量到 nm，去浮點雜訊
    xg = q(raw.x);   yg = q(raw.y);   zg = q(raw.z);
    n  = [numel(xg) numel(yg) numel(zg)];
    assert(prod(n) == N, 'conv_design_ws:notGrid', ...
           '.fld 不是完整規則格：%d x %d x %d = %d ~= %d 列', n(1),n(2),n(3),prod(n),N);

    o = [xg(1) yg(1) zg(1)];
    h = [(xg(end)-xg(1))/(n(1)-1), (yg(end)-yg(1))/(n(2)-1), (zg(end)-zg(1))/(n(3)-1)];

    ix = round((raw.x - o(1))/h(1)) + 1;                      % 1 起算
    iy = round((raw.y - o(2))/h(2)) + 1;
    iz = round((raw.z - o(3))/h(3)) + 1;
    % 驗證：反推的索引必須還原出原座標（不等距或有偏移就會抓到）
    dev = max([max(abs(o(1)+(ix-1)*h(1) - raw.x)), ...
               max(abs(o(2)+(iy-1)*h(2) - raw.y)), ...
               max(abs(o(3)+(iz-1)*h(3) - raw.z))]);
    assert(dev < 1e-9, 'conv_design_ws:notUniform', ...
           '格點與等距假設不符（最大偏差 %.3e m）', dev);

    lin = sub2ind(n, ix, iy, iz);
    V   = zeros(prod(n), 3*N_I);
    for m = 1:N_I
        for cc = 1:3
            V(lin, (m-1)*3 + cc) = raw.B(:,cc,m);             % 欄序：c 快、m 慢
        end
    end

    lat = struct('n',n, 'o',o, 'h',h, 'N_I',N_I, 'V',V);
    LAT = lat;   KEY = lkey;
    fprintf('  [lattice] %d x %d x %d、格距 %.1f/%.1f/%.1f um、z 範圍 %.4f ~ %.4f mm\n', ...
            n(1), n(2), n(3), h*1e6, zg(1)*1e3, zg(end)*1e3);
end

% ============================================================================
function solver_path()
    % [MODIFIED 2026-08-27] Tree moved to matlab\Flux\; both hardcoded paths no longer
    %   exist, so addpath silently failed and rmpath could not strip the real APDL copy
    %   -> risk of the APDL model_config/solve_* shadowing the Maxwell ones.
    CAL  = fileparts(fileparts(mfilename('fullpath')));   % ...\Flux\Maxwell
    APDL = fullfile(fileparts(CAL), 'APDL', 'Calibration_using_FEM_modeling');
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(APDL,'function'));  rmpath(fullfile(APDL,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
end

% ============================================================================
function report(info)
    if ~isempty(info.Nr)
        fprintf('  [axsh] R=%.0f um：Nr=%d -> 六軸 x %d 層%s，產 %d 點；r = %.1f ... %.1f um\n', ...
                info.R*1e6, info.Nr, info.Nr, ternary_(info.center,' + 中心',''), ...
                info.N, info.radii(1)*1e6, info.radii(end)*1e6);
    end
    fprintf('  [sample] %d 點 -> 保留 %d（鐵件 %d、盒外 %d）；frame=%s\n', ...
            info.N, info.npts_kept, info.n_iron, info.n_outbox, info.frame);
end

% ============================================================================
function v = getdef_(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% ============================================================================
function v = ternary_(c, a, b)
    if c, v = a; else, v = b; end
end
