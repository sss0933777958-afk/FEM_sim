function [P, Bstack, tri, info] = conv_design_ws(N_r, N_phi, N_theta, R, opt)
% conv_design_ws -- 工作空間（球）：決定內插點位置 + 三線性內插取場
% =========================================================================
%   [SPLIT  2026-08-23 使用者拍板] 由 conv_design.m 拆出（另一支 = conv_design_sensor）。
%     拆的理由：兩個項目的第 2、3 軸意義不同（球 = 極角/方位、圓柱 = 方位/軸向），
%     用同一組參數名必然誤導。**兩支各自自足，不共用外部引擎檔。**
%   [MERGE  2026-08-23 使用者拍板] **sphere_grid_sample.m 整支併進來**（該檔已刪除）：
%     設計階梯 + 產點 + 濾鐵 + 規則格三線性內插 + frame 轉換，全部在本檔。
%
%   ⚠ **本函式只做兩件事：①決定內插點位置 ②三線性內插取場。**
%     校正（fitting / solve_*）與收斂判斷**一律在 main.m**，迴圈也由 main.m 驅動：
%
%       [~,~,~,wi] = conv_design_ws(1,2,3, R, struct('ladder',150));   % 先取階梯表
%       for q = 1:size(wi.ladder,1)
%           t = wi.ladder(q,:);
%           [P,Bs] = conv_design_ws(t(1),t(2),t(3), R, o);   % 產點 + 取場
%           ... main 自己 fitting + solve_* ...              % 校正
%           ... main 自己累積序列、判收斂，過了就 break ...   % 判斷
%       end
%
%   ⚠ 這是**引擎**（無對應圖），被 main.m 與 5 支繪圖腳本呼叫，勿當孤兒刪除。
%
%   ── 六步流程 ────────────────────────────────────────────────
%     1. 決定三個整數     階梯 / 明給 / 由 N_target(或 c) 反解配比
%     2. 建查詢點         三個等分公式 -> 座標，z 再 + SPH_OFST 成**全域座標**
%     3. 濾鐵             filter_iron_nodes(x,y,z,cfg)，吃全域座標、零轉換
%     4. 定位             fx=(x-x0)/h -> i=floor(fx), tx=fx-i；取 8 個角點
%     5. 三線性內插       7 次 lerp（4 沿 x -> 2 沿 y -> 1 沿 z）-> B [T] -> x1e3 -> mT
%     6. 最後轉換         平移 z - SPH_OFST，再旋轉 R_act -> actuator frame
%
%   ── Step 1：設計階梯（決定「幾個點」）────────────────────────
%     小格子的三個邊（半徑 r、緯度 phi 處）：
%       徑向  Delta_r = R^3/(3 r^2 N_r)              （等體積殼 -> 正比 1/r^2）
%       南北  r*Delta_phi = 2r/(N_phi*sin(phi))      （等分 cos phi -> 帶 1/sin）
%       東西  r*sin(phi)*Delta_theta = 2*pi*r*sin(phi)/N_theta
%     「南北 = 東西」取赤道 sin(phi)=1  =>  N_theta = pi*N_phi
%     「徑向 = 南北」取體積中位半徑 (r/R)^3 = 1/2  =>  N_phi = 3*N_r
%     配比 w = N_r : N_phi : N_theta = 1 : 3 : 3*pi ~ 1 : 3 : 9.42
%     由 N_target 反解：N = 9*pi*N_r^3 => N_r = round((N_target/(9*pi))^(1/3))
%     走階梯一步：對 t./w 最小的那一軸 +1。種子慣用 [1 2 3]。
%     ⚠ 兩個基準（赤道、體積中位半徑）都是取捨：sin(phi) 與 1/r^2 讓格子不可能處處
%       方正。極區方位過密、內殼徑向細長是先天限制 —— 但**密度仍完全均勻**。
%     ⚠ N 正比於 N_r^3，階梯很粗（1824 -> 3525 -> 6156），實際 c 會偏離設定值。
%     ⚠ 球版**無退化地板**（不像圓柱版有 N_theta>=3 / N_z>=2 的硬下限）。
%
%   ── Step 2：取樣位置（三維都取格心）──────────────────────────
%     把 Jacobian 積起來、等分它的反導數，再反解：dV = r^2*sin(phi) dr dphi dtheta
%       u_k = (k-0.5)/N_r      -> r_k     = R * u_k^(1/3)          等分 r^3
%       w_j = (j-0.5)/N_phi    -> phi_j   = acos(1 - 2*w_j)        等分 cos(phi)
%       v_i = (i-0.5)/N_theta  -> theta_i = 2*pi*v_i               等分 theta
%     取格心（那個 -0.5）避開退化點（r=0、phi=0/pi）且無偏（把該層/該帶對半分）。
%     ⚠ 中點必須取在**測度**上再反解，不可取半徑或角度的算術中點
%       （最內層：u 中點 = 0.550R 對半分；r 算術中點 0.347R 只含該層 12.5% 體積）。
%
%   ── 均勻性 ──────────────────────────────────────────────────
%     OK 體積密度：每殼等體積 x 每殼點數相同(N_phi*N_theta) -> 常數
%     OK 球面密度：每帶等面積 x 每帶點數相同(N_theta)       -> 常數
%     NG 點距各向同性：(a) 極區方位擠壓 正比 sin(phi_j)
%                      (b) 各殼共用同一組方向 -> 點串在 N_phi*N_theta 根輻條上
%     ⚠ 不同 R 之間「點數相同 != 密度相同」：密度 = N/((4/3)pi R^3)。
%
%   ── 用法 ────────────────────────────────────────────────────
%     [~,~,~,i] = conv_design_ws(1,2,3, R, struct('ladder',150))   % 只出階梯表
%     [P,Bs]    = conv_design_ws(3,8,22, 150e-6)                   % 明給三元組
%     [P,Bs]    = conv_design_ws(1,2,3, R, struct('step',21))      % 階梯第 21 級
%     [P,Bs]    = conv_design_ws([],[],[], R, struct('c',2))       % 由過取樣倍率反解
%     [P,Bs]    = conv_design_ws([],[],[], [], struct('query',Q))  % 任意點取場
%
%   輸入
%     N_r,N_phi,N_theta : 設計三元組。opt.step 有給時改當**種子**用；
%                         三個都給 [] 且有 opt.N_target/opt.c 時由配比反解。
%     R                 : 取樣球半徑 [m]（query / ladder 模式外皆必給）
%     opt               : 選項 struct（皆可省略）
%       .ladder   M   **只**回 info.ladder（Mx3 階梯表），不產點（P/Bstack 空）
%       .step     q   把三元組當種子，實際用階梯**第 q 級**的設計
%       .N_target N   目標點數（三元組給 [] 時由配比反解；此時 c 被忽略）
%       .c            過取樣倍率（預設 2）：N_target = c x 該球內空氣中的原始格點數
%       .N_nodes      直接給該 R 的原始格點數（跳過讀 .fld 數點，仍用 c 算 target）
%       .query    Np x 3 **measure frame** 座標，直接指定查詢點（跳過 Step 1-2）
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
%     tri     1x3      本次實際用的設計（query 模式為 []）
%     info    struct   .tri .npts_design .npts_kept .ladder .ratio .seed .step .R
%                      .N_r .N_phi .N_theta .r_k .phi_j .theta_i .h .N_target
%                      .N_nodes .c_set .c_actual .n_iron .n_outbox .keep .frame .B
%                      （.B = Np x 3 x N_I 的原始形狀，需要時可直接取用）
%
%   ⚠ 回傳點數可能少於設計點數：落在鐵件內、或落在 .fld 格盒外的點會被剔除
%     （用 info.keep 對回原本的產點順序）。R > 401 um 時磁極會伸進取樣球。
%   ⚠ 磁場是**內插**值（匯出格距 20 um），非 FEM 節點原值。內插不增加資訊。
% =========================================================================
    if nargin < 4, R   = []; end
    if nargin < 5 || isempty(opt), opt = struct(); end
    gv = @(f,d) getdef_(opt, f, d);

    W    = [1, 3, 3*pi];                       % N_r : N_phi : N_theta（配比）
    FL   = [1 1 1];                            % 球版無退化地板
    STEP = gv('step',     []);
    LADM = gv('ladder',   []);
    QRY  = gv('query',    []);
    NTGT = gv('N_target', []);
    CVAL = gv('c',        2);
    NNOD = gv('N_nodes',  []);
    FRAME     = lower(gv('frame', 'actuator'));
    DROP_IRON = gv('drop_iron', true);
    QUIET     = logical(gv('quiet', true));
    assert(any(strcmp(FRAME,{'actuator','measure'})), 'opt.frame 必為 actuator | measure');

    seed = [N_r N_phi N_theta];
    if any(cellfun(@isempty, {N_r, N_phi, N_theta})), seed = [1 2 3]; end
    validateattributes(seed, {'numeric'}, {'vector','numel',3,'positive','integer'}, ...
                       'conv_design_ws', 'seed / 三元組');
    seed = max(seed(:).', FL);

    % ---- ladder 模式：只出階梯表、不產點 -----------------------------------
    if ~isempty(LADM)
        validateattributes(LADM, {'numeric'}, {'scalar','positive','integer'});
        P = [];   Bstack = [];   tri = [];
        info = struct('tri',[], 'npts_design',[], 'npts_kept',[], ...
                      'ladder',ladder_(W, seed, FL, LADM), 'ratio',W, 'seed',seed, ...
                      'step',[], 'R',R, 'frame',FRAME);
        return
    end

    solver_path();
    MODEL = gv('model', 'long2016_hexapole_halfcut');
    GEOM  = gv('geom',  '');
    if isempty(GEOM) && strcmp(MODEL,'long2016_hexapole_halfcut'), GEOM = 'tip40um'; end
    cfg     = model_config(MODEL, GEOM);
    VARIANT = gv('variant', cfg.default_variant);
    RAWIN   = gv('raw', []);

    info = struct('tri',[], 'npts_design',[], 'npts_kept',[], 'ladder',[], 'ratio',W, ...
                  'seed',seed, 'step',STEP, 'R',R, 'frame',FRAME, 'variant',VARIANT, ...
                  'N_r',[], 'N_phi',[], 'N_theta',[], 'r_k',[], 'phi_j',[], 'theta_i',[], ...
                  'h',[], 'N_target',[], 'N_nodes',[], 'c_set',CVAL, 'c_actual',[], ...
                  'n_iron',0, 'n_outbox',0, 'keep',[], 'N',[], 'B',[]);

    %% ==== Step 1+2：決定三個整數 -> 建查詢點（直接建在全域座標）==========
    LAD = [];   tri = [];
    if isempty(QRY)
        assert(~isempty(R) && isscalar(R) && R > 0, ...
               'conv_design_ws:noR', 'conv_design_ws 必須給取樣球半徑 R [m]');

        if ~isempty(STEP)                                  % 階梯第 q 級
            validateattributes(STEP, {'numeric'}, {'scalar','positive','integer'});
            LAD = ladder_(W, seed, FL, STEP);
            tri = LAD(STEP,:);
        elseif ~any(cellfun(@isempty, {N_r, N_phi, N_theta}))
            tri = seed;                                    % 明給三元組
        else                                               % 由 N_target / c 反解配比
            if isempty(NTGT)
                if isempty(NNOD), NNOD = count_fld_nodes(cfg, VARIANT, R); end
                NTGT = CVAL * NNOD;
            end
            assert(NTGT >= 27, 'N_target 太小（%g），至少 27 才配得出 (1,3,9)', NTGT);
            s   = (NTGT / prod(W))^(1/3);
            tri = max([round(s), round(W(2)*round(s)), round(W(3)*round(s))], FL);
        end
        if isempty(NTGT), NTGT = prod(tri); end

        % --- 三個 1D 座標（都取格心）+ 張成點雲 ---
        nr = tri(1);   np_ = tri(2);   nt = tri(3);
        u = ((1:nr)  - 0.5) / nr;     r_k     = R * u.^(1/3);     % 等分 r^3
        w = ((1:np_) - 0.5) / np_;    phi_j   = acos(1 - 2*w);    % 等分 cos(phi)
        v = ((1:nt)  - 0.5) / nt;     theta_i = 2*pi * v;         % 等分 theta

        [K, J, I] = ndgrid(1:nr, 1:np_, 1:nt);           % k 最慢、i 最快 -> 同殼相鄰
        rr = r_k(K(:));   pp = phi_j(J(:));   tt = theta_i(I(:));
        x = rr(:) .* sin(pp(:)) .* cos(tt(:));           % **全域座標**
        y = rr(:) .* sin(pp(:)) .* sin(tt(:));
        z = rr(:) .* cos(pp(:)) + cfg.SPH_OFST;          % 球心放 z = SPH_OFST

        Vol = (4/3)*pi*R^3;
        info.tri = tri;   info.npts_design = prod(tri);   info.ladder = LAD;
        info.N_r = nr;    info.N_phi = np_;   info.N_theta = nt;
        info.r_k = r_k;   info.phi_j = phi_j;  info.theta_i = theta_i;
        info.N_target = NTGT;   info.N_nodes = NNOD;   info.h = (Vol/numel(x))^(1/3);
        if ~isempty(NNOD), info.c_actual = numel(x) / NNOD; end
    else
        % 驗證/通用介面：直接給 measure frame 的查詢點
        assert(size(QRY,2) == 3, 'opt.query 必須是 Np x 3（measure frame [m]）');
        x = QRY(:,1);   y = QRY(:,2);   z = QRY(:,3) + cfg.SPH_OFST;
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
function LAD = ladder_(w, seed, fl, M)
% 設計階梯：每步對 t./w 最小的那一軸 +1（並守住地板 fl）。
    LAD = zeros(M,3);   t = max(seed, fl);
    for q = 1:M
        LAD(q,:) = t;
        [~,j] = min(t ./ w);   t(j) = t(j) + 1;
        t = max(t, fl);
    end
end

% ============================================================================
function N_nodes = count_fld_nodes(cfg, variant, R)
% 數該球內**空氣中**的 .fld 原始格點（座標 persistent 快取）。
    persistent XYZ KEY
    % key 加上 fld_dir：不同 model 的 variant 都叫 'maxwell'，只用 variant 當 key
    % 會讓第二個 model 拿到第一個的格點（跨模型污染）。
    ckey = [cfg.fld_dir '|' variant];
    if isempty(KEY) || ~strcmp(KEY, ckey)
        raw = extract_maxwell_data(cfg, 'all', variant);
        XYZ = [raw.x, raw.y, raw.z];              % 全域座標 [m]
        KEY = ckey;
    end
    r2  = XYZ(:,1).^2 + XYZ(:,2).^2 + (XYZ(:,3) - cfg.SPH_OFST).^2;
    in  = r2 <= R^2;
    air = filter_iron_nodes(XYZ(in,1), XYZ(in,2), XYZ(in,3), cfg);
    N_nodes = nnz(air);
    fprintf('  [nodes] R<=%.0f um 內原始格點 %d（球內 %d、扣鐵件 %d）\n', ...
            R*1e6, N_nodes, nnz(in), nnz(in)-N_nodes);
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
    if ~isempty(info.N_r)
        fprintf('  [grid] R=%.0f um：目標 %g 點 -> (N_r,N_phi,N_theta)=(%d,%d,%d)、產 %d 點\n', ...
                info.R*1e6, info.N_target, info.N_r, info.N_phi, info.N_theta, info.N);
        if ~isempty(info.c_actual)
            fprintf('         c 設定 %.2f -> 實際 %.3f（偏差 %+.2f%%）；h = %.2f um（原始格距 20 um）\n', ...
                    info.c_set, info.c_actual, (info.N/info.N_target-1)*100, info.h*1e6);
        end
        fprintf('         r_k [um] = %s\n', num2str(info.r_k*1e6, '%.1f '));
    end
    fprintf('  [sample] %d 點 -> 保留 %d（鐵件 %d、盒外 %d）；frame=%s\n', ...
            info.N, info.npts_kept, info.n_iron, info.n_outbox, info.frame);
end

% ============================================================================
function v = getdef_(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
