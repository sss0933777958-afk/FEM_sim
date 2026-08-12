function [x, y, z, B, info] = sphere_grid_sample(R, c, opt)
% sphere_grid_sample -- 球內「等測度網格」查詢點 + 該位置的磁場（規則格三線性內插）
% =========================================================================
%   給取樣半徑 R 與過取樣倍率 c，自動配出 (N_r, N_phi, N_theta)、產生查詢點、
%   從 Maxwell 匯出格內插取場，回傳座標與磁場。
%
%     [x, y, z, B, info] = sphere_grid_sample(150e-6, 2)
%
%   ── 六步流程 ──────────────────────────────────────────────────
%     1. 建查詢點（含配比）   數格點 -> N_target -> 配 (N_r,N_phi,N_theta) -> 三個
%                             等分公式 -> 座標，z 再 + SPH_OFST 直接成**全域座標**
%     2. 濾鐵                 filter_iron_nodes(x,y,z,cfg)，吃全域座標、零轉換
%     3. 定位                 fx=(x-x0)/h -> i=floor(fx), tx=fx-i；取 8 個角點
%     4. 三線性內插           7 次 lerp（4 沿 x -> 2 沿 y -> 1 沿 z）-> B [T]
%     5. 單位換算             x1e3 -> mT
%     6. 最後轉換             平移 z - SPH_OFST，再旋轉 R_act -> actuator frame
%
%   ── Step 1-a：點數怎麼定（綁在資料自己的解析度上）───────────────
%     N_target = c x N_nodes(R)
%       N_nodes = 該球內**空氣中**的 .fld 原始格點數（Maxwell 匯出格距 20 um）
%       c = 1 即「查詢密度 = FEM 匯出格密度」；c = 2 即過取樣 100%
%     因 N_nodes 正比於體積，這等價於「固定密度」：h = (V/N)^(1/3) = 20um / c^(1/3)
%     ⚠ h 只隨 c^(1/3) 縮 —— 點數 x8 才讓點距減半，報酬遞減很快。
%     ⚠ 過取樣不會壞擬合（內插器是由格點值完全決定的函數，多取點只是把同一個函數
%       多求值幾次，等於改善體積積分的品質），但**不可**把 N 說成「獨立取樣點數」、
%       也不可用 sqrt(N) 的統計 —— 獨立資訊上限就是 N_nodes。
%
%   ── Step 1-b：三個整數怎麼配（格子盡量方正）──────────────────
%     小格子的三個邊（半徑 r、緯度 phi 處）：
%       徑向  Delta_r = R^3/(3 r^2 N_r)              （等體積殼 -> 正比 1/r^2）
%       南北  r*Delta_phi = 2r/(N_phi*sin(phi))      （等分 cos phi -> 帶 1/sin）
%       東西  r*sin(phi)*Delta_theta = 2*pi*r*sin(phi)/N_theta
%
%     「南北 = 東西」，取赤道 sin(phi)=1：
%        2*pi*r*sin/N_theta = 2r/(N_phi*sin) => N_theta = pi*N_phi*sin^2 => N_theta = pi*N_phi
%     「徑向 = 南北」，取體積中位半徑 (r/R)^3 = 1/2（一半體積在內、一半在外）：
%        R^3/(3 r^2 N_r) = 2r/N_phi          => N_phi = 6*N_r*(r/R)^3 => N_phi = 3*N_r
%
%     比例固定       N_r : N_phi : N_theta = 1 : 3 : 3*pi ~ 1 : 3 : 9.42
%     比例只給形狀，大小由總點數定：三數同乘 s -> 乘積乘 s^3，故
%       N = 9*pi*N_r^3 ~ 28.27*N_r^3    =>   N_r = round( (N_target/(9*pi))^(1/3) )
%
%     ⚠ 兩個基準（赤道、體積中位半徑）都是取捨：sin(phi) 與 1/r^2 讓格子不可能處處
%       方正。極區方位過密、內殼徑向細長是先天限制 —— 但**密度仍完全均勻**。
%     ⚠ N 正比於 N_r^3，階梯很粗（1824 -> 3525 -> 6156），實際 c 會偏離設定值。
%
%   ── Step 1-c：取樣位置（三維都取格心）────────────────────────
%     把 Jacobian 積起來、等分它的反導數，再反解：dV = r^2*sin(phi) dr dphi dtheta
%       維度   Jacobian   累積測度            等分誰     反解
%       r      r^2        int r^2 dr ~ r^3    r^3        r     = R*u^(1/3)
%       phi    sin(phi)   int sin dphi=-cos   cos(phi)   phi   = acos(1-2w)
%       theta  1          int dtheta = theta  theta      theta = 2*pi*v
%
%       u_k = (k-0.5)/N_r      -> r_k     = R * u_k^(1/3)          k = 1..N_r
%       w_j = (j-0.5)/N_phi    -> phi_j   = acos(1 - 2*w_j)        j = 1..N_phi
%       v_i = (i-0.5)/N_theta  -> theta_i = 2*pi*v_i               i = 1..N_theta
%
%     取格心（那個 -0.5）避開退化點（r=0、phi=0/pi）且無偏（把該層/該帶對半分）。
%     ⚠ 中點必須取在**測度**上再反解，不可取半徑或角度的算術中點
%       （最內層：u 中點 = 0.550R 對半分；r 算術中點 0.347R 只含該層 12.5% 體積）。
%
%   ── 均勻性 ─────────────────────────────────────────────────────
%     OK 體積密度：每殼等體積 x 每殼點數相同(N_phi*N_theta) -> 常數
%     OK 球面密度：每帶等面積 x 每帶點數相同(N_theta)       -> 常數
%     NG 點距各向同性：(a) 極區方位擠壓 正比 sin(phi_j)
%                      (b) 各殼共用同一組方向 -> 點串在 N_phi*N_theta 根輻條上
%     ⚠ 不同 R 之間「點數相同 != 密度相同」：密度 = N/((4/3)pi R^3)，R 加倍而 N 不變
%       則密度掉 8 倍。本函式固定 c（= 固定密度），N 會隨 R^3 長。
%
%   ── 用法 ───────────────────────────────────────────────────────
%     [x,y,z,B]      = sphere_grid_sample(150e-6)               % c 預設 2
%     [x,y,z,B]      = sphere_grid_sample(150e-6, 1)            % 不過取樣
%     [x,y,z,B,info] = sphere_grid_sample(150e-6, 2, struct('frame','measure'))
%     [x,y,z,B]      = sphere_grid_sample([],[], struct('query',Q))  % 驗證用：自訂查詢點
%
%   輸入
%     R    取樣球半徑 [m]
%     c    過取樣倍率（預設 2）
%     opt  選項 struct（皆可省略）
%       .frame      'actuator'（預設）| 'measure'   輸出座標與 B 的 frame
%       .drop_iron  true/false  是否剔除落在鐵件內的點（預設 true）
%       .variant    .fld 變體（預設 cfg.default_variant = 'maxwell'）
%       .N_target   直接指定目標點數（跳過數格點；此時 c 被忽略）
%       .N_nodes    直接給該 R 的原始格點數（跳過讀 .fld 數點，仍用 c 算 target）
%       .query      Np x 3 **measure frame** 座標，直接指定查詢點（跳過 Step 1 的
%                   配比與產點）。驗證用（例：拿 .fld 格點當查詢點，應原值重現）。
%
%   輸出
%     x,y,z  Np x 1 座標 [m]（frame 依 opt.frame；'actuator' 時球心為原點）
%     B      Np x 3 x 6 磁場 [mT]，第 3 維 = 6 個激發（paper 序 P1..P6）
%     info   struct：.R .c_set .c_actual .N_target .N_nodes .N .N_r .N_phi .N_theta
%                    .r_k .phi_j .theta_i .h .n_kept .n_iron .n_outbox .keep
%
%   ⚠ 回傳點數可能少於 N：落在鐵件內、或落在 .fld 格盒外的點會被剔除（用 info.keep
%     對回原本的產點順序）。R > 401 um 時磁極會伸進取樣球。
%   ⚠ 磁場是**內插**值（匯出格距 20 um），非 FEM 節點原值。內插不增加資訊。
% =========================================================================
    if nargin < 2 || isempty(c),   c   = 2;         end
    if nargin < 3 || isempty(opt), opt = struct();  end
    gv = @(f,d) getdef(opt, f, d);
    FRAME     = gv('frame',     'actuator');
    DROP_IRON = gv('drop_iron', true);
    assert(any(strcmpi(FRAME,{'actuator','measure'})), 'opt.frame 必為 actuator | measure');

    solver_path();
    cfg     = model_config('long2016_hexapole_halfcut','tip40um');
    VARIANT = gv('variant', cfg.default_variant);
    QRY     = gv('query', []);

    info = struct('R',R, 'c_set',c, 'variant',VARIANT, 'frame',lower(FRAME), ...
                  'N_nodes',[], 'N_target',[], 'N_r',[], 'N_phi',[], 'N_theta',[], ...
                  'r_k',[], 'phi_j',[], 'theta_i',[], 'h',[], 'c_actual',[], ...
                  'n_iron',0, 'n_outbox',0);

    %% ==== Step 1：配比 + 建查詢點（直接建在全域座標）====================
    if isempty(QRY)
        validateattributes(R, {'numeric'}, {'scalar','positive'});

        % --- 1-a 目標點數：c x 該球內空氣中的原始格點數 ---
        N_target = gv('N_target', []);
        N_nodes  = gv('N_nodes',  []);
        if isempty(N_target)
            if isempty(N_nodes), N_nodes = count_fld_nodes(cfg, VARIANT, R); end
            N_target = c * N_nodes;
        end
        assert(N_target >= 27, 'N_target 太小（%g），至少 27 才配得出 (1,3,9)', N_target);

        % --- 1-b 配比 N_r : N_phi : N_theta = 1 : 3 : 3pi ---
        %     比例三項乘積 = 1*3*3pi = 9pi；三數同乘 s 則乘積乘 s^3
        s       = (N_target / (9*pi))^(1/3);
        N_r     = max(1, round(s));
        N_phi   = max(1, round(3    * N_r));
        N_theta = max(1, round(3*pi * N_r));

        % --- 1-c 三個 1D 座標（都取格心）+ 張成點雲 ---
        u = ((1:N_r)     - 0.5) / N_r;      r_k     = R * u.^(1/3);     % 等分 r^3
        w = ((1:N_phi)   - 0.5) / N_phi;    phi_j   = acos(1 - 2*w);    % 等分 cos(phi)
        v = ((1:N_theta) - 0.5) / N_theta;  theta_i = 2*pi * v;         % 等分 theta

        [K, J, I] = ndgrid(1:N_r, 1:N_phi, 1:N_theta);   % k 最慢、i 最快 -> 同殼相鄰
        rr = r_k(K(:));   pp = phi_j(J(:));   tt = theta_i(I(:));
        x = rr(:) .* sin(pp(:)) .* cos(tt(:));           % **全域座標**
        y = rr(:) .* sin(pp(:)) .* sin(tt(:));
        z = rr(:) .* cos(pp(:)) + cfg.SPH_OFST;          % 球心放 z = SPH_OFST

        N  = numel(x);   V = (4/3)*pi*R^3;
        info.N_nodes = N_nodes;   info.N_target = N_target;
        info.N_r = N_r;   info.N_phi = N_phi;   info.N_theta = N_theta;
        info.r_k = r_k;   info.phi_j = phi_j;   info.theta_i = theta_i;
        info.h   = (V/N)^(1/3);
        if ~isempty(N_nodes), info.c_actual = N / N_nodes; end
    else
        % 驗證/通用介面：直接給 measure frame 的查詢點
        assert(size(QRY,2)==3, 'opt.query 必須是 Np x 3（measure frame [m]）');
        x = QRY(:,1);   y = QRY(:,2);   z = QRY(:,3) + cfg.SPH_OFST;
    end
    info.N = numel(x);

    %% ==== Step 2：濾鐵（filter_iron_nodes 吃全域座標）====================
    keep = true(numel(x),1);
    if DROP_IRON
        keep = filter_iron_nodes(x, y, z, cfg);   keep = keep(:);
        info.n_iron = nnz(~keep);
    end

    %% ==== Step 3+4：定位 + 三線性內插 ====================================
    lat = lattice_cached(cfg, VARIANT);           % 61^3 規則格（persistent 快取）
    [Bt, inbox] = trilerp(lat, x(keep), y(keep), z(keep));       % [T]
    info.n_outbox = nnz(~inbox);
    idx = find(keep);   keep(idx(~inbox)) = false;               % 盒外一併剔除

    %% ==== Step 5：單位換算 T -> mT ========================================
    B = 1e3 * Bt;                                 % Bt 已只含 in-box 的列

    %% ==== Step 6：平移 z - SPH_OFST，再旋轉 R_act ========================
    x = x(keep);   y = y(keep);   z = z(keep) - cfg.SPH_OFST;    % (1) 平移
    if strcmpi(FRAME,'actuator')                                 % (2) 旋轉
        P = (cfg.R_act * [x, y, z].').';
        x = P(:,1);   y = P(:,2);   z = P(:,3);
        for m = 1:size(B,3)
            B(:,:,m) = (cfg.R_act * B(:,:,m).').';               % B 只旋轉
        end
    end

    info.keep = keep;   info.n_kept = numel(x);
    report(info);
end

% ============================================================================
function N_nodes = count_fld_nodes(cfg, variant, R)
% Step 1-a：數該球內**空氣中**的 .fld 原始格點（座標 persistent 快取）。
    persistent XYZ KEY
    if isempty(KEY) || ~strcmp(KEY, variant)
        raw = extract_maxwell_data(cfg, 'all', variant);
        XYZ = [raw.x, raw.y, raw.z];              % 全域座標 [m]
        KEY = variant;
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
% Step 3+4：定位 + 三線性內插（7 次 lerp）。xq/yq/zq 為 Np x 1 **全域座標** [m]。
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

    % --- 三輪 lerp：8 -> 4 -> 2 -> 1（t 為 Np x 1，隱式擴張到 Np x 18）---
    lp = @(A,Bb,t) A + t .* (Bb - A);
    B00 = lp(V000, V100, tx);   B01 = lp(V001, V101, tx);     % 沿 x，4 次
    B10 = lp(V010, V110, tx);   B11 = lp(V011, V111, tx);
    C0  = lp(B00,  B10,  ty);   C1  = lp(B01,  B11,  ty);     % 沿 y，2 次
    Vq  = lp(C0,   C1,   tz);                                  % 沿 z，1 次

    Bq = reshape(Vq, [], 3, lat.N_I);            % 欄序 (c 快, m 慢) 與 lat.V 一致
end

% ============================================================================
function lat = lattice_cached(cfg, variant)
% 讀 .fld -> 組成規則格查表（persistent 快取，同 session 只讀一次）。
%   ⚠ 用座標反推 (i,j,k) 再散射，**與檔案列順序無關** —— 不假設 z 最快。
    persistent LAT KEY
    if ~isempty(KEY) && strcmp(KEY, variant), lat = LAT;  return; end

    raw = extract_maxwell_data(cfg, 'all', variant);          % [T]、Maxwell frame
    N   = numel(raw.x);   N_I = size(raw.B,3);

    q  = @(vv) unique(round(vv*1e9))/1e9;                     % 量到 nm，去浮點雜訊
    xg = q(raw.x);   yg = q(raw.y);   zg = q(raw.z);
    n  = [numel(xg) numel(yg) numel(zg)];
    assert(prod(n) == N, 'sphere_grid_sample:notGrid', ...
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
    assert(dev < 1e-9, 'sphere_grid_sample:notUniform', ...
           '格點與等距假設不符（最大偏差 %.3e m）', dev);

    lin = sub2ind(n, ix, iy, iz);
    V   = zeros(prod(n), 3*N_I);
    for m = 1:N_I
        for cc = 1:3
            V(lin, (m-1)*3 + cc) = raw.B(:,cc,m);             % 欄序：c 快、m 慢
        end
    end

    lat = struct('n',n, 'o',o, 'h',h, 'N_I',N_I, 'V',V);
    LAT = lat;   KEY = variant;
    fprintf('  [lattice] %d x %d x %d、格距 %.1f/%.1f/%.1f um、z 範圍 %.4f ~ %.4f mm\n', ...
            n(1), n(2), n(3), h*1e6, zg(1)*1e3, zg(end)*1e3);
end

% ============================================================================
function solver_path()
    APDL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    CAL  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(APDL,'function'));  rmpath(fullfile(APDL,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
end

% ============================================================================
function v = getdef(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% ============================================================================
function report(info)
    if ~isempty(info.N_r)
        fprintf('  [grid] R=%.0f um：目標 %g 點 -> (N_r,N_phi,N_theta)=(%d,%d,%d)、產 %d 點\n', ...
                info.R*1e6, info.N_target, info.N_r, info.N_phi, info.N_theta, info.N);
        fprintf('         c 設定 %.2f -> 實際 %.3f（偏差 %+.2f%%）；h = %.2f um（原始格距 20 um）\n', ...
                info.c_set, info.c_actual, (info.N/info.N_target-1)*100, info.h*1e6);
        fprintf('         r_k [um] = %s\n', num2str(info.r_k*1e6, '%.1f '));
    end
    fprintf('  [sample] %d 點 -> 保留 %d（鐵件 %d、盒外 %d）；frame=%s\n', ...
            info.N, info.n_kept, info.n_iron, info.n_outbox, info.frame);
end
