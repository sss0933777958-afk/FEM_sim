function [x, y, z, B, tri, info] = conv_design_sensor(N_r, N_theta, N_z, opt)
% conv_design_sensor -- 感測器（圓柱）：決定內插點位置 + 三線性內插取場
% =========================================================================
%   [SPLIT 2026-08-23 使用者拍板] 由 conv_design.m 拆出（另一支 = conv_design_ws），
%     並把 Sensor_grid_sample.m **整支併進來**（該檔已刪除）。
%   [MERGE 2026-08-23 使用者拍板] 再併入 sphere_grid_sample.m 的**三線性內插**部分，
%     讓兩支 conv 完全對稱 —— 使用者原話：「sensor 多了 V 那段，其他跟 ws 一樣」。
%     **兩支各自自足，不共用外部引擎檔**（三線性那段兩邊各一份，是這個決定的代價）。
%
%   ⚠ **本函式只做兩件事：①決定內插點位置 ②三線性內插取場。**
%     組 V（投影到 n̂ 再取體積平均）交給 build_V_matrix —— 它逐 sensor 呼叫本函式
%     拿「點與場」，自己只做 V(i,kc) = S_hall * mean(b·n̂)。
%     收斂判斷與迴圈在 main.m：
%
%       [~,~,~,~,~,si] = conv_design_sensor(1,3,2, struct('ladder',150, ...
%                            'sensor_r',sr, 'axial_tol',at));       % 先取階梯表
%       Vref = build_V_matrix(..., [10 100 10], ...);               % 密參考
%       for q = 1:size(si.ladder,1)
%           V = build_V_matrix(..., si.ladder(q,:), ...);           % 內部呼叫本函式
%           if max(abs(V(:)-Vref(:))./abs(Vref(:))) < TOL, break; end   % 判斷在 main
%       end
%
%   ⚠ 這是**引擎**（無對應圖），勿當孤兒刪除。
%
%   ── Step 1：設計階梯（決定「幾個點」）────────────────────────
%     小格子的三個邊（半徑 r 處）：徑向 R^2/(2 r N_r)、方位 2*pi*r/N_theta、軸向 H/N_z
%     取**中位面積半徑** r = R/sqrt(2)（一半面積在內、一半在外）令三者相等：
%       徑向 = 方位：R^2/(2 r N_r) = 2*pi*r/N_theta  =>  N_theta = 2*pi*N_r
%       方位 = 軸向：2*pi*r/N_theta = H/N_z          =>  N_z = sqrt(2)*(H/R)*N_r
%     配比    w = N_r : N_theta : N_z = 1 : 2*pi : sqrt(2)*(H/R)
%     總點數  N = 2*pi*sqrt(2)*(H/R) * N_r^3  =>  N_r = round((N_target/prod(w))^(1/3))
%     走階梯一步：對 t./w 最小的那一軸 +1。種子慣用 [1 3 2]。
%
%   ── ⚠ 兩個退化地板（實測 2026-08-23，內建強制、沿階梯也不會掉回去）─────
%     N_z >= 2      N_z=1 時場平均誤差卡在 0.6% 的地板，**加再多徑向/方位點都沒用**
%                   （(4,48,1) 灌到 192 點仍是 0.399%）；一跨到 2 立刻掉 2.8 倍。
%                   成因：軸向中點法對曲率的殘差。
%     N_theta >= 3  等間距中點法在週期變數上是**譜精度**：N_theta 個點精確積掉所有
%                   m 不是 N_theta 倍數的 Fourier 模態，誤差來自 m = N_theta。
%                   實測（N_r=6, N_z=10 固定）：1 -> 6.36%（連 m=1 線性梯度都積錯）、
%                   2 -> 1.94%、3 -> 0.44%、6 -> 0.15%、33 -> 0.004%。
%     ⚠ [BEHAVIOUR CHANGE] 舊 Sensor_grid_sample 的 .N_target 反解只有 N_z 有下限 2，
%       N_theta 下限是 1；本檔把地板統一成 [1 3 2]（含 .N_target 模式）。
%
%   ── Step 2：取樣位置（三維都取格心）──────────────────────────
%     把 Jacobian 積起來、等分它的反導數，再反解：dV = r dr dtheta dz
%       u_k = (k-0.5)/N_r      -> r_k     = R*sqrt(u_k)     等分 r^2
%       v_i = (i-0.5)/N_theta  -> theta_i = 2*pi*v_i        等分 theta
%       w_j = (j-0.5)/N_z      -> z_j     = H*w_j           等分 z
%     取格心（那個 -0.5）避開退化點（r=0）且無偏（把該環/該層對半分）。
%     ⚠ 感測器圓柱是**單邊**的：底面在 sensor_pos、沿 +n 延伸 H，幾何中心在局部
%       z = H/2。預設 span='base' 即此慣例，與 build_V_matrix 的 a = axial_tol*rand 一致。
%
%   ── Step 3+4：定位 + 三線性內插（與 conv_design_ws 同一套實作）────────
%     規則格定位 -> 8 角點 -> 7 次 lerp -> B [T] -> x1e3 -> mT。
%     ⚠ 感測器這條**不濾鐵、不做 R_act 旋轉**：sensor 貼在磁極表面外，法線 n̂ 也是
%       全域框的量 -> 場留在全域框，與撒點 1:1 對應（少一個點都會讓體積平均失真）。
%
%   ── 用法 ────────────────────────────────────────────────────
%     [~,~,~,~,~,i] = conv_design_sensor(1,3,2, struct('ladder',150, ...
%                          'sensor_r',R,'axial_tol',H))          % 只出階梯表
%     [x,y,z]       = conv_design_sensor(3,20,3, o)              % 只產點（不取場）
%     [x,y,z,B]     = conv_design_sensor(3,20,3, o)              % 產點 + 取場
%
%   輸入
%     N_r,N_theta,N_z : 設計三元組。opt.step 有給時改當**種子**用；
%                       三個都給 [] 且有 opt.N_target 時由配比反解。
%     opt             : 選項 struct
%       .sensor_r   圓柱半徑 R [m]（預設 0.15e-3）
%       .axial_tol  圓柱高   H [m]（預設 0.10e-3）
%       .step     q   把三元組當種子，實際用階梯**第 q 級**的設計
%       .ladder   M   **只**回 info.ladder（Mx3 階梯表），不產點（x/y/z 空）
%       .N_target N   由配比反解三元組（三元組給 [] 時用）
%       .center   3x1 圓柱**底面**中心（全域座標 [m]）；預設 [0;0;0]
%       .axis     3x1 圓柱軸方向（內部單位化）；預設 [0;0;1]
%       .span     'base'（預設，局部 z in [0,H]）| 'center'（z in [-H/2,H/2]）
%       .model .geom .variant  資料路由（只在要取場時用）
%       .raw      直接給已載入的 raw（強烈建議：sensor 場有 13.9M 格點）
%
%   輸出
%     x,y,z   N x 1 座標 [m]（全域座標；ladder 模式下為空）
%     B       N x 3 x N_I 磁場 [mT]（全域框，與 x/y/z 1:1 對應）
%             ⚠ 只有 nargout >= 4 才會實際內插；只要點的話別要第 4 個輸出。
%     tri     1x3   本次實際用的設計
%     info    struct .tri .npts .ladder .ratio .seed .step .R .H .aspect .V .h
%                    .r_k .theta_i .z_j .center .axis .span .N_target .n_outbox
% =========================================================================
    if nargin < 4 || isempty(opt), opt = struct(); end
    gv = @(f,d) getdef_(opt, f, d);

    R    = gv('sensor_r',  0.15e-3);
    H    = gv('axial_tol', 0.10e-3);
    validateattributes(R, {'numeric'}, {'scalar','positive','finite'});
    validateattributes(H, {'numeric'}, {'scalar','positive','finite'});

    W    = [1, 2*pi, sqrt(2)*H/R];             % N_r : N_theta : N_z（配比）
    FL   = [1, 3, 2];                          % ⚠ 退化地板（見檔頭）
    STEP = gv('step',     []);
    LADM = gv('ladder',   []);
    NTGT = gv('N_target', []);
    cen  = gv('center',   [0;0;0]);
    axv  = gv('axis',     [0;0;1]);
    SPAN = gv('span',     'base');
    assert(any(strcmpi(SPAN,{'base','center'})), ...
           'conv_design_sensor:span', 'opt.span 必為 ''base'' | ''center''');

    seed = [N_r N_theta N_z];
    if any(cellfun(@isempty, {N_r, N_theta, N_z})), seed = [1 3 2]; end
    validateattributes(seed, {'numeric'}, {'vector','numel',3,'positive','integer'}, ...
                       'conv_design_sensor', 'seed / 三元組');
    seed = max(seed(:).', FL);

    % ---- ladder 模式：只出階梯表、不產點 -----------------------------------
    if ~isempty(LADM)
        validateattributes(LADM, {'numeric'}, {'scalar','positive','integer'});
        x = [];  y = [];  z = [];  B = [];  tri = [];
        info = mkinfo([], [], R, H, W, seed, [], ladder_(W, seed, FL, LADM), ...
                      NTGT, cen, axv, SPAN, [], [], [], 0);
        return
    end

    % ---- Step 1：決定本次設計 ----------------------------------------------
    LAD = [];
    if ~isempty(STEP)
        validateattributes(STEP, {'numeric'}, {'scalar','positive','integer'});
        LAD = ladder_(W, seed, FL, STEP);
        tri = LAD(STEP,:);
    elseif any(cellfun(@isempty, {N_r, N_theta, N_z})) && ~isempty(NTGT)
        s   = (NTGT / prod(W))^(1/3);          % 三數同乘 s -> 乘積乘 s^3
        tri = max(round(W * s), FL);           % 地板統一 [1 3 2]（見檔頭）
    else
        tri = seed;
    end
    if isempty(NTGT), NTGT = prod(tri); end

    % ---- Step 2：產點（局部格心）→ 擺位（局部 -> 全域）---------------------
    [x, y, z, r_k, th_i, z_j] = mkpts(R, H, tri, SPAN, cen, axv);

    % ---- Step 3+4：定位 + 三線性內插（只在要場時才做）----------------------
    B = [];   nout = 0;
    if nargout >= 4
        solver_path();
        MODEL = gv('model', 'long2016_hexapole_halfcut');
        GEOM  = gv('geom',  '');
        if isempty(GEOM) && strcmp(MODEL,'long2016_hexapole_halfcut'), GEOM = 'tip40um'; end
        cfg     = model_config(MODEL, GEOM);
        VARIANT = gv('variant', cfg.default_variant);

        lat = lattice_cached(cfg, VARIANT, gv('raw', []));
        [Bt, inbox] = trilerp(lat, x, y, z);                     % [T]，全域座標
        nout = nnz(~inbox);
        assert(nout == 0, 'conv_design_sensor:outbox', ...
            ['%d / %d 個撒點落在 .fld 盒外 —— 場的匯出框沒蓋到感測器。' ...
             '感測器的體積平均需要每個點都有值，不可默默少算。'], nout, numel(x));
        B = 1e3 * Bt;                                            % T -> mT
    end

    info = mkinfo(tri, numel(x), R, H, W, seed, STEP, LAD, NTGT, cen, axv, SPAN, ...
                  r_k, th_i, z_j, nout);
end

% ============================================================================
function [x, y, z, r_k, theta_i, z_j] = mkpts(R, H, t, SPAN, cen, axv)
% Step 2：由三個整數產點（等測度格心）。
    N_r = t(1);   N_th = t(2);   N_z = t(3);
    u  = ((1:N_r)  - 0.5) / N_r;     r_k     = R * sqrt(u);     % 等分 r^2
    v  = ((1:N_th) - 0.5) / N_th;    theta_i = 2*pi * v;        % 等分 theta
    ww = ((1:N_z)  - 0.5) / N_z;     z_j     = H * ww;          % 等分 z
    if strcmpi(SPAN, 'center'), z_j = z_j - H/2; end

    [K, I, J] = ndgrid(1:N_r, 1:N_th, 1:N_z);        % k 最快（同 (theta,z) 的不同
    rr = r_k(K(:));                                  %   半徑相鄰）；順序不影響結果
    tt = theta_i(I(:));
    zl = z_j(J(:));
    xl = rr(:) .* cos(tt(:));                        % 圓柱局部座標
    yl = rr(:) .* sin(tt(:));
    zl = zl(:);

    %   ⚠ 正交基底與 build_V_matrix 的亂數版 **逐字相同**（含 flat 面 n=[0;0;1]
    %     的退化 fallback）—— 這樣網格點雲與隨機點雲落在同一個局部框，
    %     兩者可以逐點比較；改了就不可比。
    a  = axv(:) / norm(axv(:));
    t1 = [-a(2); a(1); 0];   if norm(t1) < 1e-9, t1 = [1;0;0]; end
    t1 = t1 / norm(t1);      t2 = cross(a, t1);

    P = cen(:).' + zl.*a.' + xl.*t1.' + yl.*t2.';    % N x 3
    x = P(:,1);   y = P(:,2);   z = P(:,3);
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
function [Bq, ok] = trilerp(lat, xq, yq, zq)
% Step 3+4：定位 + 三線性內插（7 次 lerp）。xq/yq/zq 為 Np x 1 **全域座標** [m]。
%   ⚠ 與 conv_design_ws 的同名 local function **逐字相同**（使用者拍板：兩支 conv
%     各自自足、不共用外部引擎檔）。改一邊要同步另一邊。
    nx = lat.n(1);   ny = lat.n(2);   nz = lat.n(3);

    fx = (xq - lat.o(1)) / lat.h(1);
    fy = (yq - lat.o(2)) / lat.h(2);
    fz = (zq - lat.o(3)) / lat.h(3);
    % ⚠ 容差不可省：落在**最外面那一面**的點算出來會是 60.0000000000001 這種值，
    %   無容差就會被誤判成盒外。TOL 以「格」為單位，1e-9 格 = 2e-14 m。
    TOL = 1e-9;
    ok = fx >= -TOL & fx <= nx-1+TOL & fy >= -TOL & fy <= ny-1+TOL ...
       & fz >= -TOL & fz <= nz-1+TOL;
    fx = fx(ok);   fy = fy(ok);   fz = fz(ok);

    i0 = min(max(floor(fx), 0), nx-2);   tx = min(max(fx - i0, 0), 1);
    j0 = min(max(floor(fy), 0), ny-2);   ty = min(max(fy - j0, 0), 1);
    k0 = min(max(floor(fz), 0), nz-2);   tz = min(max(fz - k0, 0), 1);

    a = i0 + 1;   b = j0 + 1;   c = k0 + 1;
    sz = [nx ny nz];
    L  = @(da,db,dc) sub2ind(sz, a+da, b+db, c+dc);
    V  = lat.V;                                   % (nx*ny*nz) x (3*N_I)
    V000 = V(L(0,0,0),:);   V100 = V(L(1,0,0),:);
    V010 = V(L(0,1,0),:);   V110 = V(L(1,1,0),:);
    V001 = V(L(0,0,1),:);   V101 = V(L(1,0,1),:);
    V011 = V(L(0,1,1),:);   V111 = V(L(1,1,1),:);

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
%   ⚠ 與 conv_design_ws 的同名 local function **逐字相同**，但**各有自己的 persistent
%     快取** —— 這反而是好事：ws 吃 WP 細格（226981 點）、sensor 吃 voltage 粗格
%     （13.9M 點），兩份分開存就不會互相踢掉（單一快取時 main.m 每輪都要重建）。
    persistent LAT KEY
    if nargin < 3, rawin = []; end
    if isempty(rawin)
        lkey = [cfg.fld_dir '|' variant];
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
    assert(prod(n) == N, 'conv_design_sensor:notGrid', ...
           '.fld 不是完整規則格：%d x %d x %d = %d ~= %d 列', n(1),n(2),n(3),prod(n),N);

    o = [xg(1) yg(1) zg(1)];
    h = [(xg(end)-xg(1))/(n(1)-1), (yg(end)-yg(1))/(n(2)-1), (zg(end)-zg(1))/(n(3)-1)];

    ix = round((raw.x - o(1))/h(1)) + 1;
    iy = round((raw.y - o(2))/h(2)) + 1;
    iz = round((raw.z - o(3))/h(3)) + 1;
    dev = max([max(abs(o(1)+(ix-1)*h(1) - raw.x)), ...
               max(abs(o(2)+(iy-1)*h(2) - raw.y)), ...
               max(abs(o(3)+(iz-1)*h(3) - raw.z))]);
    assert(dev < 1e-9, 'conv_design_sensor:notUniform', ...
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
    APDL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    CAL  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(APDL,'function'));  rmpath(fullfile(APDL,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));  addpath(fullfile(CAL,'common_path'));
end

% ============================================================================
function info = mkinfo(tri, N, R, H, w, seed, step, LAD, Ntgt, cen, axv, SPAN, rk, th, zj, nout)
    V = pi * R^2 * H;
    h = [];   if ~isempty(N) && N > 0, h = (V/N)^(1/3); end
    info = struct('tri',tri, 'npts',N, 'R',R, 'H',H, 'aspect',H/R, 'V',V, 'h',h, ...
                  'ratio',w, 'seed',seed, 'step',step, 'ladder',LAD, 'N_target',Ntgt, ...
                  'r_k',rk, 'theta_i',th, 'z_j',zj, 'n_outbox',nout, ...
                  'center',cen(:), 'axis',axv(:)/norm(axv(:)), 'span',lower(SPAN));
end

% ============================================================================
function v = getdef_(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
