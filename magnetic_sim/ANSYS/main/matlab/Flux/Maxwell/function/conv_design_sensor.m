function [x, y, z, B, tri, info] = conv_design_sensor(Nx, Ny, opt)
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
%     設計現在是**定值**（使用者 2026-08-28 定案 Nx=Ny=100），main.m 不再跑 sensor 收斂迴圈：
%
%       V = build_V_matrix(..., [Nx Ny], ...);      % 內部呼叫本函式
%
%   ⚠ 這是**引擎**（無對應圖），勿當孤兒刪除。
%
%   ── Step 1：設計 = 手填 (Nx, Ny)（配比與階梯都已廢除）──────────
%     [REMOVED 2026-08-28 使用者拍板]「配比不用了，廢掉」-> 原本的
%       w = N_r : N_theta : N_z 以及 .ladder / .step / .N_target 全部移除。
%     [MODIFIED 2026-08-28 使用者拍板] 三軸圓柱設計再改成**中心面二維設計**，
%     參數由三個變兩個：只給 (Nx, Ny)。傳 .ladder/.step/.N_target 會直接報錯。
%
%   ── ⚠ 退化情形 ─────────────────────────────────────────────
%     (Nx,Ny) = (1,1) -> 候選只有四個角 (±R,±R)，距離 R*sqrt(2) 全在圓外 -> **0 點**，
%     本函式會明報錯。最小可用設計是 (2,2) -> 5 點。
%     另外**奇數等分沒有軸上的點**（Nx=3 的格線是 -R,-R/3,+R/3,+R，沒有 0），
%     所以點數不是單調的：Nx=Ny= 2->5、3->4、4->13、6->29、10->81、100->7845。
%
%   ── Step 2：取樣位置（中心面笛卡兒格）────────────────────────
%     [MODIFIED 2026-08-28 使用者拍板] 在感測區域的**中心面**上鋪笛卡兒格：
%       x_i = -R + 2R*i/Nx     i = 0..Nx      （格線 Nx+1 條）
%       y_j = -R + 2R*j/Ny     j = 0..Ny
%       保留 x^2 + y^2 <= R^2；全部落在 z = H/2（span='center' 時 z = 0）
%     -> 點數 ~ (pi/4)*(Nx+1)*(Ny+1)。定案 Nx=Ny=100 -> 7845 點、點距 3 um。
%     ✅ 這個擺法修掉了前一版（三軸等分）的兩個加權偏差：中心面 = 軸向中點法
%        （前版 z=H*j/Nz 是右端點、整層偏頂面）、笛卡兒格 = 面積權重天生均勻
%        （前版等分半徑 + 未加權平均會偏內圈）。
%     ⚠ 但單一平面**不是體積平均**：與整個圓柱的真值差一個軸向曲率項
%        (~H^2/24 * d2f/dz2)。實測（long2016、soff 4.572mm、對 40 萬點均勻亂數）
%        這個地板是 **0.45%**，且**與面內密度無關**：Nx=Ny=8(49點) 0.442%、
%        20(317點) 0.455%、100(7845點) 0.4514%、400(125629點) 0.4517%。
%        面內本身的收斂則很快：離「面內真值(Nx=Ny=400)」在 Nx=Ny=8 就 0.043%、
%        20 是 0.015%、100 是 0.001%。要壓那 0.45% 得動軸向（多切幾片），不是加密面內。
%     ⚠ 感測器圓柱是**單邊**的：底面在 sensor_pos、沿 +n 延伸 H，幾何中心在局部
%       z = H/2。預設 span='base' 即此慣例，與 build_V_matrix 的 a = axial_tol*rand 一致。
%
%   ── Step 3+4：定位 + 三線性內插（與 conv_design_ws 同一套實作）────────
%     規則格定位 -> 8 角點 -> 7 次 lerp -> B [T] -> x1e3 -> mT。
%     ⚠ 感測器這條**不濾鐵、不做 R_act 旋轉**：sensor 貼在磁極表面外，法線 n̂ 也是
%       全域框的量 -> 場留在全域框，與撒點 1:1 對應（少一個點都會讓體積平均失真）。
%
%   ── 用法 ────────────────────────────────────────────────────
%     [x,y,z]       = conv_design_sensor(100,100, o)             % 只產點（7845 點）
%     [x,y,z,B]     = conv_design_sensor(100,100, o)             % 產點 + 取場
%     [x,y,z]       = conv_design_sensor(2,2, o)                 % 最小可用設計 -> 5 點
%     opt             : 選項 struct
%       .sensor_r   圓柱半徑 R [m]（預設 0.15e-3）
%       .axial_tol  圓柱高   H [m]（預設 0.10e-3）
%       .center   3x1 圓柱**底面**中心（全域座標 [m]）；預設 [0;0;0]
%       .axis     3x1 圓柱軸方向（內部單位化）；預設 [0;0;1]
%       .span     'base'（預設，局部 z in [0,H]）| 'center'（z in [-H/2,H/2]）
%       .model .geom .variant  資料路由（只在要取場時用）
%       .raw      直接給已載入的 raw（強烈建議：sensor 場有 13.9M 格點）
%
%   輸出
%     x,y,z   N x 1 座標 [m]（全域座標），N = 圓內格點數 ~ (pi/4)*(Nx+1)*(Ny+1)
%     B       N x 3 x N_I 磁場 [mT]（全域框，與 x/y/z 1:1 對應）
%             ⚠ 只有 nargout >= 4 才會實際內插；只要點的話別要第 4 個輸出。
%     tri     1x2   本次實際用的設計 [Nx Ny]
%     info    struct .tri .npts .R .H .aspect .A .h(等效面內點距)
%                    .x_i .y_j .z0 .center .axis .span .n_outbox
%             （.ratio/.seed/.step/.ladder/.N_target 與 .r_k/.theta_i/.z_j 都已移除）
% =========================================================================
    if nargin < 3 || isempty(opt), opt = struct(); end
    gv = @(f,d) getdef_(opt, f, d);

    R    = gv('sensor_r',  0.15e-3);
    H    = gv('axial_tol', 0.10e-3);
    validateattributes(R, {'numeric'}, {'scalar','positive','finite'});
    validateattributes(H, {'numeric'}, {'scalar','positive','finite'});

    cen  = gv('center',   [0;0;0]);
    axv  = gv('axis',     [0;0;1]);
    SPAN = gv('span',     'base');
    assert(any(strcmpi(SPAN,{'base','center'})), ...
           'conv_design_sensor:span', 'opt.span 必為 ''base'' | ''center''');

    % [REMOVED 2026-08-28 使用者拍板]「配比廢掉」-> .ladder / .step / .N_target 全部失效。
    for f = {'ladder','step','N_target'}
        if isfield(opt, f{1}) && ~isempty(opt.(f{1}))
            error('conv_design_sensor:ratioScrapped', ...
                ['opt.%s 已於 2026-08-28 隨「配比」一起廢除；' newline ...
                 '設計現在只接受 (Nx, Ny) 兩個整數（中心面笛卡兒格）。'], f{1});
        end
    end

    tri = [Nx Ny];
    assert(~any(cellfun(@isempty, {Nx, Ny})), 'conv_design_sensor:needNxNy', ...
           '(Nx, Ny) 必須明確給定。');
    validateattributes(tri, {'numeric'}, {'vector','numel',2,'positive','integer'}, ...
                       'conv_design_sensor', '(Nx, Ny)');
    tri = tri(:).';

    % ---- Step 2：產點（中心面笛卡兒格）→ 擺位（局部 -> 全域）--------------
    [x, y, z, x_i, y_j, z0] = mkpts(R, H, tri, SPAN, cen, axv);
    assert(~isempty(x), 'conv_design_sensor:emptyDesign', ...
           ['(Nx,Ny) = (%d,%d) 沒有任何格點落在圓內（Nx=Ny=1 只有四個角、全在圓外）。' ...
            newline '請至少給 (2,2)（-> 5 點）。'], tri(1), tri(2));

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

    info = mkinfo(tri, numel(x), R, H, cen, axv, SPAN, x_i, y_j, z0, nout);
end

% ============================================================================
function [x, y, z, x_i, y_j, z0] = mkpts(R, H, t, SPAN, cen, axv)
% Step 2：中心面笛卡兒格。[MODIFIED 2026-08-28 使用者拍板]
%   x 切 Nx 等分、y 切 Ny 等分（格線各 Nx+1 / Ny+1 條、範圍 [-R,R]），
%   只留 x^2 + y^2 <= R^2 的節點；全部落在中心面 z = H/2（span='center' 時 z = 0）。
    Nx = t(1);   Ny = t(2);
    x_i = linspace(-R, R, Nx+1);
    y_j = linspace(-R, R, Ny+1);
    z0  = H/2;   if strcmpi(SPAN, 'center'), z0 = 0; end

    [X, Y] = ndgrid(x_i, y_j);
    keep = hypot(X(:), Y(:)) <= R*(1 + 1e-12);       % 容差：邊界上的點算圓內
    xl = X(keep);   yl = Y(keep);   zl = z0 * ones(nnz(keep), 1);

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
% [REMOVED 2026-08-28 使用者拍板] local ladder_（每步對 t./w 最小的那軸 +1）隨「配比」
%   一起廢除。新的階梯規則待定；決定之前 conv_design_sensor 只接受手填三元組。

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
function info = mkinfo(tri, N, R, H, cen, axv, SPAN, x_i, y_j, z0, nout)
% [MODIFIED 2026-08-28] 配比廢除 -> .ratio/.seed/.step/.ladder/.N_target 移除；
%   取樣改中心面笛卡兒格 -> .r_k/.theta_i/.z_j 換成 .x_i/.y_j/.z0、.A 換算面積。
    A = pi * R^2;
    h = [];   if ~isempty(N) && N > 0, h = sqrt(A/N); end     % 等效點距（面內）
    info = struct('tri',tri, 'npts',N, 'R',R, 'H',H, 'aspect',H/R, 'A',A, 'h',h, ...
                  'x_i',x_i, 'y_j',y_j, 'z0',z0, 'n_outbox',nout, ...
                  'center',cen(:), 'axis',axv(:)/norm(axv(:)), 'span',lower(SPAN));
end

% ============================================================================
function v = getdef_(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
