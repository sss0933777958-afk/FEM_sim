function [V, exc_sign, sensor_pos, sensor_n, dbg] = build_V_matrix(cfg, variant, raw, S_hall, SOFF_upper, n_uniform, sensor_r, axial_tol, sensor_override, V_METHOD, SOFF_lower, FACE_lower)
%   [ADDED 2026-08-06] 第 5 個輸出 dbg（僅 V_METHOD='scattered'）：逐點診斷用
%     dbg.ang{i,kc}  sensor i、激發 kc 下每個撒點的 ∠(b, n̂⁺) [deg]
%     dbg.bn{i,kc}   同上的 b·n̂⁺ [mT]；dbg.bmag{i,kc} 同上的 ‖b‖ [mT]
%     dbg.samp{i}    該 sensor 的撒點座標（ANSYS 框 [m]）——與算 V 用的是**同一批點**
%   [ADDED] SOFF_lower（可選，預設 4.572e-3）：下極 sensor 沿錐面距極尖 [m]。
%   原本下極寫死 4.572e-3、只有上極可調；要「兩層同時移動」時必須一起給。
%   [ADDED 2026-08-05] FACE_lower（可選，預設 'cone'）：下極 sensor 貼哪個面
%     'cone' = 半切件**底錐面**（原行為，法線朝下出鋼）
%     'flat' = 半切**平切上表面**（法線 +z；半切面是通過極軸的水平面 → 沿水平極軸走 SOFF、離面 AIR）
%   上極永遠貼自己的完整錐面，不受此參數影響。
%BUILD_V_MATRIX  電壓提取：sensor 幾何 → 每 sensor 圓柱撒點 → 內插 raw 場 → 6×6 V[mV]。
%   [V, exc_sign] = BUILD_V_MATRIX(cfg, variant, raw, S_hall, SOFF_upper, n_uniform, sensor_r, axial_tol, sensor_override, V_METHOD)
%   sensor 幾何來源優先序：
%     ① sensor_override（可選 struct .pos/.n，escape hatch）
%     ② cfg.sensor_pos / cfg.sensor_n（config「提供」的幾何，如 tip400um）
%     ③ utils/pole_sensor_geometry（**sensor 幾何唯一來源**；CAD 實測錐體 + 真實氣隙）
%   V_METHOD：'csv-tet'（預設；讀 sensor_local CSV 建 tet 重心內插——需 CSV 與 solve mesh 同網格）
%             'scattered'（座標式 scatteredInterpolant 直接對 raw solve 場取樣——CSV≠solve mesh 時用，如 tip400）。
%             'grid'（**建議**；2026-08-23 加）：Maxwell 的 .fld 是規則格 —— 直接「定位格子
%               + 三線性內插」，與工作空間那條（conv_design_ws）用**完全同一套**內插。
%               比 'scattered' 快得多：不必在 13.9M 格點上做 6 次鄰域距離搜尋、也不必建
%               108 次 Delaunay。⚠ 兩者都是分片線性但**不逐位相同**（Delaunay 把立方格
%               切四面體、三線性在立方格內做張量積），換過去 V 會有同階的微小變動。
%   需 ansys_path + readmatrix + triangulation/pointLocation（csv-tet）/ scatteredInterpolant（scattered）在 path。
    if nargin < 5  || isempty(SOFF_upper), SOFF_upper = 4.572e-3; end
    if nargin < 6  || isempty(n_uniform),  n_uniform  = 1e4;      end
    if nargin < 7  || isempty(sensor_r),   sensor_r   = 0.15e-3;  end
    if nargin < 8  || isempty(axial_tol),  axial_tol  = 0.10e-3;  end
    if nargin < 9,  sensor_override = []; end
    if nargin < 10 || isempty(V_METHOD),   V_METHOD   = 'csv-tet'; end
    if nargin < 11 || isempty(SOFF_lower), SOFF_lower = 4.572e-3;  end   % [ADDED]
    if nargin < 12 || isempty(FACE_lower), FACE_lower = 'cone';    end   % [ADDED 2026-08-05]
    N_I = cfg.N_I;

    % ① sensor 幾何（WP 框）：override > config 提供 > 內建 baseline 錐面
    if ~isempty(sensor_override)
        sensor_pos = sensor_override.pos;  sensor_n = sensor_override.n;
    elseif isfield(cfg,'sensor_pos') && isfield(cfg,'sensor_n') && ~isempty(cfg.sensor_pos)
        sensor_pos = cfg.sensor_pos;  sensor_n = cfg.sensor_n;
    else
        % [MODIFIED 2026-08-08] 改呼叫 utils/pole_sensor_geometry（sensor 幾何唯一來源），
        %   原本的 local sensor_geometry 已移除。
        [sensor_pos, sensor_n] = pole_sensor_geometry(cfg, struct( ...
            'soff_upper', SOFF_upper, 'soff_lower', SOFF_lower, 'face_lower', FACE_lower));
    end

    % ② 每 sensor 圓柱內取 n_uniform 點（ANSYS 框）——兩內插法共用
    % [MODIFIED 2026-08-21 使用者拍板] n_uniform 現在兩用：
    %   **純量** -> 舊的 Monte Carlo 亂數撒點（rng(0) 可重現）；既有結果逐位不變。
    %   **1x3**  -> [N_r N_theta N_z] 等測度網格（conv_design_sensor），確定性、
    %              密度精確均勻。下游完全不必改：兩種撒法都是「體積均勻」，
    %              所以後面那個未加權 mean 依然是合法的體積平均。
    %   ⚠ 網格版的正交基底與下面亂數版**逐字相同**（conv_design_sensor 內同一段），
    %     兩者落在同一個局部框，可逐點對照。
    GRID = ~isscalar(n_uniform);
    if ~GRID, rng(0); end
    cen = zeros(3,6);  samp = cell(1,6);  sfld = cell(1,6);
    for i = 1:6
        ci = sensor_pos(:,i) + [0;0;cfg.SPH_OFST];    % WP 框 → ANSYS 框
        cen(:,i) = ci;
        ni = sensor_n(:,i);
        if GRID
            % [MODIFIED 2026-08-23 使用者拍板] Sensor_grid_sample 與 sphere_grid_sample
            %   都已併入 conv_design_sensor（兩檔已刪除）。分工：conv_design_sensor
            %   做「產點 + 三線性內插取場」，本檔**只組 V**（投影到 n̂ 再取體積平均）
            %   —— 即「build_V_matrix 從 conv_design_sensor 拿資料」。
            [sx, sy, sz, sB] = conv_design_sensor(n_uniform(1), n_uniform(2), n_uniform(3), ...
                struct('sensor_r', sensor_r, 'axial_tol', axial_tol, ...
                       'center', ci, 'axis', ni, 'span', 'base', 'raw', raw, ...
                       'variant', variant, 'model', cfg.model, 'geom', cfg.geom));
            samp{i} = [sx, sy, sz];   sfld{i} = sB;
        else
            t1 = [-ni(2); ni(1); 0]; if norm(t1) < 1e-9, t1 = [1;0;0]; end
            t1 = t1/norm(t1);  t2 = cross(ni, t1);
            a  = axial_tol * rand(n_uniform,1);
            rr = sensor_r  * sqrt(rand(n_uniform,1));     % √U → 面積均勻
            th = 2*pi      * rand(n_uniform,1);
            samp{i} = ci.' + a.*ni.' + (rr.*cos(th)).*t1.' + (rr.*sin(th)).*t2.';
        end
    end
    if GRID
        fprintf('  [V] sensor 取樣＝等測度網格 (%d,%d,%d) = %d 點/sensor\n', ...
                n_uniform, size(samp{1},1));
    else
        fprintf('  [V] sensor 取樣＝亂數 %d 點/sensor（rng(0)）\n', n_uniform);
    end

    % ③ 內插 raw 場到撒點 → V=S_hall·⟨B·n̂⟩[mV]（依 V_METHOD 分兩法）
    dbg = struct('ang',{cell(6,N_I)}, 'bn',{cell(6,N_I)}, 'bmag',{cell(6,N_I)}, 'samp',{samp});
    switch V_METHOD
        case 'csv-tet'
            V = vmat_csv_tet(cfg, variant, raw, samp, sensor_n, S_hall, N_I);
        case 'scattered'
            R_loc = 1.5e-3;
            if isfield(cfg,'sensor_r_loc') && ~isempty(cfg.sensor_r_loc), R_loc = cfg.sensor_r_loc; end
            [V, dbg] = vmat_scattered(raw, samp, cen, sensor_n, S_hall, N_I, R_loc, dbg);
        case 'grid'
            [V, dbg] = vmat_grid(sfld, sensor_n, S_hall, N_I, dbg);
        otherwise
            error('build_V_matrix: V_METHOD 必為 ''csv-tet'' | ''scattered'' | ''grid''（得 ''%s''）', V_METHOD);
    end

    % ④ all-source：翻下極激發欄（= cfg.s_source 對應激發極；與 build_actuator_data 一致）
    exc_sign = cfg.s_source(cfg.apdl_to_paper_idx);   % 1×N_I
    V = V .* exc_sign;
end

% ---- csv-tet：sensor_local CSV 建 tet、重心內插（需 CSV 與 solve mesh 同網格）----
function V = vmat_csv_tet(cfg, variant, raw, samp, sensor_n, S_hall, N_I)
    csvdir = ansys_path(cfg.model, 'csv', variant);
    Nn = readmatrix(fullfile(csvdir, 'sensor_local_nodes.csv'));   % [id x y z]
    Ee = readmatrix(fullfile(csvdir, 'sensor_local_elems.csv'));   % [eid s1..s8]
    nid = Nn(:,1);  Pn = Nn(:,2:4);
    g2l = zeros(max(nid),1);  g2l(nid) = 1:numel(nid);
    Sl = Ee(:,2:9);  tets = zeros(size(Sl,1),4);  kk = 0;
    for r = 1:size(Sl,1)
        u = unique(Sl(r,:), 'stable');            % SOLID96 8 槽 → tet 4 相異
        if numel(u) == 4, kk = kk+1; tets(kk,:) = g2l(u); end
    end
    tets = tets(1:kk,:);
    v1 = Pn(tets(:,2),:)-Pn(tets(:,1),:);  v2 = Pn(tets(:,3),:)-Pn(tets(:,1),:);  v3 = Pn(tets(:,4),:)-Pn(tets(:,1),:);
    bad = dot(v1, cross(v2,v3,2), 2) < 0;  tets(bad,[3 4]) = tets(bad,[4 3]);   % 正體積
    TR = triangulation(tets, Pn);
    m2  = max(max(nid), max(raw.node_id));
    id2 = zeros(m2,1);  id2(raw.node_id) = 1:numel(raw.node_id);
    li  = zeros(numel(nid),1);  inb = nid <= m2;  li(inb) = id2(nid(inb));
    if any(li == 0)
        error('build_V_matrix: sensor_local 節點 ID 對不上 raw（網格不一致；此變體改用 V_METHOD=''scattered''）');
    end
    V = zeros(6, N_I);
    for kc = 1:N_I
        Bnode = 1e3 * raw.B(li, :, kc);               % sensor_local 節點 B [mT]（T→×1e3）
        for i = 1:6
            pts = samp{i};  ni = sensor_n(:,i);
            ti  = pointLocation(TR, pts);  good = ~isnan(ti);
            bc  = cartesianToBarycentric(TR, ti(good), pts(good,:));
            conn = TR.ConnectivityList(ti(good),:);
            Bp = zeros(nnz(good), 3);
            for c = 1:4, Bp = Bp + bc(:,c).*Bnode(conn(:,c),:); end   % B(p)=Σλ·B
            V(i,kc) = S_hall * mean(Bp * ni);         % ⟨B·n̂⟩[mT]×S_hall[mV/mT] = [mV]
        end
    end
end

% ---- scattered：對 raw solve 場在 sensor 鄰域建 scatteredInterpolant（CSV≠solve mesh 時用）----
function [V, dbg] = vmat_scattered(raw, samp, cen, sensor_n, S_hall, N_I, R_loc, dbg)
    X = [raw.x, raw.y, raw.z];                        % ANSYS 框 [m]
    V = zeros(6, N_I);
    for i = 1:6
        pts = samp{i};  ni = sensor_n(:,i);  ci = cen(:,i);
        rb = R_loc;  ng = vecnorm(X - ci.', 2, 2) < rb;      % 鄰域源節點（不夠就放大）
        while nnz(ng) < 40 && rb < 8e-3, rb = rb*1.5; ng = vecnorm(X - ci.', 2, 2) < rb; end
        Xs = X(ng,:);
        for kc = 1:N_I
            Bk = 1e3 * raw.B(ng,:,kc);               % mT
            FX = scatteredInterpolant(Xs(:,1),Xs(:,2),Xs(:,3),Bk(:,1),'linear','none');
            FY = scatteredInterpolant(Xs(:,1),Xs(:,2),Xs(:,3),Bk(:,2),'linear','none');
            FZ = scatteredInterpolant(Xs(:,1),Xs(:,2),Xs(:,3),Bk(:,3),'linear','none');
            Bp = [FX(pts), FY(pts), FZ(pts)];
            bn = Bp * ni;                            % 逐點 b·n̂⁺ [mT]
            V(i,kc) = S_hall * mean(bn, 'omitnan');
            % [ADDED 2026-08-06] 逐點診斷：∠(b, n̂⁺) = acos( (b·n̂)/‖b‖ )
            bm = vecnorm(Bp, 2, 2);
            dbg.ang{i,kc}  = acosd(min(1, max(-1, bn./bm)));
            dbg.bn{i,kc}   = bn;
            dbg.bmag{i,kc} = bm;
        end
    end
end

% ---- grid：場已由 conv_design_sensor 內插好，本函式**只組 V** ----------------
%   [MODIFIED 2026-08-23 使用者拍板] 原本是呼叫 conv_design_ws 的 query 模式
%   自己內插；現在三線性已併進 conv_design_sensor（與 conv_design_ws 對稱），
%   這裡拿到的 sfld{i} 就是該 sensor 撒點上的場 [mT]（全域框，與撒點 1:1）。
%   盒外檢查也已移進 conv_design_sensor（它 assert 撒點必須全部落在 .fld 盒內）。
function [V, dbg] = vmat_grid(sfld, sensor_n, S_hall, N_I, dbg)
    V = zeros(6, N_I);
    for i = 1:6
        ni = sensor_n(:,i);
        for kc = 1:N_I
            Bp = sfld{i}(:,:,kc);                            % [mT]
            bn = Bp * ni;
            V(i,kc) = S_hall * mean(bn, 'omitnan');
            bm = vecnorm(Bp, 2, 2);
            dbg.ang{i,kc}  = acosd(min(1, max(-1, bn./bm)));
            dbg.bn{i,kc}   = bn;
            dbg.bmag{i,kc} = bm;
        end
    end
end
