function [V, exc_sign] = build_V_matrix(cfg, variant, raw, S_hall, SOFF_upper, n_uniform, sensor_r, axial_tol, sensor_override, V_METHOD)
%BUILD_V_MATRIX  電壓提取：sensor 幾何 → 每 sensor 圓柱撒點 → 內插 raw 場 → 6×6 V[mV]。
%   [V, exc_sign] = BUILD_V_MATRIX(cfg, variant, raw, S_hall, SOFF_upper, n_uniform, sensor_r, axial_tol, sensor_override, V_METHOD)
%   sensor 幾何來源優先序：
%     ① sensor_override（可選 struct .pos/.n，escape hatch）
%     ② cfg.sensor_pos / cfg.sensor_n（config「提供」的幾何，如 tip400um）
%     ③ utils/pole_sensor_geometry（**sensor 幾何唯一來源**；CAD 實測錐體 + 真實氣隙）
%   V_METHOD：'csv-tet'（預設；讀 sensor_local CSV 建 tet 重心內插——需 CSV 與 solve mesh 同網格）
%             'scattered'（座標式 scatteredInterpolant 直接對 raw solve 場取樣——CSV≠solve mesh 時用，如 tip400）。
%   需 ansys_path + readmatrix + triangulation/pointLocation（csv-tet）/ scatteredInterpolant（scattered）在 path。
    if nargin < 5  || isempty(SOFF_upper), SOFF_upper = 4.572e-3; end
    if nargin < 6  || isempty(n_uniform),  n_uniform  = 1e4;      end
    if nargin < 7  || isempty(sensor_r),   sensor_r   = 0.15e-3;  end
    if nargin < 8  || isempty(axial_tol),  axial_tol  = 0.10e-3;  end
    if nargin < 9,  sensor_override = []; end
    if nargin < 10 || isempty(V_METHOD),   V_METHOD   = 'csv-tet'; end
    N_I = cfg.N_I;

    % ① sensor 幾何（WP 框）：override > config 提供 > 內建 baseline 錐面
    if ~isempty(sensor_override)
        sensor_pos = sensor_override.pos;  sensor_n = sensor_override.n;
    elseif isfield(cfg,'sensor_pos') && isfield(cfg,'sensor_n') && ~isempty(cfg.sensor_pos)
        sensor_pos = cfg.sensor_pos;  sensor_n = cfg.sensor_n;
    else
        % [MODIFIED 2026-08-08] 改呼叫 utils/pole_sensor_geometry（sensor 幾何唯一來源），
        %   原本的 local sensor_geometry 已移除。下極一律貼底錐面（本分支無 FACE_lower）。
        [sensor_pos, sensor_n] = pole_sensor_geometry(cfg, struct('soff_upper', SOFF_upper));
    end

    % ② 每 sensor 圓柱內均勻撒 n_uniform 點（ANSYS 框、rng 可重現）——兩法共用
    rng(0);
    cen = zeros(3,6);  samp = cell(1,6);
    for i = 1:6
        ci = sensor_pos(:,i) + [0;0;cfg.SPH_OFST];    % WP 框 → ANSYS 框
        cen(:,i) = ci;
        ni = sensor_n(:,i);
        t1 = [-ni(2); ni(1); 0]; if norm(t1) < 1e-9, t1 = [1;0;0]; end
        t1 = t1/norm(t1);  t2 = cross(ni, t1);
        a  = axial_tol * rand(n_uniform,1);
        rr = sensor_r  * sqrt(rand(n_uniform,1));     % √U → 面積均勻
        th = 2*pi      * rand(n_uniform,1);
        samp{i} = ci.' + a.*ni.' + (rr.*cos(th)).*t1.' + (rr.*sin(th)).*t2.';
    end

    % ③ 內插 raw 場到撒點 → V=S_hall·⟨B·n̂⟩[mV]（依 V_METHOD 分兩法）
    switch V_METHOD
        case 'csv-tet'
            V = vmat_csv_tet(cfg, variant, raw, samp, sensor_n, S_hall, N_I);
        case 'scattered'
            R_loc = 1.5e-3;
            if isfield(cfg,'sensor_r_loc') && ~isempty(cfg.sensor_r_loc), R_loc = cfg.sensor_r_loc; end
            V = vmat_scattered(raw, samp, cen, sensor_n, S_hall, N_I, R_loc);
        otherwise
            error('build_V_matrix: V_METHOD 必為 ''csv-tet'' | ''scattered''（得 ''%s''）', V_METHOD);
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
function V = vmat_scattered(raw, samp, cen, sensor_n, S_hall, N_I, R_loc)
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
            V(i,kc) = S_hall * mean(Bp * ni, 'omitnan');
        end
    end
end
