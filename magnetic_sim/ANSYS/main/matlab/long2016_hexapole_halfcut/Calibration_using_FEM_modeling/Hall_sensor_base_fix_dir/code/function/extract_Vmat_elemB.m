function res = extract_Vmat_elemB(cnst, results_root, variant)
% EXTRACT_VMAT_ELEMB  sensor 電壓抽取（未平滑「元素常數 B」法）。
% -------------------------------------------------------------------------
% 【建構中 — 依使用者逐步描述累加；目前進度：第 4 步】
%
% 目標（最終）：依序 6 個感測器，在各自「感測面 R0.15mm 圓盤(無高度)」上取樣，
%   用已匯出的「未平滑元素常數 B」(每個一階 tet 一個常數) 算每顆 sensor 的電壓，
%   組成 Vmat（與 extract_Vmat / extract_Vmat_interp 並列、做比較）。
%   資料來源：data/coilN/<variant>/coilN_{elemB,elemconn,nodes}_R05.dat（R0.5mm 球、ANSYS 框）。
%
% 第 1 步：對每個 sensor，從圓心 + 法向量建「與感測面平行的無限平面」(ax+by+cz+d=0)。
% 第 2 步：粗篩——R0.5mm 球內每 tet 的 4 節點代入平面式，符號全同=沒切過(丟)、
%         混合=切過(留)；coil 無關(同網格)，算一次。
% 第 3 步：被切 tet 的 6 邊，端點異號=被切，線性內插求交點；3~4 點排序成截面多邊形
%         (投到平面 2D 基底、繞質心角度排)。
% 第 4 步：感測面本地 2D 座標——參考 u_hat=[1 0 0](取其平面內分量正規化)、v_hat=n×u；
%         截面頂點投到 (u,v)，圓盤圓心在 (0,0)。res.cross2d 即此 2D 多邊形。
%
% 輸入：
%   cnst          幾何常數（mt_constants()）。
%   results_root  (選填) data 根目錄，預設 canonical（…/ANSYS_data/<model>/data）。
%   variant       (選填) FEM 變體子夾，預設 'gap200um_mueq'。
% 輸出（暫時，隨步驟成長）：
%   res.plane       4×6   各 sensor 平面 [a;b;c;d]（WP 框）
%   res.sensor_pos  3×6   各 sensor 圓心 [m]（WP 框）
%   res.sensor_n    3×6   各 sensor 法線 n+（單位）
%   res.cut         1×6 cell  各 sensor「被平面切過的 tet eid」清單（第2步粗篩）
%   res.conn        Map  eid → [n1 n2 n3 n4]（區域 tet 連接性，coil 無關）
%   res.node_wp     Map  nid → [x y z]（WP 框 [m]）
% -------------------------------------------------------------------------
    if nargin < 2 || isempty(results_root)
        results_root = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\' ...
                        'ANSYS_data\long2016_hexapole_halfcut\data'];
    end
    if nargin < 3 || isempty(variant), variant = 'gap200um_mueq'; end

    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);   % 6 sensor 圓心 + 法線 n+
    NS = size(sensor_pos, 2);                               % = 6

    %% ===== 第 1 步：各 sensor 建感測面平行無限平面 =====
    plane = zeros(4, NS);
    for i = 1:NS
        plane(:,i) = make_plane(sensor_pos(:,i), sensor_n(:,i));
    end

    %% ===== 載區域幾何（coil1；連接性/節點 coil 無關）=====
    gdir = fullfile(results_root, 'coil1', variant);
    node_wp = load_nodes_wp(fullfile(gdir,'coil1_nodes_R05.dat'), cnst.SPH_OFST);  % nid → xyz(WP,m)
    conn    = load_conn(fullfile(gdir,'coil1_elemconn_R05.dat'));                  % eid → [4 nodes]

    %% ===== 第 2+3 步：粗篩被切 tet（第2）+ 求截面交點(內插)並排序成多邊形（第3）=====
    RSPH = 0.5e-3;                                          % 匯出球半徑（= ESEL 用值）
    EDG  = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];                 % tet 6 條邊（頂點對）
    eids = cell2mat(keys(conn));                           % 區域所有 tet eid
    cut     = cell(1, NS);                                 % 各 sensor 被切 tet eid
    cross   = cell(1, NS);                                 % 截面多邊形（排序、3D、WP）
    cross2d = cell(1, NS);                                 % 截面多邊形（2D 平面座標、給裁切）
    basis   = cell(1, NS);                                 % 平面 2D 基底
    for i = 1:NS
        pl = plane(:,i);  ci = sensor_pos(:,i);
        [u, v] = plane_basis(pl(1:3));                     % 平面內正交基底 u,v ⊥ n
        basis{i} = struct('u', u, 'v', v, 'o', ci);
        keep = zeros(1, numel(eids));  pc = {};  pc2 = {};  nk = 0;
        for e = eids
            nd = conn(e);
            X  = [node_wp(nd(1)); node_wp(nd(2)); node_wp(nd(3)); node_wp(nd(4))];  % 4×3 WP
            cc = mean(X, 1).';
            if norm(cc - ci) > RSPH, continue; end         % 限該 sensor 的球
            s = X * pl(1:3) + pl(4);                       % 4 節點平面值
            if all(s > 0) || all(s < 0), continue; end     % 同號 = 沒被切 → 跳（第2步）

            % --- 第 3 步：6 邊逐一找交點（異號邊，線性內插）---
            P = zeros(0, 3);
            for k = 1:6
                a = EDG(k,1);  b = EDG(k,2);
                if s(a)*s(b) < 0                            % 異號 = 該邊被切
                    t = s(a) / (s(a) - s(b));              % 內插參數 ∈(0,1)
                    P(end+1,:) = (1-t)*X(a,:) + t*X(b,:);  %#ok<AGROW> 交點座標
                end
            end
            % --- 排序成多邊形：投到 2D(u,v)、繞質心 atan2 角度排（CCW）---
            P2  = [(P - ci.')*u, (P - ci.')*v];            % Npt×2 平面座標
            cen = mean(P2, 1);
            [~, ord] = sort(atan2(P2(:,2)-cen(2), P2(:,1)-cen(1)));
            P = P(ord,:);  P2 = P2(ord,:);

            nk = nk + 1;  keep(nk) = e;  pc{nk} = P;  pc2{nk} = P2; %#ok<AGROW>
        end
        cut{i} = keep(1:nk);  cross{i} = pc;  cross2d{i} = pc2;
    end

    %% ===== 第 4 步：(待續) 圓盤≈64-gon ∩ 截面多邊形 → 面積 × 元素 B·n+ 加權平均 → 電壓 =====

    res = struct('plane', plane, 'sensor_pos', sensor_pos, 'sensor_n', sensor_n, ...
                 'cut', {cut}, 'cross', {cross}, 'cross2d', {cross2d}, 'basis', {basis}, ...
                 'conn', conn, 'node_wp', node_wp);
end

% ---- 感測面本地 2D 基底（參考 u_hat=[1 0 0]、皆 ⊥ n、單位）----
%   v_hat = normalize(n × [1 0 0])（使用者式 v=n×u）；
%   u_hat = normalize(v_hat × n) = [1 0 0] 在平面內的分量（正規化）。
%   注：[1 0 0] 一般不⊥n，故取其平面內分量當 u 軸，確保等距(圓盤維持圓、面積正確)。
function [u, v] = plane_basis(n)
    n = n(:) / norm(n);
    uref = [1;0;0];
    v = cross(n, uref);
    if norm(v) < 1e-6, uref = [0;1;0]; v = cross(n, uref); end  % n ∥ x 退路
    v = v / norm(v);                       % v_hat = normalize(n × [1 0 0])
    u = cross(v, n);  u = u / norm(u);     % u_hat = [1 0 0] 平面內分量（⊥v、⊥n、單位）
end

% ---- 第 1 步：點 + 法向量 → 平面通式 ----
function p = make_plane(center, normal)
    n = normal(:) / norm(normal);          % 單位法向量 = 平面 [a;b;c]
    d = -dot(n, center(:));                % 平面常數（過 center）
    p = [n; d];
end

% ---- 載 NLIST 節點檔 → containers.Map: nid → [x y z](WP, m) ----
function M = load_nodes_wp(fpath, sph_ofst)
    fid = fopen(fpath,'r');  raw = fread(fid,'*char').';  fclose(fid);
    lines = regexp(raw, '\r?\n', 'split');
    ids = []; XYZ = [];
    for k = 1:numel(lines)
        t = sscanf(lines{k}, '%f');                       % 整數 + 3 浮點的行才會湊滿 4 個
        if numel(t) >= 4
            ids(end+1,1) = round(t(1)); %#ok<AGROW>
            XYZ(end+1,:) = [t(2), t(3), t(4) - sph_ofst]; %#ok<AGROW>  z → WP
        end
    end
    M = containers.Map('KeyType','double','ValueType','any');
    for k = 1:numel(ids), M(ids(k)) = XYZ(k,:); end
end

% ---- 載 CSV 連接性檔 → containers.Map: eid → [n1 n2 n3 n4]（unique tet）----
function M = load_conn(fpath)
    A = readmatrix(fpath);                                % 逗號分隔；col1=eid, col2-9=8 node 槽
    M = containers.Map('KeyType','double','ValueType','any');
    for r = 1:size(A,1)
        eid = round(A(r,1));
        slots = round(A(r,2:9));
        u = unique(slots, 'stable');  u = u(u>0);
        if numel(u) == 4, M(eid) = u; end
    end
end
