function [Vmat, exc_sign] = extract_Vmat_interp(results_root, cnst, apdl_to_paper_idx, ...
                                                sensor_pos, sensor_n, S_hall, mesh_csv_dir, n_uniform, sensor_r, axial_tol, variant)
% EXTRACT_VMAT_INTERP  sensor 電壓矩陣 Vmat（standard 粗網格「真·FEM tet 重心內插」版）。
% -------------------------------------------------------------------------
% 動機：standard(未加密)網格在 sensor 取樣圓柱內 0 真實節點 → 無法直接平均。
%   改在 sensor 圓柱內「均勻撒 n_uniform 個點」，對 standard 場做真·FEM tet 重心內插：
%   每個取樣點落在哪顆原始四面體 → 用該 tet 4 角的重心權重 λ 內插 B(p)=Σλ_iB_i，
%   逐 sensor 對 B·n+ 平均 ×S_hall = sensor 電壓。等價 ANSYS PDEF，但用 .dat、免重解。
%
% 取樣點：各 sensor 圓柱(Ø2·sensor_r × axial_tol)內**均勻撒 n_uniform 個點**（純幾何、
%   rng(0) 可重現；不依賴加密節點位置）。圓柱：底面=sensor 面、沿 n+ 高 axial_tol、半徑 sensor_r。
% 網格連接性：mesh_csv_dir 內 sensor_local_{nodes,elems}.csv（MT_Export_SensorLocalMesh.txt 由
%   mesh_baseline.db 匯出的 6-sensor 局部 SOLID96 tet；只匯網格、無重解）。
% 場：standard coil1..6 的 .dat 節點 B（按 node_id 對齊局部節點）。
% 符號：raw B·n+ 後套 all-source（翻下極激發 P1/P3/P6），與 extract_Vmat 一致。
%
% 輸入同 extract_Vmat，外加：
%   mesh_csv_dir   局部網格 CSV 資料夾
%   n_uniform      (選填) 每 sensor 均勻取樣點數，預設 1000
%   sensor_r/axial_tol (選填) 圓柱半徑/高，預設 0.15e-3 / 0.10e-3（同 extract_Vmat）
% 輸出：Vmat(6×6, [mV]；B 於函式內 ×1e3 轉 mT、S_hall=130 mV/mT)、exc_sign(1×6)。
%   需 import_ansys_data 在 path；用 triangulation/pointLocation（base MATLAB）。
% -------------------------------------------------------------------------
    if nargin < 8 || isempty(n_uniform), n_uniform = 1000;    end
    if nargin < 9 || isempty(sensor_r),  sensor_r  = 0.15e-3; end
    if nargin < 10|| isempty(axial_tol), axial_tol = 0.10e-3; end
    if nargin < 11|| isempty(variant),   variant   = 'standard'; end  % [MODIFIED] 場來源子夾（variant 為第 11 參數；網格 CSV 仍用 standard 拓樸；gap100um 等同網格變體可直接讀）

    % ---- 局部網格 → triangulation（真 FEM tet）----
    N = readmatrix(fullfile(mesh_csv_dir,'sensor_local_nodes.csv'));   % [id x y z]
    E = readmatrix(fullfile(mesh_csv_dir,'sensor_local_elems.csv'));   % [eid s1..s8]
    nid = N(:,1); P = N(:,2:4); mxid = max(nid); g2l = zeros(mxid,1); g2l(nid) = 1:numel(nid);
    S = E(:,2:9); tets = zeros(size(S,1),4); kk = 0;
    for r = 1:size(S,1)
        u = unique(S(r,:),'stable');                                  % SOLID96 8 槽 → tet 4 相異
        if numel(u) == 4, kk = kk+1; tets(kk,:) = g2l(u); end
    end
    tets = tets(1:kk,:);
    v1 = P(tets(:,2),:)-P(tets(:,1),:); v2 = P(tets(:,3),:)-P(tets(:,1),:); v3 = P(tets(:,4),:)-P(tets(:,1),:);
    vol = dot(v1,cross(v2,v3,2),2); bad = vol < 0; tets(bad,[3 4]) = tets(bad,[4 3]);   % 正體積
    TR = triangulation(tets, P);

    % ---- 取樣點：各 sensor 圓柱內均勻撒 n_uniform 點（純幾何、ANSYS 框、rng 可重現）----
    rng(0);
    samp = cell(1,6);
    for i = 1:6
        ci = sensor_pos(:,i) + [0;0;cnst.SPH_OFST];     % sensor 中心(ANSYS 框)
        ni = sensor_n(:,i);
        t1 = [-ni(2); ni(1); 0]; if norm(t1) < 1e-9, t1 = [1;0;0]; end
        t1 = t1/norm(t1); t2 = cross(ni, t1);            % 圓柱橫切基底 ⊥ n+
        a  = axial_tol * rand(n_uniform,1);              % 軸向 U[0,H]（沿 n+）
        r  = sensor_r  * sqrt(rand(n_uniform,1));        % 徑向（√U → 面積均勻）
        th = 2*pi      * rand(n_uniform,1);              % 方位 U[0,2π]
        samp{i} = ci.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';  % n_uniform×3
    end

    % ---- 逐 coil：standard B → 對齊局部節點 → 逐取樣點重心內插 → B·n+ 平均 ----
    Vmat = zeros(6,6);
    for kc = 1:6
        cn = sprintf('coil%d', kc);
        ds = import_ansys_data(fullfile(results_root, cn, variant),'all',cn);   % [MODIFIED] standard→variant（同網格、只換場來源）
        m2 = max(max(nid), max(ds.node_id)); id2 = zeros(m2,1); id2(ds.node_id) = 1:numel(ds.node_id);
        li = zeros(numel(nid),1); inb = nid <= m2; li(inb) = id2(nid(inb));
        if any(li == 0)
            error('extract_Vmat_interp: 局部節點 ID 對不上 standard .dat（網格不一致）。');
        end
        Bnode = 1e3*[ds.bx(li), ds.by(li), ds.bz(li)];                % 每局部節點 B：Tesla → ×1e3 原生 mT（Unit Sheet）
        for i = 1:6
            pts = samp{i}; ni = sensor_n(:,i);
            ti = pointLocation(TR, pts); good = ~isnan(ti);
            bc = cartesianToBarycentric(TR, ti(good), pts(good,:));
            conn = TR.ConnectivityList(ti(good),:);
            Bp = zeros(nnz(good),3);
            for c = 1:4, Bp = Bp + bc(:,c).*Bnode(conn(:,c),:); end    % B(p)=Σλ_iB_i
            Vmat(i,kc) = S_hall * mean(Bp*ni);                        % ⟨B·n+⟩[mT]×S_hall[mV/mT] = sensor 電壓 [mV]
            if kc == 1
                fprintf('  [interp] sensor P%d: %d/%d 點找到包含 tet\n', i, nnz(good), numel(good));
            end
        end
        fprintf('  [interp] coil%d 內插完成\n', kc);
    end

    % ---- all-source：翻下極激發欄（P1/P3/P6）----
    exc_sign = ones(1,6);
    for j = 1:6
        if ismember(apdl_to_paper_idx(j), [1 3 6]), exc_sign(j) = -1; end
    end
    Vmat = Vmat .* exc_sign;
end
