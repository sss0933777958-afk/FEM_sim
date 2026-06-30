function [Vmat, exc_sign, Bn_raw] = extract_Vmat_interp_center(results_root, cnst, ...
                                        apdl_to_paper_idx, sensor_pos, sensor_n, S_hall, mesh_csv_dir, variant)
% EXTRACT_VMAT_INTERP_CENTER  sensor 電壓矩陣（standard 粗網格「單點重心內插」版）。
% -------------------------------------------------------------------------
% 取值方式（使用者拍板）：sensor 範圍仍是半徑 0.15mm、高 0.1mm 的圓柱，但**不再平均圓柱內
%   節點**，改成在「圓柱底面中心」用真·FEM tet 重心內插**單一點**的 B·n+。
%   底面 = 從磁極沿 n+ 方向 0.41mm（air-gap）第一個碰到的面；其中心 = sensor_pos（離極 0.41mm）。
%   評估點(ANSYS 框) = sensor_pos(:,i) + [0;0;SPH_OFST]。
%
% 內插：點落在哪顆原始 SOLID96 四面體 → 用該 tet 4 角重心權重 λ 內插 B(p)=Σλ_iB_i。
%   連接性取自 mesh_csv_dir 的 sensor_local_{nodes,elems}.csv（mesh_baseline.db 匯出、免重解）；
%   場取 standard coil1..6 的 .dat 節點 B（按 node_id 對齊）。
% 符號：raw B·n+ 後套 all-source（翻下極激發 P1/P3/P6），與 extract_Vmat 一致。
%
% 輸出：Vmat(6×6,[V])、exc_sign(1×6)、Bn_raw(6×6,[T]) = 內插單點 raw B·n+（未乘 S、未翻號）。
%   需 import_ansys_data 在 path；用 triangulation/pointLocation（base MATLAB）。
%   variant (選填)：FEM 結果子夾名，預設 'standard'（baseline）；'gap200um_mueq' 等變體同網格可直接用。
% -------------------------------------------------------------------------
    if nargin < 8 || isempty(variant), variant = 'standard'; end

    % ---- 局部網格 → triangulation（真 FEM tet；同 extract_Vmat_interp）----
    N = readmatrix(fullfile(mesh_csv_dir,'sensor_local_nodes.csv'));   % [id x y z]
    E = readmatrix(fullfile(mesh_csv_dir,'sensor_local_elems.csv'));   % [eid s1..s8]
    nid = N(:,1); P = N(:,2:4); mxid = max(nid); g2l = zeros(mxid,1); g2l(nid) = 1:numel(nid);
    Sl = E(:,2:9); tets = zeros(size(Sl,1),4); kk = 0;
    for r = 1:size(Sl,1)
        u = unique(Sl(r,:),'stable');                                 % SOLID96 8 槽 → tet 4 相異
        if numel(u) == 4, kk = kk+1; tets(kk,:) = g2l(u); end
    end
    tets = tets(1:kk,:);
    v1 = P(tets(:,2),:)-P(tets(:,1),:); v2 = P(tets(:,3),:)-P(tets(:,1),:); v3 = P(tets(:,4),:)-P(tets(:,1),:);
    vol = dot(v1,cross(v2,v3,2),2); bad = vol < 0; tets(bad,[3 4]) = tets(bad,[4 3]);   % 正體積
    TR = triangulation(tets, P);

    % ---- 評估點：各 sensor 圓柱底面中心（= sensor_pos，ANSYS 框）----
    ctr = zeros(3,6);
    for i = 1:6, ctr(:,i) = sensor_pos(:,i) + [0;0;cnst.SPH_OFST]; end

    % ---- 逐 coil：standard B → 對齊局部節點 → 各 sensor 單點重心內插 ----
    Vmat = zeros(6,6); Bn_raw = zeros(6,6);
    for kc = 1:6
        cn = sprintf('coil%d', kc);
        ds = import_ansys_data(fullfile(results_root, cn, variant),'all',cn);
        m2 = max(max(nid), max(ds.node_id)); id2 = zeros(m2,1); id2(ds.node_id) = 1:numel(ds.node_id);
        li = zeros(numel(nid),1); inb = nid <= m2; li(inb) = id2(nid(inb));
        if any(li == 0)
            error('extract_Vmat_interp_center: 局部節點 ID 對不上 standard .dat（網格不一致）。');
        end
        Bnode = [ds.bx(li), ds.by(li), ds.bz(li)];                    % 每局部節點的 standard B
        for i = 1:6
            p = ctr(:,i).'; ni = sensor_n(:,i);
            ti = pointLocation(TR, p);                                % 含此點的 tet
            if isnan(ti)
                error('extract_Vmat_interp_center: sensor P%d 底面中心不在局部網格內（擴大匯出框）。', i);
            end
            bc   = cartesianToBarycentric(TR, ti, p);                 % 1×4 重心權重
            conn = TR.ConnectivityList(ti,:);
            Bp   = bc * Bnode(conn,:);                                % 1×3：B(p)=Σλ_iB_i
            Bn_raw(i,kc) = Bp * ni;                                   % raw 單點 B·n+ [T]
            Vmat(i,kc)   = S_hall * Bn_raw(i,kc);
        end
        fprintf('  [interp-center] coil%d 6 顆 sensor 單點內插完成\n', kc);
    end

    % ---- all-source：翻下極激發欄（P1/P3/P6）----
    exc_sign = ones(1,6);
    for j = 1:6
        if ismember(apdl_to_paper_idx(j), [1 3 6]), exc_sign(j) = -1; end
    end
    Vmat   = Vmat   .* exc_sign;
    Bn_raw = Bn_raw .* exc_sign;     % 也回傳 all-source 後的單點 B·n+（方便看 sign）
end
