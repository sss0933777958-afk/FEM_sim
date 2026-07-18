function [Vmat, exc_sign] = build_V_matrix(results_root, cnst, apdl_to_paper_idx, ...
                                          sensor_pos, sensor_n, S_hall, mesh_csv_dir, n_uniform, sensor_r, axial_tol, variant, coils)
%BUILD_V_MATRIX  NTU：6 sensor × 6 coil 電壓矩陣 V̄（實節點 Delaunay 重心內插版）。
% -------------------------------------------------------------------------
% NTU full_assembly 網格在每顆 sensor 已加密（6 顆 R0.16mm air 球，~500 節點/sensor）→ 直接用
% **實 FEM 節點建 Delaunay**（比照 field_viz/sensor_average），不需匯出 CSV tet 拓樸。
%   每 coil 讀 data/full_assembly/coilN/full_assembly_*_all.dat（world/ANSYS 框，SI）；
%   每 sensor：取中心 0.30mm 內真實節點 → delaunayTriangulation → 圓柱內面積均勻撒 n_uniform 點
%   → pointLocation + 重心內插 B(p)=Σλ_iB_i → 平均 B·n+ ×S_hall = sensor 電壓 [mV]。
%   sensor_pos 已是 world 框（無 SPH_OFST）；all-source 翻下極激發欄 P1/P3/P6。
% 輸出：Vmat(6×6,[mV])、exc_sign(1×6)。（mesh_csv_dir 參數保留相容、NTU 不用。）
% -------------------------------------------------------------------------
    if nargin < 8 || isempty(n_uniform), n_uniform = 1000;        end
    if nargin < 9 || isempty(sensor_r),  sensor_r  = cnst.SENSOR_R; end
    if nargin < 10|| isempty(axial_tol), axial_tol = cnst.SENSOR_H; end
    if nargin < 11|| isempty(variant),   variant   = 'full_assembly'; end
    if nargin < 12|| isempty(coils),     coils     = 1:6;          end
    LOC  = 0.30e-3;                                   % 局部節點球（建 Delaunay）
    OVER = 2*n_uniform;                              % 過取樣（凸包外 NaN 剔除後仍夠 n_uniform）
    rng(1);

    Vmat = zeros(6,6);
    for kc = coils
        cn = sprintf('coil%d', kc);
        ds = import_ansys_data(fullfile(results_root, variant, cn), 'all', variant);   % 前綴=full_assembly；SI
        X = [ds.x ds.y ds.z];  Bn3 = [ds.bx ds.by ds.bz];
        for i = 1:6
            c = sensor_pos(:,i);  ni = sensor_n(:,i);
            loc = vecnorm(X - c.', 2, 2) < LOC;                     % 中心 0.30mm 內真實節點
            P = X(loc,:);  BP = Bn3(loc,:);
            TR = delaunayTriangulation(P);
            Tb = null(ni.');  t1 = Tb(:,1);  t2 = Tb(:,2);          % 圓柱橫切基底 ⊥ n+
            a  = axial_tol*rand(OVER,1);  r = sensor_r*sqrt(rand(OVER,1));  th = 2*pi*rand(OVER,1);   % 面積均勻
            pts = c.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';
            ti = pointLocation(TR, pts);  good = find(~isnan(ti));
            good = good(1:min(n_uniform, numel(good)));             % 取前 n_uniform 個有效點
            bc = cartesianToBarycentric(TR, ti(good), pts(good,:));
            conn = TR.ConnectivityList(ti(good),:);
            Bint = zeros(numel(good),3);
            for cc = 1:4, Bint = Bint + bc(:,cc).*BP(conn(:,cc),:); end
            Bproj = (Bint*ni) * 1e3;                                % Tesla → mT
            Vmat(i,kc) = S_hall * mean(Bproj);                     % ⟨B·n+⟩[mT]×S_hall[mV/mT] = [mV]
            if kc == coils(1)
                fprintf('  [interp] sensor P%d: %d local nodes, %d/%d pts\n', i, nnz(loc), numel(good), OVER);
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
