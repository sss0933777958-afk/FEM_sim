function D = load_coils_actuator(model, cnst, apdl_to_paper_idx, dataset, variant)
%LOAD_COILS_ACTUATOR  NTU：載 full_assembly 6-coil 場 → WP-centered world 框（NTU 扁平板，無 magic-angle 旋轉）。
%   讀 data/full_assembly/coilN/full_assembly_*_all.dat（檔名前綴=variant='full_assembly'），
%   平移到 WP（cnst.WP）為原點、留 WP 附近節點、NEGATE B（all-source legacy），回：
%     .Pa      Nair×3   WP-centered 取樣點 [m]（不旋轉）
%     .r2      Nair×1   |p−WP|²（select_ball 用）
%     .Ba      Nair×3×6 all-source B [mT]（ANSYS T ×1e3；world 框，不旋轉）
%     .R_act   3×3      =I（NTU 無 actuator 旋轉；欄位保留供下游相容）
%     .Pc_base 3×6      單位 tip 方向（WP-centered，非正交；NTU 無 magic-angle）
%     .F       6×N_I    電流置換矩陣
%   ⚠ NTU 與 long2016 不同：不做 R_act 旋轉、不用 SPH_OFST；WP=6 tip 形心；不 filter_iron
%     （select_ball 取 WP 小球 < tip 距 0.515mm → 皆空氣）。
    if nargin < 5 || isempty(variant), variant = 'full_assembly'; end
    N_I = 6;
    WP = cnst.WP(:).';                                       % 1×3 [m]
    tipwp   = [cnst.pole_tip_x_wp; cnst.pole_tip_y_wp; cnst.pole_tip_z_wp];  % 3×6
    Pc_base = tipwp ./ vecnorm(tipwp);                       % 單位 tip 方向（3×6）
    F = zeros(6, N_I);
    for j = 1:N_I, F(apdl_to_paper_idx(j), j) = 1; end

    R_LOAD = 2e-3;                                           % 只留 WP 半徑 2mm 內節點（省記憶體）
    d1 = import_ansys_data(ansys_path(model,'data',variant,'coil1'), dataset, variant);   % 前綴=full_assembly
    P1  = [d1.x, d1.y, d1.z] - WP;                           % WP-centered
    sel = sum(P1.^2, 2) < R_LOAD^2;
    P_meas = P1(sel, :);
    r2 = sum(P_meas.^2, 2);

    Nair = size(P_meas, 1);
    Ba = zeros(Nair, 3, N_I);
    for k = 1:N_I
        if k == 1, dk = d1;
        else, cn = sprintf('coil%d', k); dk = import_ansys_data(ansys_path(model,'data',variant,cn), dataset, variant);
        end
        Bk = -1e3 * [dk.bx(sel), dk.by(sel), dk.bz(sel)];   % all-source negate；T→mT
        Ba(:,:,k) = Bk;                                     % world 框（不旋轉）
    end

    D = struct('Pa',P_meas, 'r2',r2, 'Ba',Ba, 'R_act',eye(3), 'Pc_base',Pc_base, 'F',F, 'WP',WP);
end
