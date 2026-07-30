function Bstack_t = interp_field_to_points(cfg, variant, R_act, P_target, R_loc)
%INTERP_FIELD_TO_POINTS  把某 FEM 變體的 6-coil WP 場，內插到指定 actuator-frame 點上。
%   Bstack_t = INTERP_FIELD_TO_POINTS(cfg, variant, R_act, P_target, R_loc)
%   用途：tip 變體公平比較 —— 把 <variant>（如 'tip400um'）的場，取在基準變體（tip40/graded）
%   自己的 R≤R_select 節點座標 P_target 上，達成「同數量、同分布」的取樣，再共用同一套 fit / solve。
%
%   輸入：
%     cfg      : model_config 輸出（提供 SPH_OFST / s_source / apdl_to_paper_idx / N_I）
%     variant  : FEM 變體子夾（data/<variant>/coilN），例 'tip400um'
%     R_act    : 3×3 actuator 旋轉（用「基準變體」的；tip400 的尖端到尖端方向 = baseline magic-angle，
%                故 tip400 與 tip40 同一 actuator frame，可直接在此 frame 內插）
%     P_target : Np×3 查詢點（actuator frame [m]）；通常 = 基準變體 select_ball 的 P
%     R_loc    : 內插「源節點」鄰域球半徑 [m]（須 > 查詢半徑以圍住 P_target；預設 0.6e-3）
%   輸出：
%     Bstack_t : 3Np×N_I  堆疊 [Bx;By;Bz;...] per coil（**all-source mT**，與 select_ball 的
%                Bstack 同慣例：×1e3 T→mT + cfg.s_source 翻下極 sink；下游不必再翻號）
%
%   做法：讀 <variant> raw（measure/Tesla）→ WP frame → 取 WP 鄰域源節點（tip400 尖端退到 ~1.9mm、
%   0.6mm 球內全空氣、免濾鐵）→ 套 all-source → R_act 轉 actuator → Bx/By/Bz 各建
%   scatteredInterpolant('linear','none' 不外插) → 查 P_target。R≤150µm 全在源網格內部、正常無 NaN。
%   需 extract_ansys_data + ansys_path + import_ansys_data 在 path。2026-07-22。
    if nargin < 5 || isempty(R_loc), R_loc = 0.6e-3; end

    raw  = extract_ansys_data(cfg, 'all', variant);
    zwp  = raw.z - cfg.SPH_OFST;
    Pw   = [raw.x, raw.y, zwp];                       % WP frame [m]
    near = sum(Pw.^2, 2) < R_loc^2;                   % WP 鄰域源節點（同 select_ball 的 r2 判準）
    Xs   = (R_act * Pw(near,:).').';                  % → actuator frame

    col_sign = cfg.s_source(cfg.apdl_to_paper_idx);   % all-source：翻下極 sink（與 build_actuator_data 一致）
    Np  = size(P_target, 1);
    N_I = cfg.N_I;
    Bstack_t = zeros(3*Np, N_I);
    for j = 1:N_I
        Bk = (1e3 * col_sign(j)) * raw.B(near,:,j);   % T→mT + all-source
        Ba = (R_act * Bk.').';                         % → actuator frame
        FX = scatteredInterpolant(Xs(:,1), Xs(:,2), Xs(:,3), Ba(:,1), 'linear', 'none');
        FY = scatteredInterpolant(Xs(:,1), Xs(:,2), Xs(:,3), Ba(:,2), 'linear', 'none');
        FZ = scatteredInterpolant(Xs(:,1), Xs(:,2), Xs(:,3), Ba(:,3), 'linear', 'none');
        Bq = [FX(P_target), FY(P_target), FZ(P_target)];   % Np×3
        Bstack_t(:,j) = reshape(Bq.', [], 1);              % [Bx;By;Bz;...]（同 select_ball）
    end

    nn = nnz(any(isnan(reshape(Bstack_t(:,1), 3, [])).', 2));
    if nn > 0
        warning('interp_field_to_points:hullMiss', ...
                '%d / %d 查詢點落在 %s 網格 convex hull 外（NaN）；呼叫端請剔除。', nn, Np, variant);
    end
end
