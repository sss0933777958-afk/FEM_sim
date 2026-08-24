function c = mt_constants()
%MT_CONSTANTS  zhi_peng 平面六極 **R500 變體**的 per-model 設定（**Maxwell** 求解器）。
%   幾何來源（CAD = source of truth）：
%     CAD_model/zhi-Peng/STEP/zhi_peng_R500.STEP
%   （純解析幾何 PLANE/CYLINDRICAL_SURFACE/CIRCLE，無 B-spline → OCC 直接讀面即可，
%     不受「控制點 ≠ 曲面邊界」的陷阱影響，per step-geom-extraction。）
%
%   結構比照 config/long2016_hexapole_halfcut/tip40um/mt_constants.m。
%
%   ⚠ 與 R594 變體的差別：工作半徑 594 → 500 µm；連帶極板變薄（686 → 577.35 µm）、
%     導柱位置外移（±20.825 → ±21.5 mm）。極尖幾何（40 µm 倒圓 / 178 µm 舌片厚 / 7.5° 楔角）**不變**。
%
%   ⚠ 與 long2016 的差別（別互抄）：
%     ① 磁極是**平面板**（舌片厚 178 µm、尖端 40 µm 倒圓、楔角 7.5°），不是圓錐
%        → 無 POLE_R / POLE_CONE_LEN，filter_iron_nodes 會偵測到缺欄位而跳過錐體判定（見檔尾）。
%     ② 六極**共平面製作**：下極尖在板下緣、上極尖在板上緣，工作區中心在中面。
%        雖然做法不同，magic-angle 條件仍成立（見下）。
%     ③ **極尖是垂直圓角「邊」而非「點」**（沿 z 長 178 µm）→ 最近鐵不在極尖，見 IRON_SAFE_R。

    % ===== ① 提取識別 / 資料路由 =====
    c.strategy          = 'hex_magic';               % magic-angle 六極：R_act 旋轉 + 濾鐵
    c.apdl_to_paper_idx = [1, 2, 3, 4, 5, 6];        % identity：.fld 檔名 B_p<k> = paper P<k>
                                                     % （新 model 一律 identity，per pole-coil-numbering 規則）
    c.default_variant   = 'maxwell';
    c.regions           = {'all'};
    c.R_load            = [];

    % .fld 六檔已匯出並驗證（2026-08-17 20:05–23:00，各 ~1.15 GB）：
    %   header 逐字一致：Min [-2 -2 -1.711] mm  Max [2 2 2.289] mm  Grid [0.02]³ mm
    %   → 201³ = 8,120,601 格點/檔；z 框中心 (-1.711+2.289)/2 = 0.289 mm，與 WP_Z_CAD 288.675 µm 相符。
    %   六檔在同一格點的 B 互不相同（已比對）→ 確實是六組不同激發，非重複匯出。
    %   來源專案：D:\Maxwell_sim\Zhi_peng\project\R500\R500.aedt
    c.fld_dir           = 'D:\Maxwell_sim\Zhi_peng\export\R500';
    c.fld_files         = {'B_p1.fld','B_p2.fld','B_p3.fld','B_p4.fld','B_p5.fld','B_p6.fld'};
    % [ADDED 2026-08-21] variant 路由（沿用 long2016 的 fld_files_variant 機制）：
    %   'maxwell'       = 2026-08-17 首次匯出（config 預設）
    %   'maxwell_split' = 2026-08-20 重解版（磁極已切開 + 網格加密），
    %                     **同框同步距 201^3**，但場值與舊版不同（已逐行比對確認）。
    c.fld_files_variant.maxwell       = c.fld_files;
    c.fld_files_variant.maxwell_split = {'B_p1_split.fld','B_p2_split.fld','B_p3_split.fld', ...
                                         'B_p4_split.fld','B_p5_split.fld','B_p6_split.fld'};
    % [ADDED 2026-08-24] 'maxwell_gap' = 氣隙版重解（匯在**平行資料夾** export\R500_gap\），
    %   用該夾的 *_split.fld（WP 細格，與上面兩版同框同步距：
    %   Min [-2 -2 -1.711] mm / Max [2 2 2.289] mm / Grid [0.02]^3 → 201^3 格點，header 逐字一致）。
    %   同夾另有 *_glob（±27 mm / 0.2 mm 全機粗格，另一種用途，未路由；且 p6 匯出不完整）。
    c.fld_dir_variant.maxwell_gap = 'D:\Maxwell_sim\Zhi_peng\export\R500_gap';
    c.fld_files_variant.maxwell_gap = {'B_p1_split.fld','B_p2_split.fld','B_p3_split.fld', ...
                                       'B_p4_split.fld','B_p5_split.fld','B_p6_split.fld'};
    c.fld_variant_subdir = false;
    % [MODIFIED 2026-08-23 使用者拍板] 'scattered' -> 'grid'：Maxwell 的 .fld 是規則格，
    %   直接定位格子做三線性內插，與工作空間（conv_design_ws）用同一套實作。
    %   實測快 11 倍（52.5 s -> 4.9 s）、V 差異 max 0.31% / median 0.068%（分片線性的
    %   兩種切法之差，同階）。⚠ 既有 voltage .mat 用 'scattered' 算的，無法逐位重現。
    c.v_method          = 'grid';

    % ===== ② 幾何（**名目值**；長度單位 [m]）=====
    %   使用者拍板 2026-08-18：一律用名目值，實測值另存 *_cad 欄位供追溯。
    %   （比照 long2016：R_norm=500e-6、濾鐵錐體用名目 3.0/15.0 —— 名目是為了**結果可比性**，
    %     不是因為它更準。此處名目與實測差 ≤ 0.25 µm（0.05%），遠小於 long2016 那邊的 4.4%。）
    %
    %   名目值鏈自我檢核（全部對到 0.05 µm 內，證明這組名目值就是設計意圖）：
    %     板厚   2*R_norm/sqrt3            = 577.3503 µm   實測 上極 577.30 / 下極 580.00(捨入)
    %     導柱頂 板厚 + 4.63mm             =   5.20735 mm   實測 5.2073
    %     軛頂   導柱頂 + 3.2mm            =   8.40735 mm   實測 8.4073

    % --- magic-angle 六極尖端 ---
    c.R_norm    = 500e-6;                            % 工作半徑（名目；實測 499.905 µm）
    c.R_norm_xy = c.R_norm * sqrt(2/3);              % 408.2483 µm（實測 下 408.000 / 上 408.300）
    c.R_norm_z  = c.R_norm / sqrt(3);                % 288.6751 µm（實測 288.650）
    c.alpha     = atan2(c.R_norm_xy, c.R_norm_z);    % 54.7356°（公式鎖定，不可改；實測 54.7315°）
    c.R_norm_cad = 499.905e-6;                       % 實測（追溯用，不參與計算）

    % --- 6 極方位 / 層 / 命名（paper P1..P6）---
    c.pole_angles   = [0, 180, 120, 300, 60, 240];   % deg
    c.pole_labels   = {'P1','P2','P3','P4','P5','P6'};
    c.pole_is_lower = [1, 0, 1, 0, 0, 1];            % 1=下(P1/P3/P6，尖在板下緣)、0=上(P2/P4/P5，尖在板上緣)

    % ⭐ R500 的方位角實測 = 名目值（誤差 ≤ 0.003°）—— 比 R594 好得多（那版上極三根偏 +0.21°）。
    %    故 R_act 的正交性天然成立，無須「用名目值遷就正交約束」。
    c.pole_angles_cad = [0.0000, 180.0000, 120.0000, 300.0000, 60.0000, 240.0000];

    % --- 6 尖端位置（工作區框，原點 = 極板中面中心）---
    c.pole_tip_x    = c.R_norm_xy * cosd(c.pole_angles);
    c.pole_tip_y    = c.R_norm_xy * sind(c.pole_angles);
    c.pole_tip_z_wp = [-1, +1, -1, +1, +1, -1] * c.R_norm_z;   % 下極在板下緣、上極在板上緣

    % --- 磁極板本體（平板；濾鐵包絡用）---
    c.POLE_TIP_R      = 40e-6;                       % 尖端倒圓半徑（實測 R=0.0400 垂直軸圓柱面）
    c.POLE_TIP_BAND   = 178e-6;                      % ★ 磁極（舌片）厚度 —— 實測 178.0002 µm，
                                                     %   **六極極差 0.00000 µm**（使用者關切的一致性 ✓）
    c.POLE_WEDGE_HALF = 7.5;                         % 舌片楔形半角 [deg]（實測 7.5000，六極一致；
                                                     %   兩側壁平面與尖端 40 µm 圓角相切）
    c.POLE_WEDGE_R    = [c.R_norm_xy, 15.0e-3];      % 舌片徑向範圍：極尖 → 15 mm（接上矩形塊）
    c.POLE_BLOCK_R    = [15.0e-3, 25.0e-3];          % 外段矩形塊徑向範圍
    c.POLE_BLOCK_W    = 6.0e-3;                      % 外段矩形塊寬（橫向 ±3 mm）
    c.POLE_PLATE_T    = 2 * c.R_norm_z;              % 極板總厚 577.3503 µm（名目 = 2*R_norm_z）
    c.POLE_PLATE_T_cad = [580.00e-6, 577.30e-6];     % 實測 [下極, 上極]：下極三片被捨入成 0.58，
                                                     %   差 2.7 µm。**極尖位置不受影響**（下極尖在板底 z=0）。

    % --- ★ 無鐵安全半徑（平板極專屬；R594 版同樣現象）---
    %   極尖是**垂直圓角邊**（沿 z 長 178 µm），最近鐵是那條邊「離中面較近的一端」，
    %   不是極尖本身：位於 (r_xy, z_wp) = (408.25, ∓110.68) µm → 距中心 sqrt(408.25²+110.68²) = 423.0 µm。
    %   ⇒ R ≤ 423 µm 的取樣球內**完全沒有鐵**（R_norm=500 µm 的球殼會吃到約 2% 的鐵）。
    %   校正用 R ≤ 150 µm 離最近鐵還有 2.8 倍距離 → 濾鐵未實作對本半徑無影響。
    c.IRON_SAFE_R   = 423.0e-6;

    % --- 工作區中心在 CAD 座標中的位置（CAD z=0 = 極板下表面 = 下極尖）---
    c.WP_Z_CAD      = c.R_norm_z;                    % 288.6751 µm（= 中面；實測 288.650）
    %   build_actuator_data 用 `zwp = raw.z - cfg.SPH_OFST` 把工作區中心移到原點。
    %   ⚠ 號誌與 long2016 相反：那邊 SPH_OFST = -12.711 mm（CAD 原點在軛、工作區在下方），
    %     本 model CAD z=0 在極板下表面、工作區中面在 +288.7 µm。
    c.SPH_OFST      = c.WP_Z_CAD;

    % --- 極軸（單位向量 tip→base，工作區框）---
    %   平板六極：舌片沿**水平徑向**指向中心，上下極皆然（無 long2016 的上極傾角）。
    c.pole_axis = zeros(3, 6);
    for i = 1:6
        th = c.pole_angles(i) * pi/180;
        c.pole_axis(:,i) = [cos(th); sin(th); 0];
    end

    % --- 導柱（steel post，6 根）---
    c.POST_D        = 3.000e-3;                      % 直徑（實測 3.0000）
    c.POST_H        = 4.630e-3;                      % 高（實測 4.6300 / 4.6273，後者是板厚捨入的連帶）
    c.POST_Z        = [c.POLE_PLATE_T, c.POLE_PLATE_T + c.POST_H];   % [下, 上] = [577.35 µm, 5.20735 mm]
    c.post_xy       = [ -21.500, -21.500, -21.500,  21.500,  21.500,  21.500;
                        -19.500,   0.000,  19.500, -19.500,   0.000,  19.500] * 1e-3;
    c.post_xy_cad   = 21.502e-3;                     % 實測 |x| = 21.502（名目 21.5，差 2 µm）
    %   ⚠ 對齊備忘：AEDT R500 專案的線圈相對座標系用 x = ±21.5（與名目一致）。

    % --- 上軛（方框）---
    c.YOKE_OUT      = 53.0e-3;                       % 外框 53 x 53（外緣 ±26.5，四角 R3 @ ±23.5）
    c.YOKE_WIN      = 33.0e-3;                       % 內窗 33 x 33（內緣 ±16.5，四角 R3 @ ±13.5）
    c.YOKE_T        = 3.200e-3;                      % 厚
    c.YOKE_FILLET   = 3.000e-3;                      % 內外角圓角半徑
    c.YOKE_Z        = [c.POST_Z(2), c.POST_Z(2) + c.YOKE_T];   % [5.20735, 8.40735] mm

    % --- 整體外形 ---
    c.BBOX_XY       = 53.0e-3;                       % ±26.5 mm
    c.BBOX_Z        = [0, c.YOKE_Z(2)];              % [0, 8.40735] mm

    % ===== ③ actuator frame（供 build_actuator_data 消費）=====
    tip   = [c.pole_tip_x; c.pole_tip_y; c.pole_tip_z_wp];
    dhat  = tip ./ vecnorm(tip);
    c.R_act   = [dhat(:,1), dhat(:,3), dhat(:,5)].';           % P1,P3,P5 單位向量當列
    assert(abs(det(c.R_act) - 1) < 1e-9, 'R_act 非正交旋轉：極軸選錯或 magic angle 沒守');
    c.Pc_base = [ 1 -1  0  0  0  0;                            % 整數 canonical
                  0  0  1 -1  0  0;
                  0  0  0  0  1 -1];

    % ===== ④ 號誌慣例 + 物理常數 =====
    c.s_source = [+1, +1, +1, +1, +1, +1];           % Maxwell 匯出場已 all-source → 不翻號
    c.S_hall   = 130;                                % Hall 靈敏度 [mV/mT]（EQ-730L）
    c.N_c      = 70;                                 % 每極匝數（AEDT 激發 70 A、stranded）
    c.mu_0     = 4*pi*1e-7;
    c.k_m      = 1e-7;

    % ===== TODO =====
    %   (a) ✅ .fld 已匯出並驗證（見 fld_dir 註）→ **current 路可跑**。
    %   (b) **sensor 貼附幾何未定**（voltage 路）：long2016 的規則是「SOFF 沿錐面母線自紅點量、
    %       氣隙 0.41 mm 沿面法向」（見 utils/pole_sensor_geometry.m）。本 model 是平板 →
    %       須先定「Hall 貼在哪個面」（舌片 7.5° 側壁？極板上/下大平面？）才能套同一套公式。
    %       未定前 build_V_matrix 不可用。
    %   (c) 濾鐵：filter_iron_nodes 偵測到缺 POLE_R/POLE_CONE_LEN 會警告並跳過錐體判定。
    %       R ≤ 423 µm（IRON_SAFE_R）內本來就沒有鐵，故 R=150 µm 的校正不受影響；
    %       若日後要取 R > 423 µm，須另寫平板包絡（舌片：楔角 7.5° + 178 µm 厚 + 40 µm 尖端圓角）。
    %   (d) coil→pole 對應：AEDT 線圈以矩形排列（±21.5, {0, ±19.5}），與極的 60° 方位不是同一套編號
    %       → 首次載入後務必用 K̄_I 檢驗（對角占優 + 對角全正），per pole-coil-numbering 規則。
    %       ⚠ 但 R594 實測顯示本設計**六極不等強**（±x 兩根強 2.2 倍），故「非對角全負」那條
    %         恆不成立，不可據此判定擬合失敗（見 conv_design 的 ki_req 旗標）。
end
