function c = mt_constants()
%MT_CONSTANTS  zhi_peng 平面六極（V2_split）的 per-model 設定（**Maxwell** 求解器）。
%   幾何全部由 CAD 實測，來源檔：
%     CAD_model/zhi-Peng/STEP/zhi_peng_hexapole_V2_split.STEP
%   （OCC/OCP 解析面量測；該 STEP 為純解析幾何 PLANE/CYLINDRICAL_SURFACE/CIRCLE，無 B-spline，
%     故尺寸可直接讀出，不受「控制點 ≠ 曲面邊界」的陷阱影響。）
%
%   結構比照 config/long2016_hexapole_halfcut/tip40um/mt_constants.m。
%
%   ⚠ 與 long2016 的差別（別互抄）：
%     ① 磁極是**平面板**（厚 0.686 mm、尖端 40µm 倒圓），不是圓錐 → 無 POLE_R / POLE_CONE_LEN，
%        濾鐵幾何需另寫（見檔尾 TODO）。
%     ② 六極**共平面製作**：下極尖在板下緣 z=0、上極尖在板上緣 z=0.686，工作區中心在中面 z=0.343。
%        雖然做法不同，magic-angle 條件仍成立（見下）。

    % ===== ① 提取識別 / 資料路由 =====
    c.strategy          = 'hex_magic';               % magic-angle 六極：R_act 旋轉 + 濾鐵
    c.apdl_to_paper_idx = [1, 2, 3, 4, 5, 6];        % identity：.fld 檔名 B_p<k> = paper P<k>
                                                     % （新 model 一律 identity，per pole-coil-numbering 規則）
    c.default_variant   = 'maxwell';
    c.regions           = {'all'};
    c.R_load            = [];

    % [MODIFIED 2026-08-18] .fld 已搬進 R594/ 子夾（因為新增了 R500 變體，兩組並存）。
    %   原路徑 'D:\Maxwell_sim\Zhi_peng\export' 現在只有 R500/ 與 R594/ 兩個子夾、沒有 .fld。
    c.fld_dir           = 'D:\Maxwell_sim\Zhi_peng\export\R594';
    c.fld_files         = {'B_p1.fld','B_p2.fld','B_p3.fld','B_p4.fld','B_p5.fld','B_p6.fld'};
    %   [MODIFIED 2026-08-15] **六個 .fld 都已匯出**，且框已由薄片改成立方體（原註解已過時）：
    %   實測 header：Min [-2mm -2mm -1.66mm]  Max [2mm 2mm 2.34mm]  Grid [0.02mm]³
    %   → 201 x 201 x 201 = 8,120,601 格點/檔（每檔 ~1.15 GB，六檔約 6.9 GB）。
    %   z 框中心 = (-1.66+2.34)/2 = 0.34 mm，與工作區中面 WP_Z_CAD = 0.343 mm 一致。
    c.fld_variant_subdir = false;
    % [MODIFIED 2026-08-23 使用者拍板] 'scattered' -> 'grid'：Maxwell 的 .fld 是規則格，
    %   直接定位格子做三線性內插，與工作空間（conv_design_ws）用同一套實作。
    %   實測快 11 倍（52.5 s -> 4.9 s）、V 差異 max 0.31% / median 0.068%（分片線性的
    %   兩種切法之差，同階）。⚠ 既有 voltage .mat 用 'scattered' 算的，無法逐位重現。
    c.v_method          = 'grid';

    % ===== ② 幾何（CAD 實測；長度單位 [m]）=====
    % --- magic-angle 六極尖端 ---
    %   實測：極尖在平面內半徑 485.0 µm、離中面 343.0 µm（= 板厚 686 µm 的一半）。
    %   驗證 magic angle：atan(485.0/343.0) = 54.73° = arctan(sqrt2) ✓
    %                     R = sqrt(485.0^2 + 343.0^2) = 594.03 µm
    %                     R*sqrt(2/3) = 485.06 ✓ ；R/sqrt(3) = 342.96 ✓
    c.R_norm    = 594.0e-6;                          % 六極尖共球半徑（CAD 實測）
    c.R_norm_xy = c.R_norm * sqrt(2/3);              % 485.1 µm（實測 485.0）
    c.R_norm_z  = c.R_norm / sqrt(3);                % 343.0 µm（實測 343.0 = 板厚/2）
    c.alpha     = atan2(c.R_norm_xy, c.R_norm_z);    % 54.74°（公式鎖定，不可改）

    % --- 6 極方位 / 層 / 命名（paper P1..P6）---
    c.pole_angles   = [0, 180, 120, 300, 60, 240];   % deg（**名目值**，見下方 pole_angles_cad）
    c.pole_labels   = {'P1','P2','P3','P4','P5','P6'};
    c.pole_is_lower = [1, 0, 1, 0, 0, 1];            % 1=下(P1/P3/P6，尖在板下緣)、0=上(P2/P4/P5，尖在板上緣)

    % ⚠ CAD 實測方位角與名目值有偏差（下極三根準確、上極三根偏 +0.21°）：
    %     P1 0.000  P3 120.000  P6 240.000   (下極，與名目一致)
    %     P2 180.210  P4 299.790  P5 60.210  (上極，偏 ±0.21° → 極尖橫向偏 ~1.8 µm)
    %   本檔的 pole_angles 用**名目值**，因為 actuator frame 需要 P1/P3/P5 兩兩正交
    %   （det(R_act)=1；hexapole 正交約束為專案硬性規定）。實測值另存供追溯 / 誤差評估。
    c.pole_angles_cad = [0.000, 180.210, 120.000, 299.790, 60.210, 240.000];

    % --- 6 尖端位置（工作區框，原點 = 極板中面中心）---
    c.pole_tip_x    = c.R_norm_xy * cosd(c.pole_angles);
    c.pole_tip_y    = c.R_norm_xy * sind(c.pole_angles);
    c.pole_tip_z_wp = [-1, +1, -1, +1, +1, -1] * c.R_norm_z;   % 下極在板下緣、上極在板上緣

    % --- 磁極板本體 ---
    c.POLE_TIP_R    = 40e-6;                         % 尖端倒圓半徑（實測 R=0.04 圓柱面）
    c.POLE_PLATE_T  = 686e-6;                        % 極板總厚（z 0 .. 0.686 mm）= 2*R_norm_z
    c.POLE_TIP_BAND = 178e-6;                        % 尖端段厚：下極 z 0~0.178、上極 z 0.508~0.686
    c.POLE_R_OUTER  = [25.2947e-3, 34.7635e-3];      % 極板外緣最遠半徑 [±x 兩根, 四根對角根]

    % --- 工作區中心在 CAD 座標中的位置（CAD z=0 為極板下表面）---
    c.WP_Z_CAD      = 343e-6;                        % 工作區中心 z（= 中面）；.fld 匯出框正是繞此
    % [ADDED 2026-08-15] build_actuator_data 用 `zwp = raw.z - cfg.SPH_OFST` 把工作區中心移到原點
    %   → 本 model 的 SPH_OFST 就是 WP_Z_CAD（CAD z=0 在極板下表面，工作區中面在 +343 µm）。
    %   ⚠ 與 long2016 號誌相反：那邊 SPH_OFST = -12.711 mm（CAD 原點在軛，工作區在下方）。
    c.SPH_OFST      = c.WP_Z_CAD;

    % --- 導柱（steel post，6 根；CAD 實測）---
    c.POST_D        = 3.000e-3;                      % 直徑
    c.POST_H        = 4.630e-3;                      % 高（z 0.686 .. 5.316 mm）
    c.POST_Z        = [686e-6, 5316e-6];             % [下, 上]
    c.post_xy       = [ -20.825, -20.825, -20.825,  20.825,  20.825,  20.825;
                         -19.500,   0.000,  19.500, -19.500,   0.000,  19.500] * 1e-3;
    %   ⚠ 對齊備忘：Maxwell AEDT 專案（D:\Maxwell_sim\Zhi_peng\project\Project1.aedt）內的線圈
    %     用 x = ±20.82，CAD 實測為 ±20.825（差 5 µm）。線圈另有 R=6.5mm/H=4.63mm 圓柱。

    % --- 上軛（方框，CAD 實測）---
    c.YOKE_OUT      = 53.0e-3;                       % 外框 53 x 53（外緣 ±26.5，四角 R3 @ ±23.5）
    c.YOKE_WIN      = 33.0e-3;                       % 內窗 33 x 33（內緣 ±16.5，四角 R3 @ ±13.5）
    c.YOKE_T        = 3.200e-3;                      % 厚（z 5.316 .. 8.516 mm）
    c.YOKE_FILLET   = 3.000e-3;                      % 內外角圓角半徑

    % --- 整體外形 ---
    c.BBOX_XY       = 53.0e-3;                       % ±26.5 mm
    c.BBOX_Z        = [0, 8.516e-3];

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
    c.N_c      = 70;                                 % 每極匝數（AEDT 激發 Coil1 = 70 A、stranded）
    c.mu_0     = 4*pi*1e-7;
    c.k_m      = 1e-7;

    % ===== TODO（尚未定，需使用者拍板；未定前 voltage 路與濾鐵不可用）=====
    %   (a) 濾鐵幾何：filter_iron_nodes 目前假設**圓錐**磁極（吃 POLE_R / POLE_CONE_LEN），
    %       本 model 是平面板 → 需改寫成板狀包絡（板厚 686µm + 尖端 40µm 倒圓 + 楔形側面）。
    %       在那之前不要沿用 long2016 的錐體參數，會把空氣點誤判成鐵、或反之。
    %   (b) sensor 貼附幾何（voltage 路）：本 model 尚未定義 Hall sensor 位置與法線
    %       → 無 sensor_r_loc / soff / face 定義，build_V_matrix 不可用。
    %   (c) p2~p6 的 .fld 尚未匯出（見 fld_files 註）。
    %   (d) coil→pole 對應：AEDT 專案的線圈以矩形排列（±20.825, {0,±19.5}），
    %       與極的 60° 方位不是同一套編號 → 首次載入後務必用 K̄_I 檢驗
    %       （對角占優 + 對角全正 + 非對角全負），per pole-coil-numbering 規則。
end
