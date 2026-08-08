function c = mt_constants()
%MT_CONSTANTS  long2016_hexapole_halfcut（**Maxwell** 求解器）的 per-model 設定。
%   結構與 APDL 版 config/long2016_hexapole_halfcut/tip40um/mt_constants.m 相同，
%   幾何 / 物理常數逐行等價（同一顆物理模型），**只有三處為 Maxwell 專屬**：
%     ① c.s_source = 全 +1 —— Maxwell 匯出的場已是 all-source，**不再翻號**（使用者拍板 2026-07-31）。
%     ② c.apdl_to_paper_idx = identity —— .fld 檔名 B_p1..B_p6 直接就是 paper 極 P1..P6。
%     ③ 新增 c.fld_dir / c.fld_files —— 供 extract_maxwell_data 定位 .fld。
%   ⚠ 座標系：.fld 為 Maxwell/STEP frame（= APDL raw frame，WP 在 z ≈ −12.711 mm = SPH_OFST），
%     故 build_actuator_data 的 zoff = −SPH_OFST 照樣成立。

    % ===== ① 提取識別 / 資料路由 =====
    c.strategy          = 'hex_magic';               % magic-angle 六極：R_act 旋轉 + filter_iron
    c.apdl_to_paper_idx = [1, 2, 3, 4, 5, 6];        % Maxwell：檔名 B_p<k> = paper P<k>（identity）
    c.default_variant   = 'maxwell';
    c.regions           = {'all'};                    % Maxwell 一個 .fld = 一個場（無 all/wp/circuit 之分）
    c.R_load            = [];                         % hex 不用 WP 預篩球

    % Maxwell .fld 路由（extract_maxwell_data 用）；**兩組匯出、用途不同**：
    %   dataset='all'（預設）→ WP 細格，供電荷擬合（R≤150µm 球內 ~1767 格點）
    %   dataset='voltage'    → sensor 粗格，供 build_V_matrix（sensor 在 WP 外 4.5mm，細格框涵蓋不到）
    c.fld_dir           = 'D:\Maxwell_sim\long2016_hexapole_halfcut\export';
    c.fld_files         = {'B_p1.fld','B_p2.fld','B_p3.fld','B_p4.fld','B_p5.fld','B_p6.fld'};
    %   ↑ 步距 0.02mm、框 x,y∈±0.6mm z∈[-13.31,-12.11]mm、61³=226,981 格點
    % [ADDED 2026-08-05] variant → WP 細格 .fld 檔名（dataset='all' 用；voltage 一律走 fld_files_voltage）。
    %   兩組**匯出格點完全相同**（header 逐字一致），差別只在 Maxwell 側 Sphere1mm 的網格尺寸：
    %     maxwell          = Sphere1mm 0.1 mm（定案基準，134,067 tets）
    %     maxwell_mesh0p06 = Sphere1mm 0.06 mm（網格收斂測試；同框同步距 → 可與上者逐點相減）
    c.fld_files_variant.maxwell = c.fld_files;
    c.fld_files_variant.maxwell_mesh0p06 = ...
        {'B_p1_0.06.fld','B_p2_0.06.fld','B_p3_0.06.fld', ...
         'B_p4_0.06.fld','B_p5_0.06.fld','B_p6_0.06.fld'};
    c.fld_files_voltage = {'B_voltage_p1.fld','B_voltage_p2.fld','B_voltage_p3.fld', ...
                           'B_voltage_p4.fld','B_voltage_p5.fld','B_voltage_p6.fld'};
    %   ↑ 步距 0.1mm、框 x∈±15.23 y∈±14 z∈[-16.1,0]mm、**305×281×162 = 13,884,210** 格點（含兩層 sensor）
    %   ⚠ x 只有 305 個（不是 306）：跨距 30.46mm 不是 0.1 的整數倍（30.46/0.1 = 304.6）→ AEDT 生到
    %     x = +15.17 就停，**+x 端被截掉 0.06mm、格點左右不對稱**（y/z 整除故對稱）。P1 是 +x 下極、
    %     錐底約 x≈15.24mm 落在截掉處之外 → 畫到極根部（s > 14.76）時會沒有資料。
    c.fld_variant_subdir = false;                     % .fld 直接放 fld_dir，不再進 variant 子夾

    % Maxwell 無 ANSYS 的 sensor_local tet CSV → V 矩陣一律走座標式 scatteredInterpolant
    c.v_method          = 'scattered';
    c.sensor_r_loc      = 1.5e-3;                     % scattered 取源鄰域半徑 [m]（0.1mm 格 → 約 1.4 萬點/sensor）

    % ===== ② 幾何 / 物理常數（與 APDL 版逐行等價）=====
    % 工作半徑（magic-angle 鎖死公式）
    c.R_norm    = 500e-6;
    c.R_norm_xy = c.R_norm * sqrt(2/3);              % ~408 µm
    c.R_norm_z  = c.R_norm / sqrt(3);                % ~289 µm
    c.alpha     = atan2(c.R_norm_xy, c.R_norm_z);    % 尖端極角 ~54.74°

    % Yoke / protrusion（供 SPH_OFST 與 upper_incline 導出）
    c.PROT_H     = 7.0e-3;
    c.YOKE_IN_R  = 84e-3/2;
    c.YOKE_OUT_R = 106e-3/2;
    c.YOKE_MID_R = (c.YOKE_IN_R + c.YOKE_OUT_R)/2;   % ~47.5 mm
    c.SPH_OFST   = -c.PROT_H - 6e-3 + c.R_norm_z;    % WP 相對原點 ~ -12.711 mm

    % 6 極方位 / 層 / 命名（paper P1..P6）
    c.pole_angles   = [0, 180, 120, 300, 60, 240];   % deg
    c.pole_labels   = {'P1','P2','P3','P4','P5','P6'};
    c.pole_is_lower = [1, 0, 1, 0, 0, 1];            % 1=下(P1/P3/P6)、0=上(P2/P4/P5)

    % 6 尖端位置（WP 框）= 由 R_norm + 方位 + 上下號誌「現算」
    c.pole_tip_x    = c.R_norm_xy * cosd(c.pole_angles);
    c.pole_tip_y    = c.R_norm_xy * sind(c.pole_angles);
    c.pole_tip_z_wp = [-1, +1, -1, +1, +1, -1] * c.R_norm_z;

    % 濾鐵用錐體幾何（**名目值**，勿改：filter_iron_nodes 用它篩鐵，動了會改變電荷擬合）
    c.POLE_TIP_R    = 40e-6;                          % 尖端倒圓半徑
    c.POLE_R        = 3e-3;                           % 底半徑（6mm 直徑之半）
    c.POLE_CONE_LEN = 15e-3;                          % 錐長

    % [ADDED 2026-08-08] sensor 貼附面用的**真實** per-layer 錐體斜率 tan(beta)（CAD STEP 實測）
    %   下極 R(s) = 0.0327 + 0.20330*s  (beta = 11.4916 deg)
    %   上極 R(s) = 0.0330 + 0.19463*s  (beta = 11.0138 deg)   [mm]，s 自極尖點沿極軸量
    %   名目值 3.0/15.0 (11.310 deg) 對兩層都不對 -> 只在 build_V_matrix 的 sensor 幾何改用真值；
    %   截距 R0 不寫死，由 POLE_TIP_R 經虛擬錐頂公式導出（見 sensor_geometry），故 tip 變體自動跟著走。
    c.pole_cone_slope = [0.20330, 0.19463];           % [下極, 上極]

    % [ADDED 2026-08-08] 尖端倒圓的**軸向推進量**（CAD STEP 實測，六根極一致）。
    %   = 錐面與 40µm 倒圓球的「切點」距極尖的軸向距離。切點落在球面 90deg-beta 處（非赤道），
    %   故軸向只推進 0.0322mm、**不是 0.04**（實測交界圓 t=0.03216 / r=0.0391~0.0394，
    %   反推 r_f = r/cos(beta) = 0.0398~0.0402 => 倒圓確實是 40µm）。理論值 r_f*(1-sin beta)
    %   = 0.03203(下)/0.03236(上)，與實測差 0.1~0.2µm。缺此欄位時 sensor_geometry 用該理論式回退。
    c.pole_tip_axial = 0.0322e-3;                     % [m]

    % [ADDED 2026-08-08] 真實錐長（自極尖點沿極軸；CAD STEP 實測）。供畫磁極輪廓用。
    %   底半徑可由 R(cone_len) = R0 + slope*cone_len 導出（下 3.0470 / 上 3.0000 mm）。
    c.pole_cone_len_real = [14.8267e-3, 15.2443e-3];  % [下極, 上極]

    % 上極傾角（幾何導出）→ 供 pole_axis
    YOKE_H      = 2e-3;
    end_upper_r = c.YOKE_MID_R - 11.5e-3;             % 36 mm
    end_upper_z = YOKE_H + c.PROT_H + 5e-3;           % 14 mm
    tip_upper_z = -c.PROT_H - 6e-3 + 2*c.R_norm_z;
    c.upper_incline = atan2(end_upper_z - tip_upper_z, end_upper_r - c.R_norm_xy);  % ~36.6°

    % 極軸（單位向量 tip→base，WP 框）：下極水平徑向、上極沿傾角
    c.pole_axis = zeros(3, 6);
    for i = 1:6
        th = c.pole_angles(i) * pi/180;
        if c.pole_is_lower(i)
            c.pole_axis(:,i) = [cos(th); sin(th); 0];
        else
            inc = c.upper_incline;
            c.pole_axis(:,i) = [cos(inc)*cos(th); cos(inc)*sin(th); sin(inc)];
        end
    end

    % actuator frame（供 build_actuator_data「消費」）
    tip   = [c.pole_tip_x; c.pole_tip_y; c.pole_tip_z_wp];
    dhat  = tip ./ vecnorm(tip);
    c.R_act   = [dhat(:,1), dhat(:,3), dhat(:,5)].';
    c.Pc_base = [ 1 -1  0  0  0  0;                   % 整數 canonical（避免浮點漂移）
                  0  0  1 -1  0  0;
                  0  0  0  0  1 -1];

    % 號誌慣例 + 物理常數
    c.s_source = [+1, +1, +1, +1, +1, +1];           % Maxwell 場已 all-source → 不翻號
    c.S_hall   = 130;                                % Hall 靈敏度 [mV/mT]（EQ-730L；voltage 路用）
    c.N_c      = 70;                                  % 每極匝數
    c.mu_0     = 4*pi*1e-7;
    c.k_m      = 1e-7;                                % mu_0/(4π)
end
