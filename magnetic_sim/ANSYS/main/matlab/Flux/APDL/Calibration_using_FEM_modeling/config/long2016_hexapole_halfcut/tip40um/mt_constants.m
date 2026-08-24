function c = mt_constants()
%MT_CONSTANTS  long2016_hexapole_halfcut 的 per-model 設定（供 function/model_config.m 呼叫）。
%   回一個 struct c，含兩塊：
%     ① 提取識別 / 資料路由（strategy / map / 變體 / regions / 預篩）
%     ② 幾何 / 物理常數（primitives；6 極 tip、pole_axis 由 primitives「現算」，不存 18 數）
%   ③ per-run 調參（R_select / ell0 / USE_BIAS…）不在此，留 main.m。
%   物理常數值 = backup/hexapole-long2016/analysis/mt_constants.m 的忠實搬移（保持等價）。

    % ===== ① 提取識別 / 資料路由 =====
    c.strategy          = 'hex_magic';               % magic-angle 六極：R_act 旋轉 + filter_iron
    c.apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];        % coil k 激發此 paper 極（per-model，不可互抄）
    c.default_variant   = 'graded';
    c.regions           = {'all', 'wp'};
    c.R_load            = [];                         % hex 不用 WP 預篩球

    % ===== ② 幾何 / 物理常數 =====
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
    c.SPH_OFST   = -c.PROT_H - 6e-3 + c.R_norm_z;    % WP 相對 APDL 原點 ~ -12.711 mm

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

    % 上極傾角（APDL 幾何導出）→ 供 pole_axis
    YOKE_H      = 2e-3;
    end_upper_r = c.YOKE_MID_R - 11.5e-3;             % 36 mm
    end_upper_z = YOKE_H + c.PROT_H + 5e-3;           % 14 mm (APDL 座標)
    tip_upper_z = -c.PROT_H - 6e-3 + 2*c.R_norm_z;    % APDL 座標
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

    % actuator frame（供 build_actuator_data「消費」；與其原式一致：magic-angle 尖端 → canonical 電荷格）
    tip   = [c.pole_tip_x; c.pole_tip_y; c.pole_tip_z_wp];
    dhat  = tip ./ vecnorm(tip);
    c.R_act   = [dhat(:,1), dhat(:,3), dhat(:,5)].';
    c.Pc_base = [ 1 -1  0  0  0  0;                   % 整數 canonical（避免浮點漂移；baseline 逐位元不變）
                  0  0  1 -1  0  0;
                  0  0  0  0  1 -1];

    % 號誌慣例 + 物理常數
    c.s_source = [-1, +1, -1, +1, +1, -1];           % raw FEM→all-source（翻下極 sink P1/P3/P6）
    c.S_hall   = 130;                                % Hall 靈敏度 [mV/mT]（EQ-730L；voltage 路用）
    c.N_c      = 70;                                  % 每極匝數
    c.mu_0     = 4*pi*1e-7;
    c.k_m      = 1e-7;                                % mu_0/(4π)
end
