function c = mt_constants()
%MT_CONSTANTS  hung_hexapole 的 per-model 設定（供 function/model_config.m 呼叫）。
%   物理常數 = hung **自己的**（backup/hung/analysis/core/mt_constants.m 忠實搬移）：
%   POLE_R=3.175mm / CONE_LEN=15.875mm / SPH_OFST=0（WP 在原點）。
%   ⚠ 這修正了 hung 舊 pipeline（誤 addpath long2016 analysis → 用 long2016 常數 3.0/15.0/−12.711）。
%   → hung 校正結果會變、**不再對舊 loader 等價**（舊的用錯常數）。

    % ===== ① 提取識別 / 資料路由 =====
    c.strategy          = 'hex_magic';
    c.apdl_to_paper_idx = [1, 2, 3, 4, 5, 6];        % hung = identity（deck 照 paper 序）
    c.default_variant   = 'R500';                    % [MODIFIED] R500 = 原始設計（R_sphere=0.5mm），與下方 R_norm=500e-6 對齊
                                                     %   R300/R500/R700 六 coil mesh 各自一致；no_gap/gap_200um 不一致→crash
    c.regions           = {'all', 'wp'};
    c.R_load            = [];

    % ===== ② 幾何 / 物理常數（hung 自己的）=====
    c.R_norm    = 500e-6;
    c.R_norm_xy = c.R_norm * sqrt(2/3);
    c.R_norm_z  = c.R_norm / sqrt(3);
    c.alpha     = atan2(c.R_norm_xy, c.R_norm_z);
    c.SPH_OFST  = 0;                                 % hung WP = 原點（R700/R300 資料 WP-centered）

    c.pole_angles   = [0, 180, 120, 300, 60, 240];   % deg
    c.pole_labels   = {'P1','P2','P3','P4','P5','P6'};
    c.pole_is_lower = [1, 0, 1, 0, 0, 1];

    % 6 尖端（WP 框）= 精確 magic-angle 格「現算」（R_norm_xy/z 精確公式；
    %   ⚠ hung 原 mt_constants 用 magic_angle=54.74 圓整值 → tip 偏離正交格 ~1e-4、R_act assert 過不了。
    %   精確磁角 = atan(√2) 是普世鎖死幾何，hung 歷史上（借 long2016 常數時）本就用精確 tip；此處沿用。）
    c.pole_tip_x    = c.R_norm_xy * cosd(c.pole_angles);
    c.pole_tip_y    = c.R_norm_xy * sind(c.pole_angles);
    c.pole_tip_z_wp = [-1, +1, -1, +1, +1, -1] * c.R_norm_z;   % SPH_OFST=0 → z_wp = z

    % 濾鐵用錐體幾何（hung 自己尺寸）
    c.POLE_TIP_R    = 40e-6;
    c.POLE_R        = 3.175e-3;
    c.POLE_CONE_LEN = 15.875e-3;

    % 極軸（hung 慣例：單位徑向穿過 tip）
    c.pole_axis = [c.pole_tip_x; c.pole_tip_y; c.pole_tip_z_wp] / c.R_norm;

    % actuator frame（供 build_actuator_data「消費」；magic-angle 尖端 → canonical 電荷格，與其原式一致）
    tip   = [c.pole_tip_x; c.pole_tip_y; c.pole_tip_z_wp];
    dhat  = tip ./ vecnorm(tip);
    c.R_act   = [dhat(:,1), dhat(:,3), dhat(:,5)].';
    c.Pc_base = [ 1 -1  0  0  0  0;                   % 整數 canonical（baseline 逐位元不變）
                  0  0  1 -1  0  0;
                  0  0  0  0  1 -1];

    % 號誌慣例 + 物理常數
    c.s_source = [-1, -1, -1, -1, -1, -1];           % hung raw 六極全 sink → 翻全部成 source
    c.N_c      = 70;
    c.mu_0     = 4*pi*1e-7;
    c.k_m      = 1e-7;
end
