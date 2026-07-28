function c = mt_constants()
%MT_CONSTANTS  NTU_hexapole 的 per-model 設定（供 function/model_config.m 呼叫）。
%   扁平淚滴板、兩平行層、6 tip 匯於 WP（與 long2016 錐尖-球面-magic angle 不同）。
%   物理常數 = NTU 自己的（原 current_base/code/main_function/mt_constants.m 忠實搬移）。

    % ===== ① 提取識別 / 資料路由 =====
    c.strategy          = 'ntu_flat';                % 扁平板：不旋轉、不濾鐵、WP 形心置中 + R_load 預篩
    c.apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];
    c.default_variant   = 'full_assembly';
    c.regions           = {'all'};                   % NTU 只匯出 _all
    c.R_load            = 2e-3;                       % WP 預篩球半徑 [m]

    % ===== ② 幾何 / 物理常數 =====
    c.pole_labels   = {'P1','P2','P3','P4','P5','P6'};
    c.pole_angles   = [0, 180, 120, 300, 60, 240];   % deg, P1..P6
    c.pole_is_lower = [1, 0, 1, 0, 0, 1];

    % 兩層平板 z（板底 z_b）[m]
    ZB_LO = 8.80e-3;  ZB_UP = 9.55e-3;  TH = 0.25e-3;
    c.plate_zb = ZB_LO*c.pole_is_lower + ZB_UP*(~c.pole_is_lower);
    c.plate_TH = TH;

    % 6 tip（world）[m]：tip = TIP_R·(cosθ,sinθ)，z=面向 WP 的板面
    c.TIP_R = 0.353e-3;
    th = c.pole_angles * pi/180;
    c.pole_tip_x = c.TIP_R * cos(th);
    c.pole_tip_y = c.TIP_R * sin(th);
    c.pole_tip_z = c.plate_zb + c.plate_TH * c.pole_is_lower;   % 下 +0.25→9.05；上 +0→9.55

    % WP = 6 tip 形心 [m] ≈ (0,0,9.30mm)
    c.WP = [mean(c.pole_tip_x); mean(c.pole_tip_y); mean(c.pole_tip_z)];

    % tip 於 WP 框（定義 ntu_flat 的 Pc_base）
    c.pole_tip_x_wp = c.pole_tip_x - c.WP(1);
    c.pole_tip_y_wp = c.pole_tip_y - c.WP(2);
    c.pole_tip_z_wp = c.pole_tip_z - c.WP(3);
    c.R_norm = mean(vecnorm([c.pole_tip_x_wp; c.pole_tip_y_wp; c.pole_tip_z_wp]));   % ~0.515mm

    % 6 S2 sensor（板上方）+ Hall
    SX = 13.15e-3;  SY = -0.99e-3;  SZOFF = 0.66e-3;
    c.sensor_pos = zeros(3,6);  c.sensor_n = zeros(3,6);
    for i = 1:6
        t = th(i);
        c.sensor_pos(:,i) = [ SX*cos(t) - SY*sin(t);
                              SX*sin(t) + SY*cos(t);
                              c.plate_zb(i) + SZOFF ];
        c.sensor_n(:,i)   = [0; 0; 1];               % S2：n+ = +z
    end
    c.SENSOR_R = 0.15e-3;  c.SENSOR_H = 0.10e-3;  c.S_hall = 130;

    % 號誌慣例 + 物理常數
    c.s_source = [+1, +1, +1, +1, +1, +1];           % NTU raw 六極全 source（不需翻）
    c.N_c      = 50;
    c.mu_0     = 4*pi*1e-7;
    c.k_m      = 1e-7;
end
