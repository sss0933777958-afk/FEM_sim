function c = mt_constants()
% MT_CONSTANTS  NTU_hexapole 幾何常數（扁平淚滴板、兩平行層；供 calibration 管線用）。
%   來源：apdl/NTU_hexapole/geom/export/{full_assembly_params.md,pole_params.md} + 使用者實測 tip。
%   ⚠ 與 long2016（錐尖-球面-magic angle）不同：NTU 是扁平板，6 tip 匯於 WP。
%   單位一律 SI（m、T→內部；場對外 mT）。
%
%   6 極（paper P1..P6）方位/層（apdl coil→paper = [1,3,6,5,2,4]）：
%     P1 0°下、P2 180°上、P3 120°下、P4 300°上、P5 60°上、P6 240°下。

    % ---- pole 命名 / 對應 ----
    c.pole_labels   = {'P1','P2','P3','P4','P5','P6'};
    c.apdl_to_paper = {'P1','P3','P6','P5','P2','P4'};      % coil1→P1, coil2→P3, ...
    c.pole_angles   = [0, 180, 120, 300, 60, 240];          % deg, indexed P1..P6
    c.pole_is_lower = [1, 0, 1, 0, 0, 1];                   % P1,P3,P6 下層; P2,P4,P5 上層

    % ---- 兩層平板 z（板底 z_b）[m] ----
    ZB_LO = 8.80e-3;  ZB_UP = 9.55e-3;   TH = 0.25e-3;      % 板厚 0.25mm
    c.plate_zb = ZB_LO*c.pole_is_lower + ZB_UP*(~c.pole_is_lower);   % P1..P6 板底 z
    c.plate_TH = TH;

    % ---- 6 tip（world）[m]：tip = 0.353·(cosθ,sinθ)，z=「面向 WP 的板面」----
    %   使用者拍板（配對例）：下極 tip (0.353,0,9.05)[下板頂]、上極 tip (−0.353,0,9.55)[上板底]，
    %   對徑配對連線中點 → WP=(0,0,9.30)。故 下極 tip z = 板底+板厚(9.05)、上極 tip z = 板底(9.55)。
    c.TIP_R = 0.353e-3;                                     % tip 到 z 軸半徑
    th = c.pole_angles * pi/180;
    c.pole_tip_x = c.TIP_R * cos(th);                       % 1x6 [m]
    c.pole_tip_y = c.TIP_R * sin(th);
    c.pole_tip_z = c.plate_zb + c.plate_TH * c.pole_is_lower;   % 下極 +0.25→9.05；上極 +0→9.55

    % ---- WP（工作點）= 6 tip 形心 [m] = (0,0,9.30mm)（對徑配對連線中點）----
    c.WP = [mean(c.pole_tip_x); mean(c.pole_tip_y); mean(c.pole_tip_z)];

    % ---- tip 於 WP-centered 框 [m]（下游 charge model 用）----
    c.pole_tip_x_wp = c.pole_tip_x - c.WP(1);
    c.pole_tip_y_wp = c.pole_tip_y - c.WP(2);
    c.pole_tip_z_wp = c.pole_tip_z - c.WP(3);               % 下 −0.375 / 上 +0.375 mm
    tipwp = [c.pole_tip_x_wp; c.pole_tip_y_wp; c.pole_tip_z_wp];  % 3x6
    c.R_norm = mean(vecnorm(tipwp));                        % ~0.515mm（6 tip 等距 WP）

    % ---- 6 S2 sensor（板上方）：pole-local (13.15,−0.99)、world z=z_b+0.66、n+=+z ----
    %   0.66 = 板厚 0.25 + 氣隙 0.41（板頂上方 0.41mm）。比照 field_viz/sensor_average + plot_poles_sensors_3d。
    SX = 13.15e-3;  SY = -0.99e-3;  SZOFF = 0.66e-3;
    c.sensor_pos = zeros(3,6);  c.sensor_n = zeros(3,6);
    for i = 1:6
        t = th(i);
        c.sensor_pos(:,i) = [ SX*cos(t) - SY*sin(t);
                              SX*sin(t) + SY*cos(t);
                              c.plate_zb(i) + SZOFF ];
        c.sensor_n(:,i)   = [0; 0; 1];                      % S2：n+ = +z（離板朝上）
    end
    c.SENSOR_R = 0.15e-3;   c.SENSOR_H = 0.10e-3;           % 取樣圓柱 半徑/高
    c.S_hall   = 130;                                       % Hall 靈敏度 [mV/mT]

    % ---- 物理常數 ----
    c.N_c  = 50;                                            % 匝數（NTU：1A×50 turns）
    c.mu_0 = 4*pi*1e-7;
    c.k_m  = 1e-7;                                          % mu_0/(4π)
end
