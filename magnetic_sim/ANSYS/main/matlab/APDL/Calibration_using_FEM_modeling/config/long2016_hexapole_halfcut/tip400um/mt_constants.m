function c = mt_constants()
%MT_CONSTANTS  long2016_hexapole_halfcut / tip400um（CNC 400µm 鈍尖、正交修正 scenario-A）幾何 config。
%   相對 tip40um：尖端沿錐軸後退 → 用「尖端到尖端」R_act、Pc_base 非 canonical、sensor 幾何後退 +
%   上極 CONE_ANG=35.49° + 鈍尖 gap；voltage 用 scattered（sensor_local CSV 與 solve mesh 不符）、
%   WP fit 內插到 tip40um 點雲。build_actuator_data/build_V/main 只「消費」本檔提供的 R_act/Pc_base/sensor/旗標。
%   幾何常數 = deck MT_Mesh_Graded.txt @POLE_TIP_R=400e-6 評值、實測對齊 solved mesh 0µm（見 memory
%   project_long2016_tip400_sensor_calib）。跑法：main.m 設 GEOM='tip400um', BASE='voltage'（interp_to/v_method 自動）。

    % ===== ① 提取識別 / 資料路由 =====
    c.strategy          = 'hex_magic';
    c.apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];
    c.default_variant   = 'tip400um';
    c.regions           = {'all'};
    c.R_load            = [];

    % ===== ② primitives（量級同 tip40um → β / SPH_OFST 一致）=====
    c.R_norm    = 500e-6;
    c.R_norm_xy = c.R_norm * sqrt(2/3);
    c.R_norm_z  = c.R_norm / sqrt(3);
    c.PROT_H    = 7.0e-3;
    c.SPH_OFST  = -c.PROT_H - 6e-3 + c.R_norm_z;      % ~ -12.711 mm
    c.pole_angles   = [0, 180, 120, 300, 60, 240];
    c.pole_labels   = {'P1','P2','P3','P4','P5','P6'};
    c.pole_is_lower = [1, 0, 1, 0, 0, 1];
    c.POLE_TIP_R    = 400e-6;                          % 鈍尖倒圓
    c.POLE_R        = 3e-3;
    c.POLE_CONE_LEN = 15e-3;

    % ===== ③ 後退尖端（WP frame [m]，paper P1..P6）=====
    tips = 1e-3 * [ 1.8839, -0.9001, -0.9419,  0.4501,  0.4501, -0.9419;    % x
                    0.0000,  0.0000,  1.6315, -0.7795,  0.7795, -1.6315;    % y
                   -0.2887,  1.6799, -0.2887,  1.6799,  1.6799, -0.2887];   % z_wp
    c.pole_tip_x    = tips(1,:);
    c.pole_tip_y    = tips(2,:);
    c.pole_tip_z_wp = tips(3,:);

    % ===== ④ actuator frame（尖端到尖端 對連線 P1-P2/P3-P4/P5-P6；實測正交、= baseline magic-angle）=====
    a1 = tips(:,1)-tips(:,2);  a3 = tips(:,3)-tips(:,4);  a5 = tips(:,5)-tips(:,6);
    c.R_act = [a1/norm(a1), a3/norm(a3), a5/norm(a5)].';
    assert(abs(det(c.R_act)-1) < 1e-4, 'tip400 tip-to-tip R_act 非正確旋轉');
    assert(max(abs(c.R_act*c.R_act.' - eye(3)), [], 'all') < 1e-4, 'tip400 R_act 非正交');
    dhat      = tips ./ vecnorm(tips);
    c.Pc_base = c.R_act * dhat;                        % 電荷格：後退尖端方向進 actuator frame（非 canonical）

    % ===== ⑤ sensor 幾何（tip + 4.572·e2 + (gap+0.41)·n̂；上極 CONE_ANG、鈍尖 gap=rf(1−sinβ)）=====
    dir3     = @(el,az) [cos(el)*cos(az); cos(el)*sin(az); sin(el)];
    beta     = atan2(c.POLE_R, c.POLE_CONE_LEN);       % 半錐角 ~11.31°
    CONE_ANG = 35.4945*pi/180;                          % 上極正交修正後錐軸仰角
    rf = 400e-6;  gap = rf*(1 - sin(beta));  AIR = 0.41e-3;  SOFF = 4.572e-3;
    c.pole_axis  = zeros(3,6);
    c.sensor_pos = zeros(3,6);   c.sensor_n = zeros(3,6);
    for i = 1:6
        th = c.pole_angles(i)*pi/180;
        if c.pole_is_lower(i)
            e2 = dir3(-beta, th);          nh = dir3(-beta-pi/2, th);          ax = dir3(0, th);
        else
            e2 = dir3(CONE_ANG+beta, th);  nh = dir3(CONE_ANG+beta+pi/2, th);  ax = dir3(CONE_ANG, th);
        end
        c.sensor_pos(:,i) = tips(:,i) + SOFF*e2 + (gap+AIR)*nh;
        c.sensor_n(:,i)   = nh;
        c.pole_axis(:,i)  = ax;
    end

    % ===== ⑥ 方法旗標（main.m 讀，特例不寫死在 code）=====
    c.v_method     = 'scattered';    % sensor_local CSV ≠ solve mesh → 座標式對 solve 場取樣
    c.interp_to    = 'tip40um';      % WP 電荷 fit 內插到 tip40 的同一組點雲（公平比較）
    c.r_loc        = 0.6e-3;         % interp 源鄰域半徑
    c.sensor_r_loc = 1.5e-3;         % scattered sensor 鄰域半徑（與 r_loc 是兩個不同半徑）

    % ===== ⑦ 號誌 / 物理常數 =====
    c.s_source = [-1, +1, -1, +1, +1, -1];
    c.S_hall   = 130;
    c.N_c      = 70;
    c.mu_0     = 4*pi*1e-7;
    c.k_m      = 1e-7;
end
