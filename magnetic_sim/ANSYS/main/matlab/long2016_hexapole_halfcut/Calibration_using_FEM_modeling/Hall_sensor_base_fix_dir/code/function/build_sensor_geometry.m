function [sensor_pos, sensor_n, disc_u, disc_v, disc_local, Ndisc] = build_sensor_geometry(cnst)
% BUILD_SENSOR_GEOMETRY  每極 Hall-sensor 擺放 + Ø0.3mm disc 取樣格點。
%   下極=milled flat;上極=natural cone;n+ 出鋼。沿用 B_bar 慣例。
%   回傳:
%     sensor_pos 3x6  各極 sensor 全域位置 [m]
%     sensor_n   3x6  各極 sensor 法線 n+(出鋼)
%     disc_u/v   3x6  disc 平面正交基底
%     disc_local Ndisc x2  disc 同心環取樣(中心 + 15 環,共 721 點)
%     Ndisc      取樣點數
%   (ported verbatim from calib_fem.m PAGE 2)

    % ---- sensor 幾何(下極=milled flat;上極=natural cone;n+ 出鋼)----
    beta = atan2(3.0, 15.0);                                     % 錐半角 11.31°(ANSYS POLE_R 3.0 / CONE_LEN 15.0)
    sensor_pl_lower = [4.572e-3; 0.41e-3; 0];                    % 下極 pole-local:4.572 沿軸、0.41 ⊥ milled flat
    n_pl_lower      = [0; 1; 0];                                 % 下極法線 n+ = +up(⊥ 平面、出鋼)
    sensor_pl_upper = [4.572e-3*cos(beta) - 0.41e-3*sin(beta); ...% 上極 pole-local:4.572 沿錐斜邊
                       4.572e-3*sin(beta) + 0.41e-3*cos(beta); ...%   + 0.41 ⊥ 錐面外
                       0];
    n_pl_upper      = [-sin(beta); cos(beta); 0];               % 上極法線 n+ ⊥ 錐面(出鋼)

    sensor_pos = zeros(3,6); sensor_n = zeros(3,6);            % 全域 sensor 位置 / 法線
    disc_u = zeros(3,6); disc_v = zeros(3,6);                  % disc 平面兩個基底向量
    for i = 1:6                                                % 逐極 P1..P6
        th = cnst.pole_angles(i)*pi/180;                      % 方位角
        if cnst.pole_is_lower(i)                              % --- 下極 ---
            pole_axis = [cos(th); sin(th); 0];               % 極軸 = 水平徑向
            up_hat = [0;0;1]; tip_z = -cnst.R_norm_z;        % up = +z;尖端 z
            s_pl = sensor_pl_lower; n_pl = n_pl_lower;        % 下極 offset / 法線
        else                                                 % --- 上極 ---
            inc = cnst.upper_incline;                        % 上極傾角(~36.6°,沿用 mt_constants)
            pole_axis = [cos(inc)*cos(th); cos(inc)*sin(th); sin(inc)];  % 傾斜極軸
            up_un = [0;0;1] - sin(inc)*pole_axis; up_hat = up_un/norm(up_un);  % up = ẑ 去掉沿軸分量
            tip_z = +cnst.R_norm_z; s_pl = sensor_pl_upper; n_pl = n_pl_upper; % 上極 offset / 法線
        end
        tip  = [cnst.R_norm_xy*cos(th); cnst.R_norm_xy*sin(th); tip_z];  % 極尖座標(WP 框)
        side = cross(pole_axis, up_hat); side = side/norm(side);         % 側向(⊥ axis & up)
        Rg   = [pole_axis, up_hat, side];                    % pole-local → global 旋轉矩陣
        sensor_pos(:,i) = tip + Rg*s_pl;                     % sensor 全域位置
        sensor_n(:,i)   = Rg*n_pl;                           % sensor 全域法線 n+
        u = up_hat - dot(up_hat,sensor_n(:,i))*sensor_n(:,i);% disc 基底 u = up_hat 投到 ⊥n 平面
        if norm(u)<1e-9, u = pole_axis - dot(pole_axis,sensor_n(:,i))*sensor_n(:,i); end  % 退化時改用 axis
        u = u/norm(u); disc_u(:,i) = u; disc_v(:,i) = cross(sensor_n(:,i),u);  % 正交基底 u,v
    end

    % ---- Ø0.3mm Hall disc 取樣點(同心環,721 點;面積平均近似)----
    sensor_r = 0.15e-3; n_rings = 15; ring_dr = sensor_r/n_rings; % 半徑 0.15mm、15 環、環距 0.01mm
    disc_local = [0 0];                                          % 中心點
    for kk = 1:n_rings                                          % 每環點數 6k(正比半徑)
        rk = ring_dr*kk; nk = 6*kk; phi = (0:nk-1)*2*pi/nk;
        disc_local = [disc_local; rk*cos(phi(:)), rk*sin(phi(:))]; %#ok<AGROW>
    end
    Ndisc = size(disc_local,1);                                % = 721
end
