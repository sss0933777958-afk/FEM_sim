function [sensor_pos, sensor_n] = build_sensor_geometry(cnst)
% BUILD_SENSOR_GEOMETRY  算出 6 顆 Hall sensor 的「中心位置 + 法線 n+」。
% -------------------------------------------------------------------------
% 用途：
%   電壓抽取(extract_Vmat, real-node 版)只需要每顆 sensor 的中心點與法線：
%   以中心、沿 n+ 開一個圓柱去挑真實 FEM 節點、對 B·n+ 平均得 sensor 電壓。
%   （舊版還會回傳 disc 取樣格點 disc_u/v/local/Ndisc 給「內插盤面平均」用，
%     real-node 圓柱選點不需要，故已移除。）
%
% 擺放慣例（沿用 B_bar 矩陣設定）：
%   - 下極(P1/P3/P6)：尖端磨平(milled flat)，sensor 貼平面、n+ ⊥平面朝外(出鋼)。
%   - 上極(P2/P4/P5)：保留錐面(natural cone)，sensor 貼錐面、n+ ⊥錐面朝外。
%   - sensor 中心 = 距極尖 4.572mm（沿極軸/錐斜邊）、再離面 0.41mm(air-gap)。
%   - n+ 一律「出鋼」方向（決定 B·n+ 的正負號慣例，見 extract_Vmat）。
%
% 輸入：
%   cnst        幾何常數（pole_angles / pole_is_lower / upper_incline / R_norm_xy / R_norm_z）
% 輸出：
%   sensor_pos  3×6  各極 sensor 中心的全域座標 [m]（WP 框）
%   sensor_n    3×6  各極 sensor 法線 n+（單位向量，出鋼）
%
%   (改寫自 calib_fem.m PAGE 2：保留 pose 計算、移除 disc 取樣)
% -------------------------------------------------------------------------

    % ---- sensor pole-local 偏移與法線（下極=milled flat；上極=natural cone）----
    beta = atan2(3.0, 15.0);                                     % 錐半角 β = atan(POLE_R/CONE_LEN) = atan(3.0/15.0) ≈ 11.31°
    sensor_pl_lower = [4.572e-3; 0.41e-3; 0];                    % 下極 pole-local 偏移：沿軸 4.572mm、垂直磨平面 0.41mm
    n_pl_lower      = [0; 1; 0];                                 % 下極 pole-local 法線 n+ = +「up」（⊥磨平面、出鋼）
    sensor_pl_upper = [4.572e-3*cos(beta) - 0.41e-3*sin(beta); ...% 上極 pole-local 偏移：4.572 沿錐斜邊 …
                       4.572e-3*sin(beta) + 0.41e-3*cos(beta); ...%   … + 0.41 沿錐面外法線（把(沿軸,離面)旋轉 β 到錐座標）
                       0];
    n_pl_upper      = [-sin(beta); cos(beta); 0];               % 上極 pole-local 法線 n+（⊥錐面、出鋼）

    sensor_pos = zeros(3,6); sensor_n = zeros(3,6);            % 預配置：全域 sensor 位置 / 法線
    for i = 1:6                                                % 逐極 P1..P6 建立 pole-local→global 座標
        th = cnst.pole_angles(i)*pi/180;                      % 該極方位角 θ [rad]
        if cnst.pole_is_lower(i)                              % ===== 下極 =====
            pole_axis = [cos(th); sin(th); 0];               % 極軸 = 水平徑向（xy 平面內）
            up_hat = [0;0;1]; tip_z = -cnst.R_norm_z;        % 「up」取 +z；下極尖端 z = −R_norm_z
            s_pl = sensor_pl_lower; n_pl = n_pl_lower;        % 採用下極本地偏移 / 法線
        else                                                 % ===== 上極 =====
            inc = cnst.upper_incline;                        % 上極傾角（~36.6°，由 mt_constants 提供）
            pole_axis = [cos(inc)*cos(th); cos(inc)*sin(th); sin(inc)];  % 傾斜極軸（含 +z 分量）
            up_un = [0;0;1] - sin(inc)*pole_axis; up_hat = up_un/norm(up_un);  % 「up」= ẑ 去掉沿軸分量後正規化（⊥極軸）
            tip_z = +cnst.R_norm_z; s_pl = sensor_pl_upper; n_pl = n_pl_upper; % 上極尖端 z = +R_norm_z；採上極偏移/法線
        end
        tip  = [cnst.R_norm_xy*cos(th); cnst.R_norm_xy*sin(th); tip_z];  % 極尖全域座標（WP 框）
        side = cross(pole_axis, up_hat); side = side/norm(side);         % 第三軸 side = axis×up（⊥兩者），組右手系
        Rg   = [pole_axis, up_hat, side];                    % pole-local→global 旋轉矩陣（欄 = [軸, up, side]）
        sensor_pos(:,i) = tip + Rg*s_pl;                     % sensor 全域中心 = 極尖 + 旋轉後的本地偏移
        sensor_n(:,i)   = Rg*n_pl;                           % sensor 全域法線 n+ = 旋轉後的本地法線
    end
end
