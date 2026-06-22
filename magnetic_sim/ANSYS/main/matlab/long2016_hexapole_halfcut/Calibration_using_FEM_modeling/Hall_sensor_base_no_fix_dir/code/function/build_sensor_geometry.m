function [sensor_pos, sensor_n] = build_sensor_geometry(cnst)
% BUILD_SENSOR_GEOMETRY  6 顆 Hall sensor 的「中心位置 + 法線 n+」，全域(WP 中心)座標。
% -------------------------------------------------------------------------
% 用途：
%   電壓抽取(extract_Vmat)只需每顆 sensor 的中心點與法線 n+：沿 n+ 開圓柱挑真實 FEM
%   節點、對 B·n+ 平均得 sensor 電壓。本檔把 6 顆 sensor 直接用全域座標算出（無 local frame）。
%
% 兩種極（皆全域座標、共用單位向量 dir(e,a)=[cos(e)cos(a); cos(e)sin(a); sin(e)]）：
%   - 下極(P1/P3/P6，halfcut 磨平極)：磨平面在頂，sensor 沿水平徑向偏移 + 垂直抬升：
%       tip        = [R_norm_xy cosθ; R_norm_xy sinθ; −R_norm_z]
%       sensor_pos = tip + SOFF·[cosθ; sinθ; 0] + AIR·[0;0;1]
%       sensor_n   = [0;0;1]（n+ = +z，⊥磨平面、出鋼朝 WP）
%   - 上極(P2/P4/P5，自然錐面)：軸用 CAD 實際傾角 inc_up = cnst.upper_incline ≈ 36.59°（極尖→底端連線，非理想魔術角）：
%       tip        = [R_norm_xy cosθ; R_norm_xy sinθ; +R_norm_z]
%       沿錐面斜邊 slant = dir(inc_up+β, θ)；錐面外法線 nrm = dir(inc_up+β+90°, θ)（β=半錐角≈11.31°）
%       sensor_pos = tip + SOFF·slant + AIR·nrm ;  sensor_n = nrm（n+ 出鋼、指向 sensor）
%   SOFF = 4.572mm（沿表面距極尖）、AIR = 0.41mm（離面 air-gap）。
%
% 輸入：
%   cnst        幾何常數（R_norm_xy/R_norm_z、POLE_R、POLE_CONE_LEN、pole_angles、pole_is_lower）
% 輸出：
%   sensor_pos  3×6  各極 sensor 中心的全域座標 [m]（WP 框）
%   sensor_n    3×6  各極 sensor 法線 n+（單位向量，出鋼、指向 sensor）
% -------------------------------------------------------------------------
    beta   = atan2(cnst.POLE_R, cnst.POLE_CONE_LEN);         % 半錐角 β = atan(3/15) ≈ 11.31° [rad]（上極用）
    inc_up = cnst.upper_incline;                             % 上極軸傾角 = CAD 實際值 ≈ 36.59°（極尖→底端連線，非理想魔術角）[rad]
    SOFF  = 4.572e-3;                                        % sensor 沿表面距極尖的距離 [m]
    AIR   = 0.41e-3;                                         % sensor 離面 air-gap [m]
    dir = @(e,a) [cos(e)*cos(a); cos(e)*sin(a); sin(e)];     % 仰角 e、方位 a 的單位向量（球→直角）

    sensor_pos = zeros(3,6); sensor_n = zeros(3,6);          % 預配置：全域 sensor 位置 / 法線（3 座標 × 6 極）
    for i = 1:6                                              % 逐極 P1..P6
        th = cnst.pole_angles(i)*pi/180;                    % 該極方位角 θ [rad]
        cz = cos(th); sz = sin(th);                         % 方位 cos/sin
        if cnst.pole_is_lower(i)                            % ===== 下極：磨平面（n+ = +z）=====
            tip = [cnst.R_norm_xy*cz; cnst.R_norm_xy*sz; -cnst.R_norm_z];  % 下極尖端（z = −R_norm_z）
            sensor_n(:,i)   = [0;0;1];                                     % 法線 n+ = +z（⊥磨平面、出鋼朝 WP）
            sensor_pos(:,i) = tip + SOFF*[cz;sz;0] + AIR*[0;0;1];          % 沿水平徑向 SOFF + 垂直抬 AIR
        else                                               % ===== 上極：魔術角錐面 =====
            tip = [cnst.R_norm_xy*cz; cnst.R_norm_xy*sz; +cnst.R_norm_z];  % 上極尖端（z = +R_norm_z）
            slant = dir(inc_up+beta,        th);                           % 沿錐面斜邊（仰角 inc_up+β，朝 base/離 WP）
            nrm   = dir(inc_up+beta+pi/2,   th);                           % 錐面外法線（再轉 90°，出鋼、指向 sensor）
            sensor_n(:,i)   = nrm;
            sensor_pos(:,i) = tip + SOFF*slant + AIR*nrm;                  % 沿錐面 SOFF + 沿法線 AIR
        end
    end
end
