function [sensor_pos, sensor_n] = build_sensor_geometry(cnst, SOFF_upper)
% BUILD_SENSOR_GEOMETRY  6 顆 Hall sensor 的「中心位置 + 法線 n+」，全域(WP 中心)座標。
% -------------------------------------------------------------------------
% 用途：
%   電壓抽取(extract_Vmat)只需每顆 sensor 的中心點與法線 n+：沿 n+ 開圓柱挑真實 FEM
%   節點、對 B·n+ 平均得 sensor 電壓。本檔把 6 顆 sensor 直接用全域座標算出（無 local frame）。
%
% 嚴格照文件 doc/sensor_location/pdf/hall_sensor_position.pdf（相對 WP）：
%   ψ = 仰角（極尖相對 WP 的仰角）= atan2(R_norm_z, R_norm_xy) ≈ 35.26°（magic-angle 鎖定）；
%       上極 +ψ、下極 −ψ。β = 磁極半錐角 = atan2(POLE_R, POLE_CONE_LEN) ≈ 11.31°。
%       θ = 方位角。ℓ̂ = 工作空間半徑 = R_norm。dir(e,a)=[cos(e)cos(a);cos(e)sin(a);sin(e)]。
%   e1 = dir(ψ,θ)（→極尖）；e2 = dir(ψ+β,θ)（極尖→sensor，沿錐面）；
%   e3 = dir(ψ+β+90°,θ)（上極法線 n+，z=+cos(ψ+β)）。
%   上極（自然斜錐）：e2 = dir(ψ+β)、n+ = e3；P = ℓ̂·e1 + 4.572·e2 + 0.41·e3。
%   下極（FEM 水平半切 → 錐軸水平）：**e2/外法線改用 −β（不綁 ψ）**：
%     e2 = dir(−β,θ)（底錐面 slant）、n+ = dir(−β−90°,θ) = [−sinβ cosθ; −sinβ sinθ; −cosβ]（朝下出鋼）；
%     P = ℓ̂·e1 + 4.572·dir(−β) + 0.41·dir(−β−90°)。ψ 只進 e1（極尖在 magic-angle 位置 z=−R_norm_z）。
%   SOFF = 4.572mm（沿錐面距極尖）、AIR = 0.41mm（離面）。
%   說明：PDF 原式把錐面綁在 ψ，但 FEM 下極半切後錐軸是水平的 → 下極改用 −β 才使 sensor 真正
%     離實體底錐面 0.41mm（修正 PDF 的 ψ+β）。上極 PDF(ψ=35.26°) vs FEM 錐軸(36.59°) 差 ~0.1mm，暫不改。
%
% 輸入：cnst（R_norm/R_norm_xy/R_norm_z、POLE_R、POLE_CONE_LEN、pole_angles、pole_is_lower）
%       SOFF_upper（選填）：上極沿錐面距極尖 [m]，預設 4.572e-3；可外移上極 sensor（下極固定 4.572e-3）
% 輸出：sensor_pos 3×6 全域(WP 框)[m]；sensor_n 3×6 法線 n+（出鋼、指向 sensor）
% -------------------------------------------------------------------------
    beta = atan2(cnst.POLE_R, cnst.POLE_CONE_LEN);          % 半錐角 β ≈ 11.31° [rad]
    psi0 = atan2(cnst.R_norm_z, cnst.R_norm_xy);            % 仰角 ψ ≈ 35.26° [rad]（= WP→tip 仰角，magic-angle）
    ell  = cnst.R_norm;                                     % ℓ̂ = 工作空間半徑 [m]
    if nargin < 2 || isempty(SOFF_upper), SOFF_upper = 4.572e-3; end  % [ADDED] 上極沿錐面距極尖（選填，預設原值 4.572mm）
    SOFF_lower = 4.572e-3;                                  % 下極沿錐面距極尖 [m]（固定）
    AIR  = 0.41e-3;                                         % 離面 air-gap [m]
    dir = @(e,a) [cos(e)*cos(a); cos(e)*sin(a); sin(e)];    % 仰角 e、方位 a 的單位向量

    sensor_pos = zeros(3,6); sensor_n = zeros(3,6);
    for i = 1:6                                             % 逐極 P1..P6
        th = cnst.pole_angles(i)*pi/180;                   % 方位角 θ
        if cnst.pole_is_lower(i), psi = -psi0; else, psi = +psi0; end   % 下極 −ψ、上極 +ψ
        e1 = dir(psi, th);                                 % → 極尖（ψ=±仰角）
        if cnst.pole_is_lower(i)
            % 下極：FEM 為水平半切 → 錐軸水平，底錐面 slant/外法線用 −β（不綁 ψ）
            e2   = dir(-beta,      th);                     % 底錐面 slant（水平錐軸往下 β）
            nhat = dir(-beta-pi/2, th);                    % 底錐面外法線 = [−sinβ cz; −sinβ sz; −cosβ]（n+ 朝下出鋼）
        else
            % 上極：自然斜錐，沿錐面 ψ+β、外法線 ψ+β+90°（PDF）
            e2   = dir(psi+beta,      th);                 % 極尖→sensor（沿錐面 slant）
            nhat = dir(psi+beta+pi/2, th);                 % 錐面外法線（n+ 朝上出鋼）
        end
        if cnst.pole_is_lower(i), soff = SOFF_lower; else, soff = SOFF_upper; end   % [ADDED] 上極可沿斜面外移
        sensor_pos(:,i) = ell*e1 + soff*e2 + AIR*nhat;     % P = ℓ̂·e1 + soff·e2 + 0.41·n+
        sensor_n(:,i)   = nhat;
    end
end
