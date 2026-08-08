function [pos, nhat, geo] = pole_sensor_geometry(cfg, opt)
%POLE_SENSOR_GEOMETRY  六顆 Hall sensor 的中心位置與外法線 n+（WP frame、公尺）。
%
%   [pos, nhat, geo] = pole_sensor_geometry(cfg)
%   [pos, nhat, geo] = pole_sensor_geometry(cfg, opt)
%
%   **這是 sensor 幾何的唯一來源**（2026-08-08 定案）。build_V_matrix 與所有畫 sensor /
%   磁極輪廓的腳本都必須呼叫本函式，不得再自己複製一份幾何（過去有 11 份副本各自漂移）。
%   本函式只接收模型尺寸參數（cfg，來自 mt_constants）並算幾何，不碰場資料、不畫圖。
%
% ---------------------------------------------------------------------------
% 幾何建構（psi = 繞極軸的方位角，psi=0 為 canonical sensor 方位）
%
%     pos(psi) = T + s_ax*e2 + R(s_ax)*rhat(psi) + AIR*e3(psi)
%
%     T         極尖              = R_norm * e1
%     e1        原點 -> 極尖       dir(+-psi0, th)      psi0 = 35.2644 deg（magic angle 仰角）
%     e2        極軸              dir(el, th)          el = 0（下極）/ upper_incline（上極）
%     rhat(0)   徑向、朝貼附面側   dir(el + sg*90, th)  sg = -1（下極底錐面）/ +1（上極）
%     tvec      切向、垂直子午面   [-sin(th); cos(th); 0]
%     rhat(psi) = cos(psi)*rhat(0) + sin(psi)*tvec      （繞 e2 轉 psi）
%     e3(psi)   真實錐面外法線     = cos(beta)*rhat(psi) - sin(beta)*e2
%                                   psi=0 時 = dir(el + sg*(beta+90), th)，即自 e2 轉 sg*(90+beta)
%     s_ax      軸向站位          = (SOFF - t_tan)*cos(beta) + t_tan
%     R(s)      真實錐體半徑      = (e_v + s)*tan(beta)，  e_v = r_f/sin(beta) - r_f
%
% 為什麼不是舊構法：40 µm 尖端倒圓讓**虛擬錐頂**落在極尖之外 e_v ~ 0.16 mm，極尖所在
% 軸向位置的錐面半徑已經是 ~0.033 mm 而不是 0。舊構法把貼附點錨在理想極尖、沿
% dir(+-beta) 直接走 SOFF，貼附點會埋進鋼裡 -> 對**真實**錐面的垂距只有 0.3635 mm（下極）
% / 0.4012 mm（上極），而不是設定的 0.41 mm。新構法的垂距精確等於 opt.air。
%
% SOFF 語意 = **沿貼附面量的斜面距離**（定案 4.572e-3 m = 0.18 inch）。最前面 t_tan
% 那一段是倒圓、不在錐面上，所以只有 (SOFF - t_tan) 要乘 cos(beta)。
% t_tan = **CAD STEP 實測 0.0322 mm**（錐面與倒圓球的切點距極尖的軸向距離，六根極一致）。
% ⚠ 不是 0.04 —— 切點落在球面 90deg-beta 處而非赤道，故軸向只推進 r_f*(1-sin beta)。
%
% 下極 face_lower='flat'（平切上表面）：該平面**通過極軸** -> 無徑向外推、斜面距離即
% 軸向距離、n+ = +z；psi 對平面不適用（忽略）。
% ---------------------------------------------------------------------------
%
% 輸入
%   cfg  模型設定（mt_constants）。必要欄位：
%          pole_angles(1x6,deg)、pole_is_lower(1x6)、R_norm、upper_incline、POLE_TIP_R
%        真實錐體欄位（CAD STEP 實測；缺少時以名目值/理論式回退，並在 geo.src 標示）：
%          pole_cone_slope(1x2)   [下極, 上極] tan(beta)
%          pole_tip_axial         倒圓的軸向推進量 t_tan [m]
%          pole_cone_len_real(1x2)[下極, 上極] 錐長（自極尖沿軸）[m]，供畫輪廓
%   opt  可選 struct：
%          .soff_upper  上極沿貼附面斜面距離 [m]（預設 4.572e-3）
%          .soff_lower  下極同上           [m]（預設 4.572e-3）
%          .face_lower  'cone'（預設）| 'flat' | 'upper_cone'（填滿極的上錐面，對照用）
%          .air         離面氣隙 [m]（預設 0.41e-3）
%          .psi         繞極軸旋轉 [rad]，純量或 1x6（預設 0）
%
% 輸出
%   pos   3x6  sensor 中心（WP frame, m）
%   nhat  3x6  外法線 n+（單位向量）
%   geo   struct，供畫圖/診斷共用同一組幾何：
%          .tip 3x6 極尖T   .axis 3x6 e2   .rad0 3x6 rhat(0)   .tanv 3x6 tvec
%          .ehat_r 3x6 水平徑向[cos th;sin th;0]（2D 子午面投影用：v2d=[ehat_r'*v; v(3)]）
%          .beta 1x6 真實半錐角[rad]      .s_ax 1x6 軸向站位[m]
%          .R    1x6 該站位錐體半徑[m]    .apex_off 1x6 虛擬錐頂外移 e_v[m]
%          .t_tan 1x6 倒圓軸向推進[m]     .cone_len 1x6 真實錐長[m]
%          .r_tip 1x6 倒圓半徑[m]         .face 1x6 cellstr 'cone'|'flat'
%          .src  幾何來源標記 'cad' | 'nominal'
%
% 相關：.claude/rules/ansys-cad-alignment.md（CAD = source of truth）、
%       memory long2016-pole-cone-geometry（STEP 實測數值與量法）。

    if nargin < 2 || isempty(opt), opt = struct(); end
    gv = @(f, d) getfielddef(opt, f, d);
    SOFF_U = gv('soff_upper', 4.572e-3);
    SOFF_L = gv('soff_lower', 4.572e-3);
    FACE_L = gv('face_lower', 'cone');
    AIR    = gv('air',        0.41e-3);
    PSI    = gv('psi',        0);
    if isscalar(PSI), PSI = repmat(PSI, 1, 6); end

    psi0   = atan2(cfg.R_norm_z, cfg.R_norm_xy);   % magic-angle 仰角 ~35.2644 deg（只進 e1）
    inc_up = cfg.upper_incline;                    % 上極錐軸仰角 ~36.5895 deg
    ell    = cfg.R_norm;
    r_f    = cfg.POLE_TIP_R;

    % ---- 錐體幾何：優先用 CAD 實測，缺欄位才回退名目值 ----
    if isfield(cfg, 'pole_cone_slope')
        bet2 = atan(cfg.pole_cone_slope(:).');     % [下極, 上極]
        geo_src = 'cad';
    else
        bet2 = atan2(cfg.POLE_R, cfg.POLE_CONE_LEN)*[1 1];
        geo_src = 'nominal';
    end
    if isfield(cfg, 'pole_cone_len_real')
        len2 = cfg.pole_cone_len_real(:).';
    else
        len2 = cfg.POLE_CONE_LEN*[1 1];
    end

    dir = @(el, az) [cos(el)*cos(az); cos(el)*sin(az); sin(el)];
    pos = zeros(3,6);  nhat = zeros(3,6);
    geo.tip = zeros(3,6);  geo.axis  = zeros(3,6);  geo.rad0     = zeros(3,6);
    geo.tanv= zeros(3,6);  geo.ehat_r= zeros(3,6);  geo.beta     = zeros(1,6);
    geo.s_ax= zeros(1,6);  geo.R     = zeros(1,6);  geo.apex_off = zeros(1,6);
    geo.t_tan=zeros(1,6);  geo.cone_len = zeros(1,6);
    geo.r_tip = r_f*ones(1,6);  geo.face = cell(1,6);  geo.src = geo_src;  geo.psi = PSI;

    for i = 1:6
        th  = cfg.pole_angles(i)*pi/180;
        low = logical(cfg.pole_is_lower(i));
        if low
            el = 0;       sg = -1;  bet = bet2(1);  soff = SOFF_L;  cl = len2(1);
        else
            el = inc_up;  sg = +1;  bet = bet2(2);  soff = SOFF_U;  cl = len2(2);
        end
        flat   = low && strcmpi(FACE_L, 'flat');
        upcone = low && strcmpi(FACE_L, 'upper_cone');   % [ADDED 2026-08-08] 填滿極的**上**錐面
        if upcone, sg = +1; end                          %   （對稱於下錐面；用於 lower_filled 對照）

        e1   = dir((2*(~low)-1)*psi0, th);           % 下極 -psi0 / 上極 +psi0
        T    = ell*e1;
        ax   = dir(el, th);                          % e2 = 極軸
        rad0 = dir(el + sg*pi/2, th);                % 徑向、朝貼附面側
        tv   = [-sin(th); cos(th); 0];               % 切向、垂直子午面

        % t_tan：CAD 實測優先，否則用球切錐理論式 r_f*(1-sin beta)
        if isfield(cfg, 'pole_tip_axial'), t_tan = cfg.pole_tip_axial;
        else,                              t_tan = r_f*(1 - sin(bet)); end

        if flat
            % 平切上表面通過極軸 -> 斜面距離即軸向距離、無徑向外推、n+ = +z
            s_ax = soff;  R = 0;  e_v = 0;
            nh   = [0; 0; 1];
            pos(:,i) = T + s_ax*ax + AIR*nh;
        else
            s_ax = (soff - t_tan)*cos(bet) + t_tan;  % 斜面距離 -> 軸向站位
            e_v  = r_f/sin(bet) - r_f;               % 虛擬錐頂在極尖外的距離
            R    = (e_v + s_ax)*tan(bet);            % 該站位的真實錐體半徑
            rh   = cos(PSI(i))*rad0 + sin(PSI(i))*tv;
            nh   = cos(bet)*rh - sin(bet)*ax;        % = dir(el+sg*(bet+pi/2),th) 當 psi=0
            pos(:,i) = T + s_ax*ax + R*rh + AIR*nh;
        end
        nhat(:,i) = nh;

        geo.tip(:,i)   = T;    geo.axis(:,i)  = ax;   geo.rad0(:,i) = rad0;
        geo.tanv(:,i)  = tv;   geo.ehat_r(:,i) = [cos(th); sin(th); 0];
        geo.beta(i)    = bet;  geo.s_ax(i)    = s_ax; geo.R(i)      = R;
        geo.apex_off(i)= e_v;  geo.t_tan(i)   = t_tan; geo.cone_len(i) = cl;
        if flat, geo.face{i} = 'flat'; elseif upcone, geo.face{i} = 'upper_cone';
        else,    geo.face{i} = 'cone'; end
    end
end

% ---- local ----
function v = getfielddef(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
