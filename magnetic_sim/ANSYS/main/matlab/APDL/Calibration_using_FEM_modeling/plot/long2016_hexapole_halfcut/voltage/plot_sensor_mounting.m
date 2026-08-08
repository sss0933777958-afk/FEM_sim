function plot_sensor_mounting(POLE, PREVIEW)
%PLOT_SENSOR_MOUNTING  tip40µm（設計 40µm 尖）P1/P2 的 Hall sensor 安裝位置側視圖（x-z）。
%   風格同 plot_sensor_mounting_tip400.m（十字 WP、無「WP」字、軸標無「(WP frame)」、實際輪廓 +
%   0.04mm 鈍尖倒圓弧）；幾何 = baseline magic-angle 尖端（0.5mm 共球）+ 上極錐軸 upper_incline≈36.59°。
%   sensor = tip 沿錐面 e2 直走 4.572mm → 再出 0.41mm（= build_sensor_geometry；tip40 的 gap=rf(1−sinβ)
%   ≈32µm 可忽略、校正不含此修正，故不畫 gap 段，與 tip400 的差別在此）。
%   tip40 倒圓僅 40µm、比 tip400 小 10×、全景近乎尖 → 右下角加**放大 inset** 顯示尖端倒圓。
%   幾何取自 mt_constants()：不硬寫 sensor 座標/法線。
    if nargin < 1 || isempty(POLE),    POLE = 'P2';    end
    if nargin < 2 || isempty(PREVIEW), PREVIEW = true; end
    POLE = upper(POLE);

    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）；sensor 幾何改由
    %   utils/pole_sensor_geometry 供給（唯一來源、CAD 實測錐體 + 真實氣隙 0.41mm）。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'));
    cnst = model_config('long2016_hexapole_halfcut','tip40um');
    [SPOS, SNOR, GEO] = pole_sensor_geometry(cnst);
    inc_up = cnst.upper_incline;                            % 上極錐軸仰角 ≈ 36.59°（baseline）
    rf     = cnst.POLE_TIP_R*1e3;                           % 尖端倒圓半徑 [mm]（40µm = 0.04）
    AIR    = 0.41;                                          % 離面 [mm]
    dir  = @(el,az) [cos(el)*cos(az); sin(el)];
    rot  = @(v,a)   [cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];

    switch POLE
        case 'P1'          % 0°、下極、半切、錐軸水平
            ip = 1;  az = 0;   islow = true;   axc = dir(0, az);
        case 'P2'          % 180°、上極、全錐、錐軸 upper_incline
            ip = 2;  az = pi;  islow = false;  axc = dir(inc_up, az);
        otherwise, error('POLE 必須是 P1 或 P2');
    end

    % ---- sensor / 輪廓幾何：全部取自共用 geo（子午面 y=0 → 取 (x,z)）----
    T      = [cnst.pole_tip_x(ip); cnst.pole_tip_z_wp(ip)]*1e3;
    nh     = SNOR([1 3],ip);                              % n+（外法線出鋼）
    sensor = SPOS([1 3],ip)*1e3;                          % 圓柱底面中心（離貼附面精確 0.41mm）
    foot   = sensor - AIR*nh;                             % 貼附面上落腳點
    beta   = GEO.beta(ip);                                % per-pole 真實半錐角（CAD）
    Lsl    = (GEO.cone_len(ip) - GEO.t_tan(ip))*1e3/cos(beta);   % 自切點沿母線到錐底的斜長

    % ---- 磁極截面輪廓（tip + 錐軸 + ±β 錐面 + 鈍尖倒圓）----
    C    = T + rf*axc;   axa = atan2(axc(2), axc(1));   angf = axa + pi;
    hw   = pi/2 - beta;
    fdA  = rot(axc, +beta);  fdB = rot(axc, -beta);
    if islow
        tha = linspace(angf, angf+hw, 40);  arc = C + rf*[cos(tha); sin(tha)];
        low = arc(:,end) + Lsl*fdB;   top = T + Lsl*axc;
        poly = [top, T, arc, arc(:,end), low];
    else
        tha = linspace(angf-hw, angf+hw, 60);  arc = C + rf*[cos(tha); sin(tha)];
        Cu = arc(:,1) + Lsl*fdA;   Cd = arc(:,end) + Lsl*fdB;
        poly = [Cu, arc, Cd];
    end

    %% ---- 繪圖（主視圖）----
    figure('Position',[60 60 980 900],'Color','w'); ax = axes; hold on;
    patch('XData',poly(1,:),'YData',poly(2,:),'FaceColor',[.80 .80 .84],'EdgeColor','k','LineWidth',2,'FaceAlpha',1);
    plab = T + 3.2*axc;
    text(plab(1),plab(2),sprintf(' %s pole',POLE),'FontSize',14,'FontWeight','bold','Color',[.3 .3 .35]);

    % WP（十字，無字）+ WP→tip 連線
    plot([0 T(1)],[0 T(2)],'-','Color',[.5 .5 .5],'LineWidth',1.4);
    plot(0,0,'+','MarkerSize',20,'Color','k','LineWidth',2.5);

    % 4.572mm（藍，直線 tip→foot）
    plot([T(1) foot(1)],[T(2) foot(2)],'-','Color',[.15 .35 .75],'LineWidth',1.8);
    mB=(T+foot)/2; text(mB(1),mB(2),'  4.572 mm','FontSize',12.5,'FontWeight','bold','Color',[.15 .35 .75]);

    % 0.41mm（棕，錐面→sensor 沿 n+）
    plot([foot(1) sensor(1)],[foot(2) sensor(2)],'-','Color',[.75 .35 .1],'LineWidth',1.8);
    m4=(foot+sensor)/2; text(m4(1)+0.12,m4(2),'0.41 mm','FontSize',12,'FontWeight','bold','Color',[.75 .35 .1]);
    plot(foot(1),foot(2),'.','MarkerSize',16,'Color','k');            % 錐面上落腳點

    % tip（小方塊，避免蓋住倒圓）
    plot(T(1),T(2),'s','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','k');
    text(T(1),T(2)+0.02,'  tip','FontSize',13,'FontWeight','bold');

    % sensor（綠圓 + 綠橫桿）
    tdir = [-nh(2); nh(1)]; tdir = tdir/norm(tdir);  DISC = 0.28;
    dd = sensor - DISC*tdir; de = sensor + DISC*tdir;
    plot([dd(1) de(1)],[dd(2) de(2)],'-','Color',[.1 .5 .15],'LineWidth',4);
    plot(sensor(1),sensor(2),'o','MarkerSize',9,'MarkerFaceColor',[.1 .7 .25],'MarkerEdgeColor','k','LineWidth',1);
    text(sensor(1),sensor(2),'  sensor','FontSize',13,'FontWeight','bold','Color',[.05 .4 .12]);

    % n+（紅箭頭）
    NLEN = 1.0;  a1 = sensor + NLEN*nh;
    plot([sensor(1) a1(1)],[sensor(2) a1(2)],'-','Color','w','LineWidth',6);
    plot([sensor(1) a1(1)],[sensor(2) a1(2)],'-','Color',[.92 .15 .15],'LineWidth',3.8);
    plot(a1(1),a1(2),'^','MarkerSize',14,'MarkerFaceColor',[.92 .15 .15],'MarkerEdgeColor','w');
    text(a1(1),a1(2),'  n_+','FontSize',16,'FontWeight','bold','Color',[.92 .15 .15]);

    % 軸（無「(WP frame)」；等比；風格①）
    allx = [0 T(1) sensor(1) a1(1) foot(1)];  allz = [0 T(2) sensor(2) a1(2) foot(2)];
    pad = 0.7;  xlim([min(allx)-pad max(allx)+pad]);  ylim([min(allz)-pad max(allz)+pad]);
    daspect([1 1 1]);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]); box on; grid off;
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));  yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel('x (mm)','FontSize',17,'FontWeight','bold');  ylabel('z (mm)','FontSize',17,'FontWeight','bold');
    title(sprintf('%s Hall sensor 安裝位置（tip40\\mum 幾何示意）',POLE),'FontSize',15,'FontWeight','bold','Interpreter','tex');
    ax.Toolbar.Visible = 'off';

    %% ---- 放大 inset：尖端 40µm 倒圓（全景看不見，故放大）----
    if strcmp(POLE,'P1'), ipos = [0.185 0.145 0.30 0.30]; else, ipos = [0.60 0.145 0.30 0.30]; end
    axi = axes('Position', ipos); hold(axi,'on');
    patch(axi,'XData',poly(1,:),'YData',poly(2,:),'FaceColor',[.80 .80 .84],'EdgeColor','k','LineWidth',2);
    plot(axi, arc(1,:), arc(2,:),'-','Color',[.85 .2 .1],'LineWidth',3);   % 倒圓弧強調（紅）
    plot(axi, T(1), T(2),'s','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k');
    w = 0.09;  xlim(axi,[T(1)-w T(1)+w]);  ylim(axi,[T(2)-w T(2)+w]);
    daspect(axi,[1 1 1]);  box(axi,'on');  set(axi,'XTick',[],'YTick',[],'LineWidth',1.5,'Color','w');
    title(axi, sprintf('tip 倒圓  r_f = %.0f \\mum', rf*1e3),'FontSize',11,'FontWeight','bold','Interpreter','tex');

    fprintf('%s tip40: tip=(%.3f,%.3f) foot=(%.3f,%.3f) sensor=(%.3f,%.3f)  [rf=%.0fum]\n', ...
            POLE, T, foot, sensor, rf*1e3);
    if PREVIEW
        out = fullfile(tempdir, sprintf('sensor_mounting_%s_preview.png',POLE));  exportgraphics(gcf,out,'Resolution',150);
    else
        % [MODIFIED 2026-08-08] 原路徑指向已刪除的舊 per-model 樹 → 改共用夾的 figures 佈局
        fdir = fullfile(CALROOT,'figures','long2016_hexapole_halfcut','voltage','common');
        if ~exist(fdir,'dir'), mkdir(fdir); end
        out = fullfile(fdir, sprintf('sensor_mounting_%s.png',POLE));  exportgraphics(gcf,out,'Resolution',200);
    end
    fprintf('saved %s\n', out);
end
