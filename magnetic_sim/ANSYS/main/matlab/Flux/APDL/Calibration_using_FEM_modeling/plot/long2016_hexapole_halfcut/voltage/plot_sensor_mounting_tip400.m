function plot_sensor_mounting_tip400(POLE, PREVIEW)
%PLOT_SENSOR_MOUNTING_TIP400  tip400µm（CNC 鈍尖）P1/P2 的 Hall sensor 安裝位置側視圖（x-z）。
%   幾何 = 正交修正後 scenario-A tip400（APDL deck 驗證）：tip 後退 + 上極錐軸 CONE_ANG=35.49°、
%   0.4mm 鈍尖倒圓。sensor：tip 沿錐面 e2 直走 4.572mm（直線）→ 補 gap=rf(1−sinβ) 到錐面 → 再出 0.41mm。
%   星星→十字、去「WP」字與「(WP frame)」、實際輪廓。
    if nargin < 1 || isempty(POLE),    POLE = 'P2';    end
    if nargin < 2 || isempty(PREVIEW), PREVIEW = true; end
    POLE = upper(POLE);

    % ---- tip400 幾何常數（WP frame, mm；由 MT_Geom_Export.txt 驗證）----
    beta     = atan2(3.0, 15.0);            % 半錐角 β ≈ 11.31°
    CONE_ANG = 35.4945*pi/180;              % 上極正交修正後錐軸仰角
    Lsl      = hypot(15.0, 3.0);            % 錐面斜長 ≈ 15.30 mm
    rf       = 0.4;                         % 尖端倒圓半徑 [mm]（=POLE_TIP_R 400µm）
    SOFF     = 4.572;  AIR = 0.41;          % 沿錐面直線距 / 離面 [mm]
    dir  = @(el,az) [cos(el)*cos(az); sin(el)];
    rot  = @(v,a)   [cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];

    switch POLE
        case 'P1'          % 0°、下極、半切、錐軸水平
            az = 0;  islow = true;
            T  = [1.8839; -0.2887];
            axc = dir(0, az);
            e2  = dir(-beta, az);                         % sensor 所在的下錐面 slant
            nh  = dir(-beta - pi/2, az);                  % n+（朝下出鋼）
        case 'P2'          % 180°、上極、全錐、錐軸 CONE_ANG
            az = pi;  islow = false;
            T  = [-0.9001; 1.6799];
            axc = dir(CONE_ANG, az);
            e2  = dir(CONE_ANG + beta, az);               % sensor 所在的錐面 slant
            nh  = dir(CONE_ANG + beta + pi/2, az);        % n+（出鋼）
        otherwise, error('POLE 必須是 P1 或 P2');
    end

    % ---- sensor：直線 4.572 → 補 gap 到錐面 → 出 0.41 ----
    foot_raw = T + SOFF*e2;                               % 直線端點（沿 e2 距 tip 4.572）
    gap      = rf*(1 - sin(beta));                        % 直線端點→錐面 垂距（鈍尖倒圓造成）≈0.322mm
    fsurf    = foot_raw + gap*nh;                         % 錐面上落腳點
    sensor   = fsurf + AIR*nh;                            % 貼面 + 0.41mm

    % ---- 磁極截面輪廓（後退 tip + 錐軸 + ±β 錐面 + 0.4mm 鈍尖倒圓）----
    C    = T + rf*axc;   axa = atan2(axc(2), axc(1));   angf = axa + pi;
    hw   = pi/2 - beta;
    fdA  = rot(axc, +beta);  fdB = rot(axc, -beta);
    if islow
        tha = linspace(angf, angf+hw, 24);  arc = C + rf*[cos(tha); sin(tha)];
        low = arc(:,end) + Lsl*fdB;   top = T + Lsl*axc;
        poly = [top, T, arc, arc(:,end), low];
    else
        tha = linspace(angf-hw, angf+hw, 40);  arc = C + rf*[cos(tha); sin(tha)];
        Cu = arc(:,1) + Lsl*fdA;   Cd = arc(:,end) + Lsl*fdB;
        poly = [Cu, arc, Cd];
    end

    %% ---- 繪圖 ----
    figure('Position',[60 60 980 900],'Color','w'); ax = axes; hold on;
    patch('XData',poly(1,:),'YData',poly(2,:),'FaceColor',[.80 .80 .84],'EdgeColor','k','LineWidth',2,'FaceAlpha',1);
    plab = T + 3.2*axc;
    text(plab(1),plab(2),sprintf(' %s pole',POLE),'FontSize',14,'FontWeight','bold','Color',[.3 .3 .35]);

    % WP（十字，無字）+ WP→tip 連線
    plot([0 T(1)],[0 T(2)],'-','Color',[.5 .5 .5],'LineWidth',1.4);
    plot(0,0,'+','MarkerSize',20,'Color','k','LineWidth',2.5);

    % 4.572mm（藍，直線 tip→foot_raw）
    plot([T(1) foot_raw(1)],[T(2) foot_raw(2)],'-','Color',[.15 .35 .75],'LineWidth',1.8);
    mB=(T+foot_raw)/2; text(mB(1),mB(2),'  4.572 mm','FontSize',12.5,'FontWeight','bold','Color',[.15 .35 .75]);

    % gap（紫，直線端點→錐面；鈍尖倒圓造成的垂距）
    plot([foot_raw(1) fsurf(1)],[foot_raw(2) fsurf(2)],'-','Color',[.55 .1 .6],'LineWidth',2.2);
    mg=(foot_raw+fsurf)/2; text(mg(1)+0.05,mg(2),sprintf(' %.2f mm',gap),'FontSize',11,'FontWeight','bold','Color',[.55 .1 .6]);
    plot(fsurf(1),fsurf(2),'.','MarkerSize',16,'Color','k');            % 錐面上落腳點

    % 0.41mm（棕，錐面→sensor 沿 n+）
    plot([fsurf(1) sensor(1)],[fsurf(2) sensor(2)],'-','Color',[.75 .35 .1],'LineWidth',1.8);
    text(sensor(1),sensor(2)-0.28,'0.41 mm','FontSize',12,'FontWeight','bold','Color',[.75 .35 .1]);

    % tip（黑方塊）
    plot(T(1),T(2),'s','MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','k');
    text(T(1),T(2),'  tip','FontSize',13,'FontWeight','bold');

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
    allx = [0 T(1) sensor(1) a1(1) foot_raw(1)];  allz = [0 T(2) sensor(2) a1(2) foot_raw(2)];
    pad = 0.7;  xlim([min(allx)-pad max(allx)+pad]);  ylim([min(allz)-pad max(allz)+pad]);
    daspect([1 1 1]);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]); box on; grid off;
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));  yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel('x (mm)','FontSize',17,'FontWeight','bold');  ylabel('z (mm)','FontSize',17,'FontWeight','bold');
    title(sprintf('%s Hall sensor 安裝位置（tip400\\mum 幾何示意）',POLE),'FontSize',15,'FontWeight','bold','Interpreter','tex');
    ax.Toolbar.Visible = 'off';

    fprintf('%s tip400: tip=(%.3f,%.3f) foot_raw=(%.3f,%.3f) gap=%.3f fsurf=(%.3f,%.3f) sensor=(%.3f,%.3f)\n', ...
            POLE, T, foot_raw, gap, fsurf, sensor);
    if PREVIEW
        out = fullfile(tempdir, sprintf('sensor_mounting_tip400_%s_preview.png',POLE));  exportgraphics(gcf,out,'Resolution',150);
    else
        fdir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\long2016_hexapole_halfcut\Calibration_using_FEM_modeling\voltage_base\figures\shared';
        out = fullfile(fdir, sprintf('sensor_mounting_tip400_%s.png',POLE));  exportgraphics(gcf,out,'Resolution',200);
    end
    fprintf('saved %s\n', out);
end
