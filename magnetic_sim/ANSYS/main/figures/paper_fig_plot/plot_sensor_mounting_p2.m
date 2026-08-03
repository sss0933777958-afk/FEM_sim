function plot_sensor_mounting_p2(SOFF)
%   [ADDED] SOFF = 藍線長度（tip→foot 沿錐面直線距, mm）；省略 = 4.572（定案值）。
%   非預設值時輸出檔名自動加 _soff<值>mm 後綴，不覆蓋原圖。
% plot_sensor_mounting_p2 -- Section3_A paper figure: P2 Hall-sensor mounting side
% view (x-z), tip40µm geometry (幾何取自 mt_constants，非硬寫)。Clean paper style:
%   每軸 3 根等距 tick、tick 朝外、無尺寸數字、無軸標題、無 "P2 pole"/"sensor" 字、
%   感測器畫成實際圓柱(R0.15×H0.10mm)邊視、n+ 紅箭頭(小頭、無字)、字體 36。
%   藍線 tip→foot(沿錐面 4.572mm) + 橘線 foot→sensor(0.41mm，直接相連；tip40 不含 gap 段)。
%   輸出 → figures/paper_fig/Section3_A/sensor_mounting_tip40_P2.png
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here),'paper_fig','Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants

    % ---- tip40 geometry (P2, WP frame, mm; = plot_sensor_mounting.m) ----
    cnst   = mt_constants();
    beta   = atan2(cnst.POLE_R, cnst.POLE_CONE_LEN);        % 半錐角 β ≈ 11.31°
    inc_up = cnst.upper_incline;                            % 上極錐軸仰角 ≈ 36.59°
    Lsl    = hypot(cnst.POLE_CONE_LEN, cnst.POLE_R)*1e3;    % 錐面斜長 ≈ 15.30 mm
    rf     = cnst.POLE_TIP_R*1e3;                           % 尖端倒圓半徑 [mm]（40µm = 0.04）
    if nargin < 1 || isempty(SOFF), SOFF = 4.572; end       % [MODIFIED] 沿錐面直線距(藍線長)可調
    AIR = 0.41;                                             % 離面 [mm]
    dir = @(el,az)[cos(el)*cos(az); sin(el)];
    rot = @(v,a)[cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];

    az = pi;                                                % P2:180°、上極、全錐
    T   = [cnst.pole_tip_x(2); cnst.pole_tip_z_wp(2)]*1e3;
    axc = dir(inc_up, az);
    e2  = dir(inc_up + beta, az);
    nh  = dir(inc_up + beta + pi/2, az);                    % n+（出鋼）

    foot   = T + SOFF*e2;                                   % 錐面上落腳點（直線端點）
    sensor = foot + AIR*nh;                                 % 貼面 + 0.41mm（圓柱底面中心）

    % ---- 磁極截面輪廓（tip + 錐軸 + ±β 錐面 + 40µm 鈍尖倒圓）----
    C   = T + rf*axc;  axa = atan2(axc(2),axc(1));  angf = axa + pi;  hw = pi/2 - beta;
    fdA = rot(axc,+beta);  fdB = rot(axc,-beta);
    tha = linspace(angf-hw,angf+hw,60);  arc = C + rf*[cos(tha); sin(tha)];
    Cu  = arc(:,1) + Lsl*fdA;  Cd = arc(:,end) + Lsl*fdB;
    poly = [Cu, arc, Cd];

    %% ---- 繪圖 ----
    RED = [.92 .15 .15];  GRN = [.10 .70 .25];
    figure('Position',[60 30 980 980],'Color','w'); ax = axes; hold on;
    patch('XData',poly(1,:),'YData',poly(2,:),'FaceColor',[.80 .80 .84], ...
          'EdgeColor','k','LineWidth',2.5,'FaceAlpha',1);

    % WP（十字，無字）+ WP→tip 連線
    plot([0 T(1)],[0 T(2)],'-','Color',[.5 .5 .5],'LineWidth',1.6);
    plot(0,0,'+','MarkerSize',22,'Color','k','LineWidth',3);

    % 藍線 tip→foot（沿錐面 4.572mm）+ 橘線 foot→sensor（0.41mm，直接相連、無數字）
    plot([T(1) foot(1)],[T(2) foot(2)],'-','Color',[.15 .35 .75],'LineWidth',3);
    plot([foot(1) sensor(1)],[foot(2) sensor(2)],'-','Color',[.75 .35 .1],'LineWidth',3);

    % tip（黑方塊，無字）
    plot(T(1),T(2),'s','MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','k');

    % sensor：實際圓柱邊視（2R=0.30 × H=0.10，軸沿 n+），綠填黑邊（無字）
    Rs = 0.15;  Hs = 0.10;  tdir = [-nh(2); nh(1)];
    base = sensor;  top = sensor + Hs*nh;
    rect = [base+Rs*tdir, base-Rs*tdir, top-Rs*tdir, top+Rs*tdir];
    patch('XData',rect(1,:),'YData',rect(2,:),'FaceColor',GRN,'EdgeColor','k','LineWidth',1.4);

    % n+：紅箭頭（桿 + 對齊 nh 的小三角頭，無字）
    NLEN = 1.0;  a1 = sensor + NLEN*nh;
    hl = 0.20;  hw2 = 0.10;  bk = a1 - hl*nh;
    plot([sensor(1) bk(1)],[sensor(2) bk(2)],'-','Color',RED,'LineWidth',4.2);
    patch([a1(1) bk(1)+hw2*tdir(1) bk(1)-hw2*tdir(1)], ...
          [a1(2) bk(2)+hw2*tdir(2) bk(2)-hw2*tdir(2)], RED,'EdgeColor',RED);

    % ---- 軸：等比、box、每軸 3 根等距 tick、tick 朝外、字體 36、無軸標題 ----
    allx = [0 T(1) sensor(1) foot(1) a1(1)];
    allz = [0 T(2) sensor(2) a1(2)];
    pad = 0.7;  xlim([min(allx)-pad max(allx)+pad]);  ylim([min(allz)-pad max(allz)+pad]);
    daspect([1 1 1]);  box on;  grid off;
    set(ax,'XTick',ticks3(xlim),'YTick',ticks3(ylim));
    set(ax,'FontSize',36,'FontWeight','bold','LineWidth',3.5,'TickLength',[.02 .02],'TickDir','out');
    ax.Toolbar.Visible = 'off';

    sfx = ''; if abs(SOFF-4.572) > 1e-9, sfx = sprintf('_soff%gmm', SOFF); end   % [ADDED] 非預設值另存
    out = fullfile(figdir, sprintf('sensor_mounting_tip40_P2%s.png', sfx));
    exportgraphics(gcf,out,'Resolution',200);
    fprintf('saved %s\n',out);
end

% ---- 每軸 3 根等距、nice-step tick（含在範圍內）----
function tk = ticks3(lm)
    s0 = (lm(2)-lm(1))/3;
    p  = 10^floor(log10(s0));  c = [1 2 3 4 5 10]*p;
    [~,i] = min(abs(c-s0));  s = c(i);
    ctr = round((lm(1)+lm(2))/2/s)*s;
    tk  = ctr + [-1 0 1]*s;
end
