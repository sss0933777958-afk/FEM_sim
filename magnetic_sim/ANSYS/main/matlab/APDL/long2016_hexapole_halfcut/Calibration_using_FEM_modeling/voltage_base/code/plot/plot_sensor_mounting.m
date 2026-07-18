function plot_sensor_mounting(POLE, PREVIEW)
% PLOT_SENSOR_MOUNTING  Hall sensor 安裝位置「純幾何示意圖」（單極，x-z 側視 / WP 框）。
% -------------------------------------------------------------------------
% 目的：對照 reference/hall_sensor_full.pdf 的安裝式
%       P = ℓ̂·e + λ1·e1 + λ2·e2  （λ1=4.572mm 沿錐面、λ2=0.41mm air-gap），
%   把「磁極 + sensor」的實際安裝幾何畫出來，驗證 build_sensor_geometry 是否符合 PDF。
%   純幾何（無 FEM 場）；sensor 表示法同 plot_P2sensor_Braw_P1exc（綠圓盤 + 綠 disc 短線 + 紅 n+ 三角箭頭）。
%
% POLE     'P1'（下極，半切、θ=0°）或 'P2'（上極，自然斜錐、θ=180°）。皆落在 y=0 的 x-z 面。
% PREVIEW  true=落 tempdir 預覽（150 DPI）；false/省略=存 figures/sensor_mounting_<POLE>.png（200 DPI）。
%
% 幾何全部來自 mt_constants() + build_sensor_geometry(cnst)：不硬寫 sensor 座標/法線。
% -------------------------------------------------------------------------
    if nargin < 1 || isempty(POLE),    POLE = 'P2';    end
    if nargin < 2 || isempty(PREVIEW), PREVIEW = true; end
    POLE = upper(POLE);
    switch POLE, case 'P1', ip = 1; case 'P2', ip = 2; otherwise, error('POLE 必須是 P1 或 P2'); end

    %% ---- paths + 常數 / 幾何 ----
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\voltage_base\code\main_function']);  % build_sensor_geometry
    cnst = mt_constants();
    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);              % WP 框 [m]

    beta = atan2(cnst.POLE_R, cnst.POLE_CONE_LEN);                     % 半錐角 β ≈ 11.31°
    psi  = atan2(cnst.R_norm_z, cnst.R_norm_xy);                       % 仰角 ψ ≈ 35.26°
    Lsl  = hypot(cnst.POLE_CONE_LEN, cnst.POLE_R)*1e3;                 % 錐面斜長 [mm]
    rot  = @(v,a)[cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];   % 2D 旋轉

    T  = [cnst.pole_tip_x(ip), cnst.pole_tip_z_wp(ip)]*1e3;           % 極尖 (x,z) [mm]
    sp = sensor_pos(:,ip);  ni = sensor_n(:,ip);                     % sensor 中心/法線 [m]
    sx = sp(1)*1e3; sz = sp(3)*1e3;                                  % sensor (x,z) [mm]
    nx = ni(1); nz = ni(3); nn = hypot(nx,nz); nx=nx/nn; nz=nz/nn;   % n+ 投 x-z 正規化
    foot = [sx;sz] - (0.41)*[nx;nz];                                 % 錐面上的落腳點 = sensor − 0.41mm·n+ [mm]

    %% ---- 磁極截面錐輪廓（淡灰填色 + 黑邊）----
    if ip == 2   % 上極：自然斜錐，pole axis=[-cos(inc);sin(inc)]（θ=180°），兩條 ±β 斜邊
        axp = [-cos(cnst.upper_incline); sin(cnst.upper_incline)];
        Ct = T(:)+Lsl*rot(axp,+beta);  Cb = T(:)+Lsl*rot(axp,-beta);
        poly = [T(:), Ct, Cb];                                        % 三角錐截面
    else         % 下極：FEM 水平半切 → 平頂(+x 水平) + 下斜邊(-β)，錐軸水平(θ=0°)
        top1 = T(:)+Lsl*[1;0];  low1 = T(:)+Lsl*rot([1;0],-beta);
        poly = [T(:), top1, low1];                                    % 半切楔形截面
    end

    %% ---- 圖 ----
    fig = figure('Position',[60 60 980 900],'Color','w');
    ax  = axes(fig); hold(ax,'on');
    patch('XData',poly(1,:),'YData',poly(2,:),'FaceColor',[.80 .80 .84], ...
          'EdgeColor','k','LineWidth',2,'FaceAlpha',1);               % 磁極實體
    if ip==2, axisdir=[-cos(cnst.upper_incline);sin(cnst.upper_incline)]; else, axisdir=[cos(-beta/2);sin(-beta/2)]; end
    plc = T(:) + 0.42*norm([sx;sz]-0.41*[nx;nz]-T(:))*axisdir;        % 錐體內部可見點（軸向 ~42%）
    text(plc(1),plc(2),sprintf('%s pole',POLE),'FontSize',14,'FontWeight','bold', ...
         'Color',[.3 .3 .35],'HorizontalAlignment','center','Interpreter','tex');

    %% ---- WP 星 + WP→tip 連線（tip 落在 WP 半徑 ℓ̂ 上）----
    plot([0 T(1)],[0 T(2)],'-','Color',[.5 .5 .5],'LineWidth',1.4);            % WP→tip = ℓ̂·e
    plot(0,0,'p','MarkerSize',18,'MarkerFaceColor',[1 .85 .1],'MarkerEdgeColor','k','LineWidth',1);
    text(-0.12,-0.02,'WP','FontSize',14,'FontWeight','bold','HorizontalAlignment','right','Interpreter','tex');

    %% ---- tip 方塊 + 沿錐面 λ1 + air-gap λ2（安裝偏移，對照 PDF）----
    plot([T(1) foot(1)],[T(2) foot(2)],'-','Color',[.15 .35 .75],'LineWidth',1.8);   % tip→foot 沿錐面 = λ1·e1
    lm = 0.5*(T(:)+foot);  text(lm(1),lm(2),'  \lambda_1=4.572 mm','FontSize',12.5, ...
         'FontWeight','bold','Color',[.15 .35 .75],'Interpreter','tex');
    plot([foot(1) sx],[foot(2) sz],'-','Color',[.75 .35 .1],'LineWidth',1.8);        % foot→sensor = λ2·n+ (0.41)
    gm = 0.5*([foot(1);foot(2)]+[sx;sz]);  text(gm(1)+0.15,gm(2),'0.41 mm', ...
         'FontSize',12,'FontWeight','bold','Color',[.75 .35 .1],'Interpreter','tex');
    plot(T(1),T(2),'s','MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','k');
    text(T(1),T(2),'  tip','FontSize',13,'FontWeight','bold','Interpreter','tex');

    %% ---- sensor glyph（綠圓盤 + ⊥n+ disc 短線 + 白襯紅 n+ 三角箭頭）----
    DISC = 0.28;  NLEN = 1.0;                                        % disc 半長 / n+ 箭長 [mm]
    tdir = [-nz; nx];                                               % ⊥ n+（disc 平面方向）
    dd = [sx;sz]-DISC*tdir;  de = [sx;sz]+DISC*tdir;
    plot([dd(1) de(1)],[dd(2) de(2)],'-','Color',[.1 .5 .15],'LineWidth',4);          % sensor 感測面
    nt = [sx;sz]+NLEN*[nx;nz];                                       % n+ 箭頭端
    plot([sx nt(1)],[sz nt(2)],'-','Color','w','LineWidth',6);                        % 白襯底
    plot([sx nt(1)],[sz nt(2)],'-','Color',[.92 .15 .15],'LineWidth',3.8);
    plot(nt(1),nt(2),'^','MarkerSize',14,'MarkerFaceColor',[.92 .15 .15], ...
         'MarkerEdgeColor','w','LineWidth',1.0,'MarkerFaceColor',[.92 .15 .15]);
    text(nt(1)+0.10,nt(2),'n_+','Color',[.92 .15 .15],'FontSize',16, ...
         'FontWeight','bold','Interpreter','tex');
    plot(sx,sz,'o','MarkerSize',9,'MarkerFaceColor',[.1 .7 .25],'MarkerEdgeColor','k','LineWidth',1);
    text(sx,sz,'  sensor','FontSize',13,'FontWeight','bold','Color',[.05 .4 .12],'Interpreter','tex');

    %% ---- 軸範圍（涵蓋 WP / tip / sensor / n+ 端）+ 等比 ----
    hold(ax,'off');
    allx = [0 T(1) foot(1) sx nt(1)];  allz = [0 T(2) foot(2) sz nt(2)];
    pad = 0.7;
    xlo=min(allx)-pad; xhi=max(allx)+pad; zlo=min(allz)-pad; zhi=max(allz)+pad;
    daspect(ax,[1 1 1]); xlim(ax,[xlo xhi]); ylim(ax,[zlo zhi]);
    xlabel(ax,'x (mm) (WP frame)'); ylabel(ax,'z (mm) (WP frame)');
    title(ax, sprintf('%s Hall sensor 安裝位置（幾何示意）', POLE),'Interpreter','none');
    ax.Toolbar.Visible = 'off';

    %% ---- [STYLE] 選項①：粗體框圖（figure-style.md）----
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
    yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    set([get(ax,'XLabel') get(ax,'YLabel')],'FontSize',17,'FontWeight','bold');
    set(get(ax,'Title'),'FontSize',15,'FontWeight','bold');

    %% ---- 輸出 ----
    if PREVIEW
        out = fullfile(tempdir, ['sensor_mounting_' POLE '_preview.png']);  res = 150;
    else
        out = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\voltage_base\figures\shared\' ...
               'sensor_mounting_' POLE '.png'];  res = 200;
    end
    exportgraphics(fig, out, 'Resolution', res);
    fprintf('%s sensor (x,z)=(%.3f,%.3f) mm ; n+(x,z)=(%.3f,%.3f) ; saved: %s\n', ...
            POLE, sx, sz, nx, nz, out);
end
