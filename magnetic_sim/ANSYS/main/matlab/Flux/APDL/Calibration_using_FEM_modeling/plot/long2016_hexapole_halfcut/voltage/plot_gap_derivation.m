function plot_gap_derivation(PREVIEW)
%PLOT_GAP_DERIVATION  gap = rf(1-sinβ) 幾何推導示意（鈍尖倒圓 → 直線端點離錐面的垂距）。
%   示意圖：錐面畫成水平、β 放大以利辨識；公式通用（實際 β=11.31°→gap=0.322mm）。
    if nargin < 1 || isempty(PREVIEW), PREVIEW = true; end
    rf = 1;  b = 30*pi/180;                 % 示意：rf=1、β 放大到 30°
    u  = [cos(b); sin(b)];                  % 錐軸方向（tip→base，右上）
    C  = [0; rf];                           % 倒圓圓心（離錐面 rf）
    Tb = [0; 0];                            % 切點（圓與錐面相切）
    T  = C - rf*u;                          % 尖端最前點 = 圓上最前點

    figure('Position',[60 60 1120 820],'Color','w'); ax = axes; hold on;

    % 錐面（水平）
    plot([-2.7 3.6],[0 0],'-','Color','k','LineWidth',2.6);
    text(3.05,-0.14,'錐面 flank','FontSize',15,'FontWeight','bold');
    % 倒圓（半徑 rf）
    th = linspace(0,2*pi,200);
    plot(C(1)+rf*cos(th),C(2)+rf*sin(th),'-','Color',[.2 .5 .85],'LineWidth',1.6);
    % 錐軸（過 C，虛線）+ β
    p0 = C - 2.2*u;  p1 = C + 1.4*u;
    plot([p0(1) p1(1)],[p0(2) p1(2)],'--','Color',[.45 .45 .45],'LineWidth',1.6);
    text(p1(1)+0.05,p1(2)+0.03,'錐軸 axis','FontSize',13,'FontWeight','bold','Color',[.4 .4 .4]);
    xcr = -rf/tan(b);                       % 錐軸與錐面交點 x
    aa = linspace(0,b,20);  R = 0.62;
    plot(xcr+R*cos(aa), R*sin(aa),'-','Color',[.2 .2 .2],'LineWidth',1.5);
    text(xcr+0.78,0.17,'\beta','FontSize',17,'FontWeight','bold');

    % 兩條半徑 rf（C→切點、C→T）
    plot([C(1) Tb(1)],[C(2) Tb(2)],'-','Color',[.85 .3 .1],'LineWidth',1.8);
    text(0.06,0.52,'r_f','FontSize',15,'FontWeight','bold','Color',[.85 .3 .1]);
    plot([C(1) T(1)],[C(2) T(2)],'-','Color',[.85 .3 .1],'LineWidth',1.8);
    text((C(1)+T(1))/2-0.06,(C(2)+T(2))/2+0.11,'r_f','FontSize',15,'FontWeight','bold','Color',[.85 .3 .1]);

    % C、切點、T
    plot(C(1),C(2),'o','MarkerSize',8,'MarkerFaceColor',[.2 .5 .85],'MarkerEdgeColor','k');
    text(C(1)+0.07,C(2)+0.06,'C','FontSize',15,'FontWeight','bold');
    plot(Tb(1),Tb(2),'o','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(T(1),T(2),'s','MarkerSize',12,'MarkerFaceColor','k','MarkerEdgeColor','k');
    text(T(1)-0.02,T(2)+0.15,'尖端前點 T','FontSize',13,'FontWeight','bold');

    % 垂直分解（在 x=T(1)）：上段 rf·sinβ（綠）、下段 gap（紫）
    plot([C(1) T(1)],[C(2) C(2)],':','Color',[.55 .55 .55],'LineWidth',1.2);   % C 的高度線
    plot([T(1) T(1)],[T(2) C(2)],'-','Color',[.1 .55 .3],'LineWidth',3.2);      % rf sinβ
    text(T(1)-1.02,(T(2)+C(2))/2+0.02,'r_f sin\beta','FontSize',13,'FontWeight','bold','Color',[.1 .55 .3]);
    plot([T(1) T(1)],[0 T(2)],'-','Color',[.55 .1 .6],'LineWidth',3.4);         % gap
    text(T(1)-1.72,T(2)/2,'gap = r_f(1 - sin\beta)','FontSize',15,'FontWeight','bold','Color',[.55 .1 .6]);

    % 感測直線（∥錐面，停在高度=gap）
    plot([T(1) 3.2],[T(2) T(2)],'-','Color',[.15 .35 .75],'LineWidth',2.2);
    text(1.35,T(2)+0.14,'感測直線（// 錐面）','FontSize',13,'FontWeight','bold','Color',[.15 .35 .75]);

    daspect([1 1 1]);  xlim([-2.8 3.8]);  ylim([-0.55 2.15]);
    set(ax,'FontSize',14,'FontWeight','bold','LineWidth',2); box on; grid off;
    set(ax,'XTick',[],'YTick',[]);
    title('gap = r_f(1 - sin\beta)   推導示意（\beta 放大；實際 \beta=11.31 度 \rightarrow gap=0.322 mm）', ...
          'FontSize',15,'FontWeight','bold');
    ax.Toolbar.Visible = 'off';

    fdir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\long2016_hexapole_halfcut\Calibration_using_FEM_modeling\voltage_base\figures\shared';
    if PREVIEW, out = fullfile(tempdir,'gap_derivation_preview.png'); res=150;
    else,       out = fullfile(fdir,'gap_derivation.png'); res=200; end
    exportgraphics(gcf, out, 'Resolution', res);
    fprintf('saved %s\n', out);
end
