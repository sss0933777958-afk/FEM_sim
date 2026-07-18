function plot_flux_box(ZCUT)
% plot_flux_box — 3D 標出「導柱通量積分截面」的位置 + 襯套本體。
%   [MODIFIED 2026-07-17] 原本畫 flange 段 Gauss 盒（①底面②空氣環③側壁）——該方案因 flange 區
%   網格過稀已作廢（flange 厚度只有一層元素、z=16.10 全環只剩 5 個節點、側壁環 104 節點/間距
%   ~1.3mm > 側壁高 1.0mm）。現改為標示使用者定案的**導柱截面 z=14.00**。
%
%   截面：z = ZCUT（預設 14.00）、r <= 2.5 圓盤 @ 導柱軸 (20.5, 0)、n̂ = +z → Φ = Σ B_z·ΔA。
%   位置在襯套 flange 底（15.60）**下方 1.6mm** = 磁通流向（yoke→導柱→板→尖端→WP，B_z<0
%   即往下）的**下游**，故量到的是「襯套分流之後、還留在導柱裡」的通量。
%   襯套本體（flange r7.5 z[15.60,16.60] + tube r3.5 z[16.60,20.80]，bore r2.5 貫穿）一併畫出。
%
%   出 1 檔 → figures/shared/flux_box.png（風格：粗體框選項①）

    if nargin < 1 || isempty(ZCUT), ZCUT = 14.00; end

    % ======================== 幾何（來自 mesh deck MT_Mesh_Graded_Full.txt）========================
    PC   = [20.5, 0];       % 導柱 / 襯套 軸心
    R_PO = 2.5;             % 導柱外徑 = 襯套 bore
    R_TB = 3.5;             % tube 外徑
    R_FL = 7.5;             % flange 外徑
    ZFB  = 15.60;           % flange 底
    ZFT  = 16.60;           % flange 頂 / tube 底
    ZTT  = 20.80;           % tube 頂（接 yoke）
    Z_POST = [9.05, 26.80]; % 導柱 z 範圍
    VIEW = [-38 14];
    % ==============================================================================================

    here   = fileparts(mfilename('fullpath'));
    vbase  = fileparts(fileparts(here));
    figdir = fullfile(vbase,'figures','shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    colPost=[0.55 0.58 0.62]; colSlv=[0.20 0.65 0.35]; colCut=[0.85 0.10 0.10];

    f = figure('Color','w','Position',[100 100 1250 1050]); hold on;

    % ---- 導柱 + 襯套本體（半透明；襯套是環，bore 被導柱填滿 → 畫實心視覺等價）----
    draw_tube(PC(1),PC(2), R_PO, Z_POST(1), Z_POST(2), colPost, 0.13);   % 導柱
    draw_tube(PC(1),PC(2), R_FL, ZFB, ZFT, colSlv, 0.20);                % 襯套 flange
    draw_tube(PC(1),PC(2), R_TB, ZFT, ZTT, colSlv, 0.20);                % 襯套 tube

    % ---- 積分截面 z=ZCUT、r<=R_PO（實心紅盤 + 粗邊）----
    th = linspace(0, 2*pi, 145);
    [TH,RR] = meshgrid(th, linspace(0, R_PO, 20));
    surf(PC(1)+RR.*cos(TH), PC(2)+RR.*sin(TH), ZCUT*ones(size(TH)), ...
         'FaceColor',colCut,'FaceAlpha',0.85,'EdgeColor','none');
    plot3(PC(1)+R_PO*cos(th), PC(2)+R_PO*sin(th), ZCUT*ones(size(th)), '-', ...
          'Color',colCut*0.6,'LineWidth',2.5);

    % ---- n̂ = +z 箭頭（法線）----
    quiver3(PC(1), PC(2), ZCUT, 0, 0, 1.1, 0, 'Color',colCut*0.6, ...
            'LineWidth',2.5,'MaxHeadSize',0.7);
    text(PC(1)+3.4, PC(2)-1.2, ZCUT+0.9, 'n = +z','Color',colCut*0.6, ...
         'FontSize',13,'FontWeight','bold');

    % ---- 標註 ----
    text(PC(1)-4.6, PC(2), ZCUT-0.7, sprintf('cut  z=%.2f,  r<=%.1f', ZCUT, R_PO), ...
         'Color',colCut*0.8,'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
    text(PC(1)+9.0, PC(2), ZFB+0.5, 'sleeve flange','Color',colSlv*0.8, ...
         'FontSize',13,'FontWeight','bold','HorizontalAlignment','center');
    text(PC(1)+6.0, PC(2), ZTT-1.0, 'sleeve tube','Color',colSlv*0.8, ...
         'FontSize',13,'FontWeight','bold','HorizontalAlignment','center');
    text(PC(1)-3.0, PC(2), 21.3, 'post','Color',colPost*0.7, ...      % [MODIFIED] 原 z=26.8 在 zlim(22) 外
         'FontSize',13,'FontWeight','bold','HorizontalAlignment','center');

    % ---- 風格：粗體框選項① ----
    ax = gca;
    xlim([PC(1)-10 PC(1)+10]); ylim([PC(2)-10 PC(2)+10]); zlim([12 22]);
    grid(ax,'off'); box(ax,'on'); daspect(ax,[1 1 1]); view(ax,VIEW(1),VIEW(2));
    set(ax,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    xlabel(ax,'x (mm)','FontWeight','bold'); ylabel(ax,'y (mm)','FontWeight','bold');
    zlabel(ax,'z (mm)','FontWeight','bold');
    xt=get(ax,'XTick'); if numel(xt)>3; set(ax,'XTick',xt(1:2:end)); end
    yt=get(ax,'YTick'); if numel(yt)>3; set(ax,'YTick',yt(1:2:end)); end
    ax.Toolbar.Visible='off';

    out = fullfile(figdir,'flux_box.png');
    exportgraphics(f, out, 'Resolution', 150);  close(f);
    fprintf('saved %s  (cut z=%.2f, A = %.4f mm^2)\n', out, ZCUT, pi*R_PO^2);
end

% ============================================================================
function draw_tube(cx, cy, R, z0, z1, col, alpha)
    [Xc,Yc,Zc] = cylinder(R,48);  Xc=Xc+cx; Yc=Yc+cy;  Zc = Zc*(z1-z0)+z0;
    surf(Xc,Yc,Zc,'FaceColor',col,'FaceAlpha',alpha,'EdgeColor','none');
    th = linspace(0,2*pi,49);
    for zc = [z0 z1]
        patch(cx+R*cos(th), cy+R*sin(th), zc*ones(size(th)), col, ...
              'FaceAlpha',alpha,'EdgeColor','none');
        plot3(cx+R*cos(th), cy+R*sin(th), zc*ones(size(th)), '-','Color',col*0.6,'LineWidth',1.2);
    end
end
