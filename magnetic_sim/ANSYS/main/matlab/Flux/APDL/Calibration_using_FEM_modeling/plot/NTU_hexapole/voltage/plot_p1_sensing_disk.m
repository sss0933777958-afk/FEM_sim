function plot_p1_sensing_disk()
% plot_p1_sensing_disk — NTU_hexapole P1 磁極上的「搜尋範圍」圓面示意圖。
%   P1（0°、下層板 z 8.80–9.05mm）扁平淚滴板 + 一個「水平圓盤」（∥ 板面/xy）：
%     半徑 r = 2.865mm、中心 (13.15, 0, 9.46)（= S2 sensor 平面，板頂上 0.41mm）。
%   此圓 = search range（節點搜尋範圍），非 sensing area。標名 + 中心點；低仰角 3D。
%   純幾何、世界座標 world、單位 mm。
%   風格（粗體框規則，選項①）：box on 標準 3D 框 + daspect 真實比例 +
%     粗體字/LineWidth2 + grid off + tick 減半 + 單位 ()。

    here   = fileparts(mfilename('fullpath'));
    vbase  = fileparts(fileparts(here));
    figdir = fullfile(vbase, 'figures', 'eighteen_param');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath(fullfile(vbase,'code','function'));       % [ADDED] 共用 helper：ntu_pole_profile/ntu_draw_plate

    % ---- P1 幾何（下層板、0°，不旋轉）----
    B2 = ntu_pole_profile();                          % Nx2 磁極邊界（pole-local = world @0°；共用 helper）
    ZB = 8.80;  TH = 0.25;  z0 = ZB;  z1 = ZB+TH;     % 板底/板頂
    colP = [0.30 0.55 0.95];                          % 下層板色

    % ---- search-range 圓盤（水平、∥xy）----
    CX = 13.15;  CY = 0.0;  CZ = 9.46;  RMAX = 2.865; % 中心 + 半徑（使用者拍板）
    colD = [0.90 0.35 0.15];                          % 圓盤色

    f = figure('Color','w','Position',[100 100 1180 940]); hold on;

    % 板
    ntu_draw_plate(B2(:,1), B2(:,2), z0, z1, colP, 0.28);

    % 水平圓盤（半透明 patch + 邊線）
    th = linspace(0, 2*pi, 181);
    dx = CX + RMAX*cos(th);  dy = CY + RMAX*sin(th);  dz = CZ*ones(size(th));
    patch(dx, dy, dz, colD, 'FaceAlpha',0.40, 'EdgeColor',colD*0.7, 'LineWidth',2.5);
    plot3(dx, dy, dz, '-', 'Color',colD*0.6, 'LineWidth',2.5);

    % 中心點
    plot3(CX, CY, CZ, 'o', 'MarkerSize',9, 'MarkerFaceColor',colD*0.8, 'MarkerEdgeColor','k', 'LineWidth',1.2);

    % 清楚的 label（上方白區、大字粗體，離開磁極邊緣）
    text(CX, CY, 10.15, 'search range  (r = 2.865 mm)', ...
        'FontSize',19, 'FontWeight','bold', 'Color',colD*0.7, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom');

    % ---- 樣式（低仰角、真實比例、粗體框）----
    ax = gca;
    xlim(ax,[min(B2(:,1))-1.0, max(B2(:,1))+1.0]);
    ylim(ax,[min(B2(:,2))-1.5, max(B2(:,2))+1.5]);
    zlim(ax,[z0-0.4, 10.4]);
    grid(ax,'off'); box(ax,'on'); daspect(ax,[1 1 1]);
    view(ax,-38, 12);                                 % 低仰角看清圓面
    set(ax,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    xlabel(ax,'x (mm)','FontWeight','bold'); ylabel(ax,'y (mm)','FontWeight','bold'); zlabel(ax,'z (mm)','FontWeight','bold');
    xt=get(ax,'XTick'); if numel(xt)>3; set(ax,'XTick',xt(1:2:end)); end
    yt=get(ax,'YTick'); if numel(yt)>3; set(ax,'YTick',yt(1:2:end)); end
    set(ax,'ZTick', 9:0.5:10);                        % z 軸短→稀疏 tick
    camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
    ax.Toolbar.Visible='off';

    outpng = fullfile(figdir, 'p1_sensing_disk.png');
    exportgraphics(f, outpng, 'Resolution', 150); close(f);
    fprintf('wrote %s\n', outpng);
end

% ============================================================================
% [MODIFIED 2026-07-17] build_pole_profile / arc / draw_plate 已抽到共用 helper
%   code/function/{ntu_pole_profile,ntu_draw_plate}.m（原本本檔與 plot_poles_sensors_3d.m
%   各複製一份、無法跨檔呼叫）。本檔改為呼叫共用版，行為不變。
