%% plot_misplace_overlay.m -- sensor 貼不準：P1 單顆 vs 六顆都側偏，region_error(φ) 疊圖
% =========================================================================
% 橫軸：方位偏移 φ [deg]（±90°，繞各自錐軸）；縱軸：WP 球內相對 RMS 場誤差 [%]
% 兩條曲線疊一起：P1 單顆貼歪（藍）vs 六顆都側偏（紅，最壞情況）；灰虛線 = nominal 基準。
% 讀 data/misplace_phi_graded.mat（P1）+ misplace_phi_all_graded.mat（六顆）。風格 ①粗體框圖。
% 輸出實檔 → figures/eighteen_param/misplace_phi_overlay.png（覆蓋迭代）。
% =========================================================================
clear; clc;

VARIANT = 'graded';
here = fileparts(mfilename('fullpath'));
TREE = fileparts(fileparts(here));                       % .../voltage_base
S1 = load(fullfile(TREE,'data',sprintf('misplace_phi_%s.mat',VARIANT)));       % P1
S6 = load(fullfile(TREE,'data',sprintf('misplace_phi_all_%s.mat',VARIANT)));   % 六顆
base = S1.errpct;
fig_dir = fullfile(TREE,'figures','eighteen_param'); if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

col1 = [0.00 0.20 0.60];   % P1 深藍
col6 = [0.85 0.10 0.10];   % 六顆 紅
fig = figure('Color','w','Units','inches','Position',[1 1 9.6 5.2]);
ax = axes(fig); hold(ax,'on');
h6 = plot(ax, S6.phi_deg, S6.eps_region, '-o', 'Color',col6, 'LineWidth',2.5, ...
          'MarkerFaceColor',col6, 'MarkerSize',5);
h1 = plot(ax, S1.phi_deg, S1.eps_region, '-o', 'Color',col1, 'LineWidth',2.5, ...
          'MarkerFaceColor',col1, 'MarkerSize',5);
emax = max([S1.eps_region; S6.eps_region]);
cands = [1 2 5 10 20 25 50];
ystep = cands(find(ceil(emax./cands) <= 12, 1));
ytop  = ceil(emax/ystep)*ystep;
% <2% 容忍窗（以 all-6 最壞情況為準）：縱向虛線框（高度避開頂部圖例）
px = S6.phi_deg(:);  py = S6.eps_region(:);  [~,imn] = min(py);  xmn = px(imn);
xc = local_cross(px, py, 2);
xLc = xc(xc<xmn);  if isempty(xLc), xL = min(px); else, xL = max(xLc); end
xRc = xc(xc>xmn);  if isempty(xRc), xR = max(px); else, xR = min(xRc); end
rectangle(ax,'Position',[xL 0 xR-xL 0.65*ytop],'EdgeColor',[.12 .12 .12],'LineStyle','--','LineWidth',2.2);
text(ax, xmn, 0.72*ytop, sprintf('< 2%%\n%.0f~%.0f deg',xL,xR), 'HorizontalAlignment','center', ...
     'FontSize',12.5,'FontWeight','bold','Color',[.12 .12 .12],'Interpreter','none');
hold(ax,'off');
xlim(ax,[-90 90]); ylim(ax,[0 ytop]);
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box(ax,'on'); grid(ax,'off');
set(ax,'XTick',-90:10:90); set(ax,'YTick',0:ystep:ytop);
set(ax,'XTickLabelRotation',0);
xlabel(ax,'\phi  (deg)','FontWeight','bold');
ylabel(ax,'relative RMS error  (%)','FontWeight','bold');
lg = legend([h1, h6], {'P1 misplaced', 'all 6 misplaced'}, ...
            'Location','north'); lg.FontSize=14; lg.FontWeight='bold';

out = fullfile(fig_dir,'misplace_phi_overlay.png');
exportgraphics(fig, out, 'Resolution', 150);
fprintf('已輸出 %s（<2%% 窗 = [%.1f, %.1f] deg）\n', out, xL, xR);

%% ---- local ----
function xc = local_cross(x, y, L)   % 線性內插找曲線穿越 y=L 的 x
    xc = [];
    for i = 1:numel(x)-1
        if (y(i)-L)*(y(i+1)-L) < 0
            xc(end+1,1) = x(i) + (x(i+1)-x(i))*(L-y(i))/(y(i+1)-y(i)); %#ok<AGROW>
        end
    end
end
