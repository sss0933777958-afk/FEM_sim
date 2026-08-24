%% plot_slide_soff.m -- sensor 縱向滑動：沿錐面 soff 位移（φ=0，固定 nominal 模型）→ region_error(soff)
% =========================================================================
% 橫軸：沿錐母線位置 soff [mm]（2.05→8，谷底 soff=4.572 nominal）；縱軸：WP 球內相對 RMS 場誤差 [%]
% 兩條疊一起：P1 單顆縱滑（藍）vs 六顆都縱滑（紅）。虛線框 = all-6（最壞情況）region_error < 2% 的容忍窗。
% 讀 data/slide_soff_graded.mat（main.m SLIDE_SOFF 區塊產）。風格 ①粗體框圖。
% 輸出實檔 → figures/eighteen_param/slide_soff_overlay.png（覆蓋迭代）。
% =========================================================================
clear; clc;

VARIANT = 'graded';
here = fileparts(mfilename('fullpath'));  TREE = fileparts(fileparts(here));   % .../voltage_base
S = load(fullfile(TREE,'data',sprintf('slide_soff_%s.mat',VARIANT)));
x = S.soff_mm(:);  eP1 = S.eps_P1(:);  eAll = S.eps_all(:);
keep = x <= 8+1e-9;  x = x(keep);  eP1 = eP1(keep);  eAll = eAll(keep);       % 只畫到 8mm
fig_dir = fullfile(TREE,'figures','eighteen_param'); if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

col1 = [0.00 0.20 0.60];   % P1 藍
col6 = [0.85 0.10 0.10];   % 六顆 紅
fig = figure('Color','w','Units','inches','Position',[1 1 9 5.4]);
ax = axes(fig); hold(ax,'on');
h6 = plot(ax, x, eAll, '-o', 'Color',col6, 'LineWidth',2.5, 'MarkerFaceColor',col6, 'MarkerSize',5);
h1 = plot(ax, x, eP1,  '-o', 'Color',col1, 'LineWidth',2.5, 'MarkerFaceColor',col1, 'MarkerSize',5);
emax = max([eP1; eAll]);
cands = [1 2 5 10 20 25 50];
ystep = cands(find(ceil(emax./cands) <= 12, 1));
ytop  = ceil(emax/ystep)*ystep;
% <2% 容忍窗（以 all-6 最壞情況為準）：縱向虛線框（高度避開頂部圖例）
[~,imn] = min(eAll);  xmn = x(imn);
xc = local_cross(x, eAll, 2);
xLc = xc(xc<xmn);  if isempty(xLc), xL = min(x); else, xL = max(xLc); end
xRc = xc(xc>xmn);  if isempty(xRc), xR = max(x); else, xR = min(xRc); end
rectangle(ax,'Position',[xL 0 xR-xL 0.65*ytop],'EdgeColor',[.12 .12 .12],'LineStyle','--','LineWidth',2.2);
text(ax, xmn, 0.72*ytop, sprintf('< 2%%\n%.2f-%.2f mm',xL,xR), 'HorizontalAlignment','center', ...
     'FontSize',12.5,'FontWeight','bold','Color',[.12 .12 .12],'Interpreter','none');
hold(ax,'off');
xlim(ax,[min(x) 8]); ylim(ax,[0 ytop]);
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box(ax,'on'); grid(ax,'off');
set(ax,'XTick',2:1:8); set(ax,'YTick',0:ystep:ytop);
xlabel(ax,'cone soff  (mm)','FontWeight','bold');
ylabel(ax,'relative RMS error  (%)','FontWeight','bold');
lg = legend([h1 h6], {'P1 slid', 'all 6 slid'}, 'Location','north');
lg.FontSize=14; lg.FontWeight='bold';

out = fullfile(fig_dir,'slide_soff_overlay.png');
exportgraphics(fig, out, 'Resolution', 150);
fprintf('已輸出 %s（<2%% 窗 = [%.2f, %.2f] mm）\n', out, xL, xR);

%% ---- local ----
function xc = local_cross(x, y, L)   % 線性內插找曲線穿越 y=L 的 x
    xc = [];
    for i = 1:numel(x)-1
        if (y(i)-L)*(y(i+1)-L) < 0
            xc(end+1,1) = x(i) + (x(i+1)-x(i))*(L-y(i))/(y(i+1)-y(i)); %#ok<AGROW>
        end
    end
end
