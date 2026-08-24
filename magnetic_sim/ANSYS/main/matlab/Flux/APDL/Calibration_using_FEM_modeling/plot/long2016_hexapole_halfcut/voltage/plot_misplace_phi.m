%% plot_misplace_phi.m -- sensor 貼不準：某極繞錐軸方位 φ 偏移 → region_error(φ)
% =========================================================================
% 橫軸：方位偏移 φ [deg]（0~350，繞錐軸一圈）；縱軸：WP 球內相對 RMS 場誤差 [%]
%   ε_region(φ) = 100·‖A·(H_V·V(φ)) − Bstack‖_F / ‖Bstack‖_F（固定 graded 18-param 模型）。
% 資料源：main.m 貼不準區塊存的 data/misplace_phi_<VARIANT>.mat（每顆 sensor 內插 500 點）。
% 灰虛線 = nominal（φ=0）基準；紅點標 φ=10°。風格 ①粗體框圖。
% 輸出實檔 → figures/eighteen_param/misplace_phi.png（覆蓋迭代）。
% =========================================================================
clear; clc;

VARIANT = 'graded';
TAG     = '_all';          % ''=單極（misplace_phi.png）；'_all'=六顆側偏（misplace_phi_all.png）
here = fileparts(mfilename('fullpath'));
TREE = fileparts(fileparts(here));                       % .../voltage_base
S = load(fullfile(TREE,'data',sprintf('misplace_phi%s_%s.mat',TAG,VARIANT)));
phi = S.phi_deg(:); eps_r = S.eps_region(:); base = S.errpct; poles = S.poles;
if numel(poles)==1, plab = sprintf('P%d misplaced',poles); else, plab = sprintf('all %d misplaced',numel(poles)); end
fig_dir = fullfile(TREE,'figures','eighteen_param'); if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

fig = figure('Color','w','Units','inches','Position',[1 1 9.6 5.2]);
ax = axes(fig); hold(ax,'on');
yline(ax, base, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',2);           % nominal 基準
plot(ax, phi, eps_r, '-o', 'Color',[0.00 0.20 0.60], 'LineWidth',2.5, ...
     'MarkerFaceColor',[0.00 0.20 0.60], 'MarkerSize',5);
hold(ax,'off');
cands = [1 2 5 10 20 25 50];                                          % 自適應 y 刻度：目標 ≤12 格
ystep = cands(find(ceil(max(eps_r)./cands) <= 12, 1));                 % P1(≤20)→2%；all-6(≤64)→10%
ytop  = ceil(max(eps_r)/ystep)*ystep;
xlim(ax,[min(phi) max(phi)]); ylim(ax,[0 ytop]);
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box(ax,'on'); grid(ax,'off');
set(ax,'XTick',-90:10:90); set(ax,'YTick',0:ystep:ytop);              % 橫軸每 10°
set(ax,'XTickLabelRotation',0);                                        % 標籤水平、不斜擺
xlabel(ax,'\phi  (deg)','FontWeight','bold');
ylabel(ax,'relative RMS error  (%)','FontWeight','bold');
lg = legend(ax, {'nominal (\phi=0)', plab}, ...
            'Location','north'); lg.FontSize=14; lg.FontWeight='bold';

out = fullfile(fig_dir, sprintf('misplace_phi%s.png',TAG));
exportgraphics(fig, out, 'Resolution', 150);
fprintf('已輸出 %s\n', out);
