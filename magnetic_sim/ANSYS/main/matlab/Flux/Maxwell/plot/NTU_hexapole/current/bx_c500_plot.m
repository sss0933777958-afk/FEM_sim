% bx_c500_plot.m -- |b_x| on the 17 centre samples (c500: 350..650 um in front of the P1 tip), single pole vs
% pair pole (P1 excitation) as discrete points on one figure; horizontal axis = distance from the P1 tip.
% Loads sp_axis_c500.mat and pp_axis_c500.mat only (written by sp_axis_calc.m / pp_axis_calc.m with SAMP='c500').
% Output (same figure in both folders): figures/NTU_hexapole/{single_pole,pair_pole}/bx_c500.png
% CMP = 'pole' (default, above) | 'range': single pole only, the two sampling ranges (sp_axis.mat: 20..340 um,
% sp_axis_c500.mat: 350..650 um) as two point series -> figures/NTU_hexapole/single_pole/bx_range.png
% SHIFT = true ('range' only): slide the 350..650 series so its first sample sits on the first 20..340 sample
% (shape unchanged; compares the two decay shapes on a common start).
if ~exist('CMP','var'), CMP = 'pole'; end
if ~exist('SHIFT','var'), SHIFT = false; end
clearvars -except CMP SHIFT;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
FIG = fullfile(CALROOT, 'figures', 'NTU_hexapole');
if strcmp(CMP, 'range')
    S1 = load(fullfile(DAT, 'sp_axis.mat'));
    S2 = load(fullfile(DAT, 'sp_axis_c500.mat'));
else
    S1 = load(fullfile(DAT, 'sp_axis_c500.mat'));
    S2 = load(fullfile(DAT, 'pp_axis_c500.mat'));
end
d1 = -S1.x(:)*1e6;  y1 = abs(S1.b(:,1));                  % distance from the tip [um], |b_x| [mT]
d2 = -S2.x(:)*1e6;  y2 = abs(S2.b(:,1));                  % pair: column 1 = P1 excitation
if strcmp(CMP, 'range') && SHIFT, d2 = d2 - (min(d2) - min(d1)); end

% ---- style (figure-style.md) ---------------------------------------------------------
FS = 60;  FSLEG = 45;  FSLAB = 46;  LWBOX = 5.0;  CANV = 14.5;  MS = 22;
CB = [0.05 0.10 0.95];  CR = [0.85 0.10 0.10];
if strcmp(CMP, 'range') && SHIFT, tx = [100 180 260];                   % 20 .. 340: frame starts on the first sample
elseif strcmp(CMP, 'range'), tx = [200 400 600];  else, tx = [425 500 575]; end
sx = tx(2) - tx(1);  xr = [tx(1)-sx, tx(end)+sx];      % 350 .. 650, ends labelled by text
[yr, ty] = ylim_band(min([y1; y2]), max([y1; y2]));

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
AH = 0.70;  if strcmp(CMP, 'range'), AH = 0.80; end          % legend inside for 'range'
ax  = axes(fig, 'Position', [0.20 0.15 0.73 AH]);  hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
h1 = plot(ax, d1, y1, 'o', 'Color',CB, 'MarkerFaceColor',CB, 'MarkerSize',MS, 'LineStyle','none', 'Clipping','off');
h2 = plot(ax, d2, y2, 's', 'Color',CR, 'MarkerFaceColor',CR, 'MarkerSize',MS, 'LineStyle','none', 'Clipping','off');
xlim(ax, xr);  ylim(ax, yr);
set(ax, 'XTick', tx, 'YTick', ty, 'YTickLabel', compose('%g', ty));
box(ax, 'on');  grid(ax, 'off');
set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickDir','in', 'Layer','top');
yoff = yr(1) - 0.025*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%d', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{Distance\;from\;P1\;tip\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{|b_x|\;(mT)}$', 'Interpreter','latex', 'FontSize',FSLAB);
if strcmp(CMP, 'range')
    lg = legend(ax, [h1 h2], {'20-340 \mum', '350-650 \mum'}, 'Location','northeast', 'NumColumns',1, ...
                'FontSize',FSLEG, 'FontWeight','bold');
else
    lg = legend(ax, [h1 h2], {'Single pole', 'Pair pole'}, 'Location','northoutside', 'NumColumns',2, ...
                'FontSize',FSLEG, 'FontWeight','bold');
end
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];

set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
if strcmp(CMP, 'range'), OUTS = {fullfile('single_pole', 'bx_range.png')};
else, OUTS = {fullfile('single_pole', 'bx_c500.png'), fullfile('pair_pole', 'bx_c500.png')}; end
for sub = OUTS
    out = fullfile(FIG, sub{1});
    print(fig, out, '-dpng', '-r200');
    fprintf('wrote %s | y %s ticks %s\n', out, mat2str(yr), mat2str(ty));
end
close(fig);

% ---------------------------------------------------------------------------------------
function [yl, tk] = ylim_band(lo, hi)
% Three equally spaced ticks with end gaps equal to the step; frame edges unlabelled.
    R = hi - lo;  mg = 0.02*R;  e = floor(log10(R/4));  best = inf;  yl = [];  tk = [];
    for ee = e-1:e+1
        for c = [1 1.2 1.5 2 2.5 3 4 5 6 8]
            s  = round(c*10^ee, 12);
            y0 = floor((lo - mg)/(s/2))*(s/2);
            if y0 + 4*s < hi + mg, continue; end
            if 4*s < best, best = 4*s;  yl = round([y0, y0 + 4*s], 12);  tk = round(y0 + (1:3)*s, 12); end
        end
    end
end
