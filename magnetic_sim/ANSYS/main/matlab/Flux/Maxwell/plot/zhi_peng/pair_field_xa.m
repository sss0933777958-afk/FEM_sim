%% pair_field_xa.m -- Zhi-peng pair pole, P1 excitation: magnetic-circuit arrows on y = 0 around the centre
%  (the region of the original 10-01 xa sampling), pole outlines and '+' at the centre.
%  Loads utils/data/zhi_peng/zp_slice_y0.mat.
%  Style as the 10-01 / 10-02 circuit figures: arrow colour and width binned by |b| (turbo, 28 bins).
%  Output: figures/zhi_peng/pair_pole/field_xa.png
clear; close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(here));   % .../Flux/Maxwell
S    = load(fullfile(CALROOT, 'utils', 'data', 'zhi_peng', 'zp_slice_y0.mat'));
OUT  = fullfile(CALROOT, 'figures', 'zhi_peng', 'pair_pole', 'field_xa.png');

FS = 60; FSLAB = 46; LWBOX = 5;
tk = [-300 0 300];  sx = tk(2) - tk(1);  xr = [tk(1)-sx, tk(end)+sx];  yr = xr;
% colour scale with clean, equally spaced ticks (as cbar_clean in pair_circuit.m)
cand = [10 20 50 100 200 500 1000];  bd = inf;
for s = cand
    hi = ceil(max(S.Bm)/s - 1e-9)*s;  n = round(hi/s);
    if n >= 3 && n <= 8 && abs(n-6) < bd, bd = abs(n-6);  CS = s;  CL = [0 hi]; end
end
nb = 28;  edges = linspace(CL(1), CL(2), nb+1);  cmap = turbo(nb);  lw = linspace(1.0, 3.6, nb);

PS = 11;  LM = 2.6;  BM = 2.2;  TM = 0.6;  CG = 0.45;  CW = 0.40;  CRM = 2.9;   % layout [in]
Wf = LM + PS + CG + CW + CRM;  Hf = BM + PS + TM;
fig = figure('Color','w','Units','inches','Position',[0.3 0.3 Wf Hf]);
ax  = axes(fig, 'Units','inches', 'Position', [LM BM PS PS]);  hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
% pole outlines on y = 0 (tongues; config/zhi_peng/Pair_pole, AEDT geometry), centre frame [um]
ZC  = S.C0(3);
POL = {[408.0, xr(2)], [0, 178] - ZC;                 % P1 tongue: x >= 408.0 um, z 0 .. 0.178 mm
       [xr(1), -408.3], [402, 580] - ZC};             % P2 tongue: x <= -408.3 um, z 0.402 .. 0.580 mm
for k = 1:2
    px = POL{k,1};  pz = POL{k,2};
    patch(ax, px([1 2 2 1]), pz([1 1 2 2]), [0.82 0.84 0.88], 'FaceAlpha', 0.30, ...
          'EdgeColor', [0.28 0.30 0.36], 'LineWidth', 4.0);
end
for q = 1:nb
    if q < nb, m = S.Bm >= edges(q) & S.Bm < edges(q+1); else, m = S.Bm >= edges(q); end
    if any(m)
        quiver(ax, S.Xs(m), S.Zs(m), S.Uq(m), S.Wq(m), 0, 'Color', cmap(q,:), 'LineWidth', lw(q), 'MaxHeadSize', 0.35);
    end
end
plot(ax, 0, 0, '+', 'Color', 'k', 'MarkerSize', 44, 'LineWidth', 5.5);       % centre
axis(ax, 'equal');  xlim(ax, xr);  ylim(ax, yr);
set(ax, 'XTick', tk, 'YTick', tk);
box(ax, 'on');
set(ax, 'FontSize', FS, 'FontWeight','bold', 'LineWidth', LWBOX, 'TickDir','in', 'Layer','top', ...
        'Units','inches', 'Position', [LM BM PS PS]);
colormap(ax, turbo);  clim(ax, CL);
yoff = yr(1) - 0.035*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%g', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize', FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{x\;(\mu m)}$', 'Interpreter','latex', 'FontSize', FSLAB);
ylabel(ax, '$\mathbf{z\;(\mu m)}$', 'Interpreter','latex', 'FontSize', FSLAB);
cb = colorbar(ax, 'Units','inches');
cb.Position = [LM+PS+CG, BM, CW, PS];  ax.Position = [LM BM PS PS];
cb.FontSize = FS;  cb.FontWeight = 'bold';  cb.LineWidth = LWBOX;
cb.Limits = CL;  cb.Ticks = CL(1):CS:CL(2);
cb.Label.Interpreter = 'latex';  cb.Label.String = '$\mathbf{\|b\|\;(mT)}$';  cb.Label.FontSize = FSLAB;
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 Wf Hf], 'PaperSize',[Wf Hf]);
print(fig, OUT, '-dpng', '-r150');
fprintf('saved %s\n', OUT);
