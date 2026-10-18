% pp_circuit_plot.m -- NTU pair_pole, P1 excitation: magnetic-circuit arrows on y = 0 (real FEM nodes, no
% interpolation), pole outlines as in pp_sampling_plot.m.  Loads utils/data/NTU_hexapole/pp_slice_y0.mat only
% (written by utils/scripts/NTU_hexapole/pp_slice_y0.m).
% Arrow direction = (b_x, b_z); length ~ (|b|/max)^PW; colour and width binned by |b| (turbo, 28 bins, linear).
% Output: figures/NTU_hexapole/pair_pole/circuit_pair.png (square canvas)
clearvars;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
S   = load(fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole', 'pp_slice_y0.mat'));
OUT = fullfile(CALROOT, 'figures', 'NTU_hexapole', 'pair_pole', 'circuit_pair.png');

% ---- style (figure-style.md) ----------------------------------------------------------
FS = 60;  FSLAB = 46;  LWBOX = 5;
TH = 250;  XL = -1000;                                % plate thickness, left (P2) tip [um]
tx = [-800 -500 -200];  sx = tx(2) - tx(1);  xr = [tx(1)-sx, tx(end)+sx];  % -1100 .. 100 (both tips inside)
tz = [-300 0 300];   sz = tz(2) - tz(1);  zr = [tz(1)-sz, tz(end)+sz];      % -600 .. 600
LEN = 32;  PW = 0.12;                                 % longest arrow [um] (< spacing 34 um); length ~ (|b|/max)^PW
% linear colour scale from the data (rule 14): integer step s with [0, ceil(max/s)*s] in about 4 parts,
% start and end numbered
NS  = [1 2 2.5 4 5 10 20 25 40 50 100 200 250 400 500 1000 2000 2500 4000 5000];
CS  = NS(find(4*NS >= max(S.Bm), 1));  CL = [0, ceil(max(S.Bm)/CS)*CS];  CT = CL(1):CS:CL(2);
nb  = 28;  edges = linspace(CL(1), CL(2), nb+1);  cmap = turbo(nb);  lw = linspace(1.6, 4.0, nb);
LB  = S.Bm;                                           % binning variable (linear |b|)

% square canvas and square axes (1200 x 1200 um); room for the z label + tick numbers on the left and the
% colourbar tick numbers + label on the right; the axes are centred vertically
SQ = 15.0;  LM = 3.0;  CG = 0.45;  CW = 0.40;  RM = 3.4;            % RM: colourbar tick numbers + label
AW = SQ - LM - CG - CW - RM;  AH = AW * diff(zr)/diff(xr);  BM = (SQ - AH)/2 + 0.4;
fig = figure('Color','w','Units','inches','Position',[0.3 0.3 SQ SQ]);
ax  = axes(fig, 'Units','inches', 'Position', [LM BM AW AH]);  hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
fill(ax, [0 xr(2) xr(2) 0],   [-TH/2 -TH/2 TH/2 TH/2], [0.82 0.84 0.88], 'FaceAlpha', 0.35, ...
     'EdgeColor', [0.28 0.30 0.36], 'LineWidth', 3);
fill(ax, [xr(1) XL XL xr(1)], [-TH/2 -TH/2 TH/2 TH/2], [0.82 0.84 0.88], 'FaceAlpha', 0.35, ...
     'EdgeColor', [0.28 0.30 0.36], 'LineWidth', 3);
bxz = hypot(S.Bx, S.Bz);  bxz(bxz == 0) = 1e-12;
scl = LEN * (S.Bm ./ max(S.Bm)).^PW ./ bxz;
U = S.Bx .* scl;  W = S.Bz .* scl;
for q = 1:nb
    if q < nb, m = LB >= edges(q) & LB < edges(q+1); else, m = LB >= edges(q); end
    if any(m)
        quiver(ax, S.X(m) - U(m)/2, S.Z(m) - W(m)/2, U(m), W(m), 0, 'Color', cmap(q,:), ...
               'LineWidth', lw(q), 'MaxHeadSize', 0.6);
    end
end
axis(ax, 'equal');  xlim(ax, xr);  ylim(ax, zr);
set(ax, 'XTick', tx, 'YTick', tz);
box(ax, 'on');
set(ax, 'FontSize', FS, 'FontWeight','bold', 'LineWidth', LWBOX, 'TickDir','in', 'Layer','top', ...
        'Units','inches', 'Position', [LM BM AW AH]);
colormap(ax, turbo);  clim(ax, CL);
yoff = zr(1) - 0.035*diff(zr);
for xv = xr
    text(ax, xv, yoff, sprintf('%d', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize', FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{x\;(\mu m)}$', 'Interpreter','latex', 'FontSize', FSLAB);
ylabel(ax, '$\mathbf{z\;(\mu m)}$', 'Interpreter','latex', 'FontSize', FSLAB);
cb = colorbar(ax, 'Units','inches');
cb.Position = [LM+AW+CG, BM, CW, AH];  ax.Position = [LM BM AW AH];
cb.FontSize = FS;  cb.FontWeight = 'bold';  cb.LineWidth = LWBOX;
cb.Limits = CL;  cb.Ticks = CT;  cb.TickLabels = compose('%d', CT);
cb.Label.Interpreter = 'latex';  cb.Label.String = '$\mathbf{\|b\|\;(mT)}$';  cb.Label.FontSize = FSLAB;
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 SQ SQ], 'PaperSize',[SQ SQ]);
print(fig, OUT, '-dpng', '-r200');
fprintf('saved %s | |b| %.1f .. %.1f mT, colour %s mT, ticks %s\n', OUT, min(S.Bm), max(S.Bm), mat2str(CL), mat2str(CT));
