% pp_lfix_plot.m -- NTU pair_pole (P1 excitation) or single_pole, centre samples: |b_x| data (17 points) with
% the 1-D model H = a*I/(x - l)^2 for the fixed l in L_LIST, a fitted for each l.
% CASE = 'pair_pole' (default) | 'single_pole'.
% MODE = 'fit' (default): data + model curves -> bx_lfix.png;  'res': residual y - H at the samples, one
% point series per l (same colours) joined by line segments -> bx_lfix_res.png;
% DROP = l values left out of the residual figure (default none).  LEGLOC = legend corner (default 'northeast').
% YR = fixed vertical range [lo hi] (three ticks at the quarter points; default [] = from the data).  Curve colours follow l (same in both figures); data:
% pair = red squares, single = blue circles (so the l = 1560.5 curve is black, in both figures).
% Loads utils/data/NTU_hexapole/{pp,sp}_lfix_fit.mat only (written by utils/scripts/NTU_hexapole/pp_lfix_fit.m).
% SAMP = 'c500' (default) | 'cyl' (single_pole cylinder model, 100 .. 1700 um; reads *_lfix_fit_cyl.mat).
% Output: figures/NTU_hexapole/<CASE>/bx_lfix[_cyl].png | bx_lfix_res[_cyl].png
if ~exist('CASE','var'), CASE = 'pair_pole'; end
if ~exist('MODE','var'), MODE = 'fit'; end
if ~exist('DROP','var'), DROP = []; end
if ~exist('LEGLOC','var'), LEGLOC = 'northeast'; end
if ~exist('YR','var'), YR = []; end
if ~exist('SAMP','var'), SAMP = 'c500'; end
clearvars -except CASE MODE DROP LEGLOC YR SAMP;  close all;
SFX = '';  if ~strcmp(SAMP, 'c500'), SFX = ['_' SAMP]; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
PFX = struct('pair_pole','pp', 'single_pole','sp');  PFX = PFX.(CASE);
S   = load(fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole', [PFX '_lfix_fit' SFX '.mat']));
FN  = struct('fit',['bx_lfix' SFX '.png'], 'res',['bx_lfix_res' SFX '.png']);
OUT = fullfile(CALROOT, 'figures', 'NTU_hexapole', CASE, FN.(MODE));

% ---- style (figure-style.md) ---------------------------------------------------------
FS = 60;  FSLEG = 45;  FSLAB = 46;  LWBOX = 5.0;  CANV = 14.5;  MS = 22;  LW = 5.0;
CB = [0.05 0.10 0.95];  CR = [0.85 0.10 0.10];
COL = {[0.47 0.67 0.19], [0.49 0.18 0.56], [0 0 0], [0.93 0.62 0.10], [0.5 0.5 0.5]};   % l = 500, 1000, 1560.5, 3000, (Inf)
LS  = {'-', '-', '-', '-', '--'};
% legend: l_hat needs the latex interpreter (tex has no \hat over \ell); \mathbf keeps the digits and unit bold
LEGF = @(s) ['$\mathbf{' s '}$'];
tx = [425 500 575];  if strcmp(SAMP, 'cyl'), tx = [500 1000 1500]; end      % cyl: 0 .. 2000
sx = tx(2) - tx(1);  xr = [tx(1)-sx, tx(end)+sx];      % 350 .. 650, ends labelled by text
KS = 1:numel(S.L_LIST);
if strcmp(MODE, 'res'), KS = KS(~ismember(S.L_LIST, DROP)); end
if strcmp(MODE, 'res'), [yr, ty] = ylim_band(min(min(S.R(:,KS))), max(max(S.R(:,KS))));
else,                  [yr, ty] = ylim_band(min([S.y; S.Hc(:)]), max([S.y; S.Hc(:)])); end
if ~isempty(YR), yr = YR;  ty = YR(1) + (1:3)*diff(YR)/4; end

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig, 'Position', [0.20 0.15 0.73 0.80]);  hold(ax, 'on');     % legend inside -> axes fill the canvas
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
if strcmp(CASE, 'pair_pole'), MK = 's';  CD = CR;  DLAB = 'Pair\;pole';
else,                         MK = 'o';  CD = CB;  DLAB = 'Single\;pole'; end
if strcmp(MODE, 'res'), plot(ax, xr, [0 0], '-', 'Color',[0.6 0.6 0.6], 'LineWidth',2); end   % zero reference
h = gobjects(1, numel(S.L_LIST));  lab = cell(1, numel(S.L_LIST));
for k = KS
    if strcmp(MODE, 'res')
        h(k) = plot(ax, S.d, S.R(:,k), ['-' MK], 'Color',COL{k}, 'MarkerFaceColor',COL{k}, 'MarkerSize',MS, ...
                    'LineWidth',0.6*LW, 'Clipping','off');
    else
        h(k) = plot(ax, S.dc, S.Hc(:,k), LS{k}, 'Color',COL{k}, 'LineWidth',LW);
    end
    if isinf(S.L_LIST(k)), lab{k} = LEGF('\hat{\ell}\;=\;\infty');
    else,                  lab{k} = LEGF(sprintf('\\hat{\\ell}\\;=\\;%g\\;\\mu m', S.L_LIST(k))); end
end
if strcmp(MODE, 'fit'), hd = plot(ax, S.d, S.y, MK, 'Color',CD, 'MarkerFaceColor',CD, 'MarkerSize',MS, 'LineStyle','none', 'Clipping','off');  end
xlim(ax, xr);  ylim(ax, yr);
set(ax, 'XTick', tx, 'YTick', ty, 'YTickLabel', compose('%g', ty));
box(ax, 'on');  grid(ax, 'off');
set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'TickDir','in', 'Layer','top');
yoff = yr(1) - 0.025*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%d', xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
if strcmp(SAMP, 'cyl'), XLAB = '$\mathbf{Distance\;from\;tip\;(\mu m)}$';
else,                   XLAB = '$\mathbf{Distance\;from\;P1\;tip\;(\mu m)}$'; end
xlabel(ax, XLAB, 'Interpreter','latex', 'FontSize',FSLAB);
if strcmp(MODE, 'res')
    ylabel(ax, '$\mathbf{Residual\;(mT)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    hl = h(KS);  ll = lab(KS);
else
    ylabel(ax, '$\mathbf{|b_x|\;(mT)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    hl = [hd h];  ll = [{LEGF(DLAB)}, lab];
end
lg = legend(ax, hl, ll, 'Location',LEGLOC, 'NumColumns',1, ...
            'Interpreter','latex', 'FontSize',FSLEG);
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];

set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
close(fig);
fprintf('wrote %s | y %s ticks %s\n', OUT, mat2str(yr), mat2str(ty));

% ---------------------------------------------------------------------------------------
function [yl, tk] = ylim_band(lo, hi)
% Three equally spaced ticks with end gaps equal to the step; frame edges unlabelled.
    R = hi - lo;  mg = 0.02*R;  e = floor(log10(R/4));  best = inf;  yl = [];  tk = [];
    for ee = e-1:e+1
        for c = [1 1.2 1.5 2 2.5 3 4 5 6 8]
            s  = round(c*10^ee, 12);
            y0 = floor((lo - mg)/s)*s;                         % frame bottom on a multiple of the step
            if y0 + 4*s < hi + mg, continue; end
            if 4*s < best, best = 4*s;  yl = round([y0, y0 + 4*s], 12);  tk = round(y0 + (1:3)*s, 12); end
        end
    end
end
