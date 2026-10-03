%% rbf_cond.m -- condition number of the RBF Gram matrix against rho, lam = 0.
%  Reads rbf_cond_R<R>.mat from utils/rbf_cond.m. load -> pick -> plot only.
%
%  Markers only, no legend: one series needs no key. The vertical axis is
%  logarithmic, so the figure rules are applied in log space -- an odd number of
%  equally spaced ticks with the end gaps equal to the spacing.
clearvars; clc;

RSEL  = 150;
here  = fileparts(mfilename('fullpath'));
FMX   = fileparts(fileparts(here));
MODEL = 'long2016_hexapole_halfcut';
S     = load(fullfile(FMX,'utils','data',sprintf('rbf_cond_R%d.mat',RSEL)));
OUT   = fullfile(FMX,'figures',MODEL,'current');
if ~exist(OUT,'dir'), mkdir(OUT); end
assert(numel(S.LAM) == 1 && S.LAM(1) == 0, 'this figure is the lam = 0 series');

FS = 60;  FSLAB = 44;  LWBOX = 5.0;  MS = 22;  CANV = 14.5;
BLU = [0 0 1];

% Starts at the first rho actually swept. Span 180 with three ticks and the end
% gaps equal to the spacing forces s = 45, hence 65 / 110 / 155 rather than the
% rounder 50 / 100 / 150 (which would leave gaps of 30 and 50).
XR = [20 200];     XT = [65 110 155];
YR = [1e0 1e16];   YT = [1e4 1e8 1e12];        % 3 ticks, four decades apart

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig,'Units','normalized','Position',[0.215 0.185 0.665 0.640]);
try, ax.Toolbar = []; end                                           %#ok<TRYNC>
try, ax.Interactions = []; end                                      %#ok<TRYNC>
hold(ax,'on');  box(ax,'on');
set(ax,'YScale','log');

plot(ax, S.RHO, S.cond(:,1), 'o', 'Color',BLU, 'MarkerFaceColor',BLU, ...
     'MarkerSize',MS, 'LineStyle','none');

set(ax, 'XLim',XR, 'YLim',YR, 'XTick',XT, 'YTick',YT, ...
        'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
        'TickDir','in', 'TickLength',[0.018 0.018], 'YMinorTick','off');
ax.YAxis.MinorTickValues = [];
xlabel(ax, '$\mathbf{\rho\;(\mu m)}$',  'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{cond(\Phi)}$',     'Interpreter','latex', 'FontSize',FSLAB);

% rule 4: the horizontal axis labels its endpoints by hand (they are not ticks)
ybot = 10^(log10(YR(1)) - 0.022*diff(log10(YR)));
for xv = XR
    text(ax, xv, ybot, sprintf('%g',xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end

set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, fullfile(OUT,sprintf('rbf_cond_R%d.png',RSEL)), '-dpng','-r200');
close(fig);
fprintf('saved %s\n', fullfile(OUT,sprintf('rbf_cond_R%d.png',RSEL)));
