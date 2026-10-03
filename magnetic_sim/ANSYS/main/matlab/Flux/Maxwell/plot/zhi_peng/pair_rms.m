%% pair_rms.m -- calibration RMS versus number of calibration points (one vs two excitations)
%  Loads the two pair calibrations (data/zhi_peng/.mat/calib_pair_P1_R150_xa.mat and ..._P1P2_...), each of which
%  carries the whole ladder (ladder_J, ladder_npts), and plots
%      RMS = sqrt( J / (3*N*M) )    [mT]   per field component, on the N calibration points
%  where J = sum_i ||S_i G - b_i||^2, N = number of points, M = number of excitations.
%  Horizontal axis N = 5 ... 201 (Nr = 2 ... 100 on the pole axis, N = 2*Nr+1).
%  Nine figure rules: tick numbers 60 / legend 45 / box and legend box 5.0, horizontal endpoints numbered by hand
%  and vertical ones not, three equally spaced ticks with end gaps equal to the spacing, square canvas written
%  with print + PaperPosition; tick numbers and legend are tex + bold (latex ignores FontWeight), axis titles latex.
%  Output: figures/zhi_peng/pair_rms.png
clear;
here = fileparts(mfilename('fullpath'));                   % .../plot/zhi_peng
CAL  = fileparts(fileparts(here));                         % .../Flux/Maxwell
MATD = fullfile(CAL, 'data', 'zhi_peng', '.mat');
OUT  = fullfile(CAL, 'figures', 'zhi_peng');
if ~exist(OUT, 'dir'), mkdir(OUT); end

NMIN = 5;   NMAX = 201;
R1 = load(fullfile(MATD, 'calib_pair_P1_R150_xa.mat'));    % one excitation
R2 = load(fullfile(MATD, 'calib_pair_P1P2_R150_xa.mat'));  % two excitations
[N1, y1] = rms_curve(R1, NMIN, NMAX);
[N2, y2] = rms_curve(R2, NMIN, NMAX);
fprintf('one excitation : N = %d..%d, RMS %.4f -> %.4f mT\n', N1(1), N1(end), y1(1), y1(end));
fprintf('two excitations: N = %d..%d, RMS %.4f -> %.4f mT\n', N2(1), N2(end), y2(1), y2(end));

% ---- axes: three ticks, end gaps = spacing ----------------------------------
sx = (NMAX - NMIN) / 4;   xt = NMIN + (1:3)*sx;   xr = [NMIN NMAX];       % 5 ... 201, ticks 54 103 152
ylo = 0.10;   sy = 0.05;  yt = ylo + (1:3)*sy;    yr = [ylo ylo + 4*sy];  % 0.10 ... 0.30, ticks 0.15 0.20 0.25
assert(min([y1; y2]) > yr(1) && max([y1; y2]) < yr(2), 'RMS outside the vertical window');

FS = 60;  FSLEG = 45;  FSLAB = 44;  LWBOX = 5.0;  LW = 5.0;  CANV = 14.5;
BLU = [0.05 0.10 0.95];   RED = [0.85 0.10 0.10];
LM = 3.2;  BM = 2.7;  PS = 10.2;                           % right margin 1.1 in: the endpoint number "201" is half 1.0 in wide
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV],'Visible','off');
ax  = axes(fig,'Units','inches','Position',[LM BM PS PS]);
try, ax.Toolbar = []; end                                  %#ok<TRYNC>
try, ax.Interactions = []; end                             %#ok<TRYNC>
hold(ax,'on');
h1 = plot(ax, N1, y1, '-', 'Color',BLU, 'LineWidth',LW);
h2 = plot(ax, N2, y2, '-', 'Color',RED, 'LineWidth',LW);
xlim(ax, xr);  ylim(ax, yr);
set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'Box','on', ...
        'TickDir','in', 'TickLength',[.018 .018], 'XTick',xt, 'YTick',yt, ...
        'YTickLabel',arrayfun(@(v) sprintf('%.2f',v), yt, 'UniformOutput',false), ...
        'Units','inches', 'Position',[LM BM PS PS]);
xlabel(ax, '$\mathbf{Number\;of\;points}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{RMS\;(mT)}$',          'Interpreter','latex', 'FontSize',FSLAB);
% horizontal endpoints are numbered by hand (they are not ticks); the vertical ones are not numbered
yoff = yr(1) - 0.022*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%d',xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
lg = legend(ax, [h1 h2], {'One excitation','Two excitations'}, 'Interpreter','tex', ...
            'FontSize',FSLEG, 'FontWeight','bold', 'Location','none');
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];
% upper-right corner, inside the empty area above the red curve (its bottom edge stays above the curve)
drawnow;   lg.Units = 'normalized';
% (legend positions are normalised to the FIGURE, so map the axes' upper-right corner, 3 % inside, into it)
lg.Position(1:2) = [(LM + 0.97*PS)/CANV - lg.Position(3), (BM + 0.97*PS)/CANV - lg.Position(4)];
hold(ax,'off');
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
outp = fullfile(OUT, 'pair_rms.png');
print(fig, outp, '-dpng', '-r200');
close(fig);
fprintf('wrote %s\n', outp);

% ============================================================================
function [N, y] = rms_curve(R, NMIN, NMAX)
% RMS per field component on the calibration points, for every ladder level whose N lies in [NMIN, NMAX].
    M  = numel(cellstr(R.EXC));
    N  = R.ladder_npts(:);   J = R.ladder_J(:);
    k  = N >= NMIN & N <= NMAX;
    N  = N(k);   y = sqrt(J(k) ./ (3 * N * M));
end
