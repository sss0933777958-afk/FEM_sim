%% pair_profile.m -- profile cost of l_hat: J / J_min versus l_hat, one and two excitations on one figure
%  Loads the two pair calibrations (data/zhi_peng/.mat/calib_pair_P1_R150_xa.mat and ..._P1P2_...), which carry the
%  profile written by main.m (PAIR_PROF, PAIR_PROF_LO, PAIR_PROF_HI): l_hat swept over the SAME integer range for both
%  cases (500 ... 856 um, from the initial guess upward); at every value l_hat is held fixed and e_hat is fitted again
%  (G eliminated in closed form), giving J(l_hat) on the 11 calibration points. J_min = the cost of the full fit
%  (l_hat free) of the same case, so each curve is normalised by its own minimum and touches 1 at its own l_hat
%  (the dots).
%  Nine figure rules: tick numbers 60 / legend 45 / box and legend box 5.0, horizontal endpoints numbered by hand and
%  vertical ones not, three equally spaced ticks with end gaps equal to the spacing (horizontal numbers are integers:
%  the sweep has integer end points and an integer quarter-width), square canvas written with print + PaperPosition;
%  tick numbers and legend are tex + bold (latex ignores FontWeight), axis titles latex.
%  Output: figures/zhi_peng/pair_profile.png
clear;
here = fileparts(mfilename('fullpath'));                   % .../plot/zhi_peng
CAL  = fileparts(fileparts(here));                         % .../Flux/Maxwell
MATD = fullfile(CAL, 'data', 'zhi_peng', '.mat');
OUT  = fullfile(CAL, 'figures', 'zhi_peng');
if ~exist(OUT, 'dir'), mkdir(OUT); end

cases = {'One excitation',  'calib_pair_P1_R150_xa.mat',   [0.05 0.10 0.95];     % blue
         'Two excitations', 'calib_pair_P1P2_R150_xa.mat', [0.85 0.10 0.10]};    % red
FS = 60;  FSLEG = 45;  FSLAB = 44;  LWBOX = 5.0;  LW = 5.0;  CANV = 14.5;
LM = 3.2;  BM = 2.7;  PS = 10.2;                           % margins / panel side [in]

L = cell(1,2);  Y = cell(1,2);  LF = zeros(1,2);
for q = 1:2
    r = load(fullfile(MATD, cases{q,2}));
    L{q}  = r.profile_l(:) * 1e6;                          % um
    Y{q}  = r.profile_J(:) / r.profile_Jmin;               % J / J_min
    LF(q) = r.l_hat * 1e6;                                 % fitted l_hat (valley), um
    assert(abs(min(Y{q}) - 1) < 1e-6, 'the valley of the sweep is not the full-fit minimum');
    fprintf('%-15s l_hat %.3f um | sweep %.2f .. %.2f um | J/J_min: %.4f at the low end, %.4f at the high end\n', ...
            cases{q,1}, LF(q), L{q}(1), L{q}(end), Y{q}(1), Y{q}(end));
end
assert(abs(L{1}(1) - L{2}(1)) < 1e-6 && abs(L{1}(end) - L{2}(end)) < 1e-6, 'the two sweeps do not share one range');

% horizontal: the shared sweep, integer end points and integer quarter-width, three ticks (end gaps = spacing)
xr = round([L{1}(1) L{1}(end)]);   d = diff(xr) / 4;
assert(d == round(d), 'the sweep does not have an integer quarter-width');
xt = xr(1) + d*(1:3);
% vertical: lower bound 0, three ticks (1..3)*s, upper bound 4*s, smallest nice s that covers both curves
ymax = max([Y{1}; Y{2}]);
smax = [0.5 1 1.5 2 2.5 3 4 5 6 8 10 15 20];   s = smax(find(4*smax >= 1.02*ymax, 1));
yt = (1:3) * s;   yr = [0 4*s];

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV],'Visible','off');
ax  = axes(fig,'Units','inches','Position',[LM BM PS PS]);
try, ax.Toolbar = []; end                                  %#ok<TRYNC>
try, ax.Interactions = []; end                             %#ok<TRYNC>
hold(ax,'on');
h = gobjects(1,2);
for q = 1:2
    h(q) = plot(ax, L{q}, Y{q}, '-', 'Color',cases{q,3}, 'LineWidth',LW);
end
for q = 1:2                                                % the valleys, drawn over both curves
    plot(ax, LF(q), 1, 'o', 'MarkerFaceColor',cases{q,3}, 'MarkerEdgeColor','k', 'MarkerSize',22, 'LineWidth',3);
end
xlim(ax, xr);  ylim(ax, yr);
set(ax, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'Box','on', ...
        'TickDir','in', 'TickLength',[.018 .018], 'XTick',xt, 'YTick',yt, ...
        'XTickLabel',arrayfun(@(v) sprintf('%d',v), xt, 'UniformOutput',false), ...
        'YTickLabel',arrayfun(@(v) sprintf('%g',v), yt, 'UniformOutput',false), ...
        'Units','inches', 'Position',[LM BM PS PS]);
xlabel(ax, '$\mathbf{\hat{\ell}\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{J\,/\,J_{min}}$',       'Interpreter','latex', 'FontSize',FSLAB);
% horizontal endpoints are numbered by hand (they are not ticks); the vertical ones are not numbered
yoff = yr(1) - 0.022*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%d',xv), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
lg = legend(ax, h, cases(:,1).', 'Interpreter','tex', 'FontSize',FSLEG, 'FontWeight','bold', 'Location','none');
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];
% upper-right corner inside the axes, 3 % in (legend positions are normalised to the FIGURE)
drawnow;   lg.Units = 'normalized';
lg.Position(1:2) = [(LM + 0.97*PS)/CANV - lg.Position(3), (BM + 0.97*PS)/CANV - lg.Position(4)];
hold(ax,'off');
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
outp = fullfile(OUT, 'pair_profile.png');
print(fig, outp, '-dpng', '-r200');
close(fig);
fprintf('wrote %s\n', outp);
