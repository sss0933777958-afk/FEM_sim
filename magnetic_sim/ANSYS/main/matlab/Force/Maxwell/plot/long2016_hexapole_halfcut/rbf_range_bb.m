% rbf_range_bb  U = b.b on the xa axis for a ladder of DATA RANGES, one figure.
%
%   The kernel is held fixed at the settled pair (rho = 120 um, lam = 1e-4) and
%   only the node radius R moves, so the question the figure answers is: does the
%   pair that was tuned at R = 250 um still describe the field once the data
%   reaches out towards the pole tips?
%
%   Reads only.  The .mat comes from one call to the sweep built into rbf_field:
%
%       U    = 'matlab/Force/Maxwell/utils';   addpath(U)
%       kern = struct('rho',120, 'lam',1e-4);
%       [bs, gs, infos] = rbf_field(num2cell(300:50:500), kern);
%       xa = linspace(-150,150,1001).';  Pq = [xa zeros(1001,2)];
%       Ubb(:,k) = sum(bs{k}(Pq).^2, 2);   dU = gs{k}(Pq);
%       Fx(:,k)  = 0.5*mgB*UF*dU(:,1);
%
%   R_norm = 500 um is where the steel starts, so the largest entry reaches the
%   pole tips exactly; per 50 um shell the maximum |B| runs 11.2 mT at the centre,
%   72.1 at 400-450 and 260.9 at 450-500, which is the dynamic range a single
%   global rho has to cover as R grows.
%
% Output: figures/long2016_hexapole_halfcut/rbf_range_bb.png

HERE = fileparts(mfilename('fullpath'));                      % .../plot/<model>
FMX  = fileparts(fileparts(HERE));                            % .../Force/Maxwell
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

S = load(fullfile(DAT, 'rbf_range.mat'));
xa = S.xa;   Ubb = S.Ubb;   RLIST = S.RLIST;   ok = S.okv(:).';
keep = find(ok);                                  % skip any level that failed
assert(~isempty(keep), 'no level in rbf_range.mat completed');

fprintf('rho = %g um, lam = %.0e held fixed\n\n', S.kern.rho, S.kern.lam);
fprintf('%8s %9s %12s %12s %12s\n','R [um]','Np','b.b(0)','b.b(+150)','max spread');
ref = Ubb(:, keep(1));
for k = keep
    fprintf('%8d %9d %12.4f %12.4f %12.4f\n', RLIST(k), S.Np(k), ...
            Ubb(xa==0,k), Ubb(end,k), max(abs(Ubb(:,k) - ref)));
end
fprintf('\nspread across the ladder: max over xa of (max_k - min_k) = %.4f mT^2\n', ...
        max(max(Ubb(:,keep),[],2) - min(Ubb(:,keep),[],2)));

% ================================ figure ======================================
FS = 60;  FSLEG = 45;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 6;
XR = [-150 150];
CM = turbo(numel(keep));
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
hh = gobjects(1,numel(keep));  lbl = cell(1,numel(keep));
for i = 1:numel(keep)
    k = keep(i);
    hh(i)  = plot(ax, xa, Ubb(:,k), '-', 'Color',CM(i,:), 'LineWidth',LW);
    lbl{i} = ['R = ' num2str(RLIST(k)) ' \mum'];
end

set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
       'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
       'XTickLabelRotation',0,'YTickLabelRotation',0);
sx = diff(XR)/4;
xlim(ax, XR);  set(ax,'XTick', XR(1)+(1:3)*sx);
[yr, yt] = ylim_nice_(max(Ubb(:,keep), [], 'all'), 3);
ylim(ax, yr);  set(ax,'YTick', yt);  ytop = yr(2);
for xv = XR
    text(ax, xv, -0.022*ytop, sprintf('%g', xv), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$',     'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{b\cdot b\;(mT^2)}$', 'Interpreter','latex', 'FontSize',FSLAB);
lg = legend(ax, hh, lbl, 'Interpreter','tex', 'Location','northwest', 'NumColumns',1);
lg.FontSize = FSLEG*0.8;  lg.FontWeight = 'bold';
lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
lg.ItemTokenSize = [50 22];
hold(ax,'off');

out = fullfile(FIG, 'rbf_range_bb.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');  close(fig);
fprintf(['wrote %s  (y ticks %s, top %g, fill %.0f%%)' newline], ...
        out, mat2str(yt), ytop, 100*max(Ubb(:,keep),[],'all')/ytop);

function [lim, tk] = ylim_nice_(maxv, N)
    if nargin < 2 || isempty(N), N = 3; end
    smin = 1.08*maxv/(N+1);
    u    = 10^floor(log10(smin));
    s    = ceil(smin/u - 1e-12) * u;
    lim  = [0 (N+1)*s];   tk = (1:N)*s;
end
