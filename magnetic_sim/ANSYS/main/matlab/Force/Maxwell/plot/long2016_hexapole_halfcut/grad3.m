% grad3  Gradient along each actuator axis: RBF vs the solid-harmonic model.
%
%   One figure per axis, each showing the component ALONG that axis:
%       x_a  <- P1 excited,  d(b.b)/dx_a
%       y_a  <- P3 excited,  d(b.b)/dy_a
%       z_a  <- P6 excited,  drawn as -d(b.b)/dz_a  (P6 sits on -z_a, so both models
%              are flipped together; that changes neither the gap nor any sign test)
%
%   Both models are fitted on the SAME data: all 1771 FEM nodes with r <= 150 um,
%   six excitations at once.
%       RBF   rho = 400 um, lam = 4.64e-7 -- the pair with the lowest test NMAE among
%             those whose gradient passes the smoothness test (d3 sign constant).
%       sph   L = 7 (K = 63 coefficients).
%   Units mT^2/um; the force is 0.5*mgB*UF = 22.55 times this, and is NOT plotted here.
%
%   Reads only (plot-scripts-pure); grad3.mat is written by the fitting run.
%
% Output: figures/long2016_hexapole_halfcut/grad_{rbf,sph}_{xa,ya,za}.png  (six figures)

HERE = fileparts(mfilename('fullpath'));                      % .../plot/<model>
FMX  = fileparts(fileparts(HERE));                            % .../Force/Maxwell
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

S = load(fullfile(DAT, 'grad3.mat'));
fprintf(['RBF rho = %g um, lam = %.2e   |   spherical harmonics L = 7 (K = 63)' newline], S.rhoP, S.lamP);

% [MODIFIED 2026-09-15] one figure per MODEL per axis (six in all), not two curves in
%   one frame.  A single series carries no legend (figure-style: name it in the caption)
%   and uses the project's single-series blue.  All six share one y frame so they can be
%   laid side by side.
FS = 60;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;  LW = 6;
CMOD = {[0.85 0.10 0.10], [0.05 0.10 0.95]};                  % RBF red, harmonics blue
XR = [-150 150];   sx = diff(XR)/4;
MODN = {'rbf','sph'};
[YR, YT] = ylim_fit_(min([S.Gr(:); S.Gs(:)]), max([S.Gr(:); S.Gs(:)]), 3);   % shared frame

for a = 1:3
    AN = S.AXN{a};
    fprintf(['%s (P%d): max |RBF - sph| = %.4f mT^2/um (%.2f %% of the peak)' newline], ...
            AN, S.POLE(a), max(abs(S.Gr(:,a)-S.Gs(:,a))), ...
            max(abs(S.Gr(:,a)-S.Gs(:,a)))/max(abs(S.Gs(:,a)))*100);
    for m = 1:2
        if m == 1, y = S.Gr(:,a); else, y = S.Gs(:,a); end

        fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
        ax  = axes('Parent',fig,'Position',[0.1895 0.1367 0.7155 0.7886]); hold(ax,'on');
        try, ax.Toolbar = []; end
        try, ax.Interactions = []; end
        plot(ax, S.sG, y, '-', 'Color',CMOD{m}, 'LineWidth',LW);

        set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
               'TickLength',[.015 .015],'TickDir','out','Box','on','Layer','top', ...
               'XTickLabelRotation',0,'YTickLabelRotation',0);
        xlim(ax, XR);   set(ax,'XTick', XR(1)+(1:3)*sx);
        ylim(ax, YR);   set(ax,'YTick', YT);
        for xv = XR
            text(ax, xv, YR(1)-0.022*diff(YR), sprintf('%g', xv), 'HorizontalAlignment','center', ...
                 'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
        end
        xlabel(ax, ['$\mathbf{' AN '\;(\mu m)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
        if S.SG(a) < 0
            ylabel(ax, ['$\mathbf{-\partial(b\cdot b)/\partial ' AN '\;(mT^2/\mu m)}$'], ...
                   'Interpreter','latex', 'FontSize',FSLAB);
        else
            ylabel(ax, ['$\mathbf{\partial(b\cdot b)/\partial ' AN '\;(mT^2/\mu m)}$'], ...
                   'Interpreter','latex', 'FontSize',FSLAB);
        end
        hold(ax,'off');

        out = fullfile(FIG, ['grad_' MODN{m} '_' strrep(AN,'_','') '.png']);
        set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
        print(fig, out, '-dpng', '-r200');   close(fig);
        fprintf(['  wrote %s' newline], out);
    end
end
fprintf(['shared y frame [%g %g], ticks %s' newline], YR, mat2str(YT));

% ---- local helper ---------------------------------------------------------------
function [lim, tk] = ylim_fit_(lo, hi, N)
% N inner ticks, equally spaced, gaps to both frame edges equal to the tick spacing
% (rules 4 and 5), and the tightest such frame a nice step allows.
    if nargin < 3 || isempty(N), N = 3; end
    half = (N+1)/2;
    need = (hi - lo)/(N+1);
    for kk = floor(log10(need)) : floor(log10(need))+1
        for m = [1 1.5 2 2.5 3 4 5 6 8]
            s = m*10^kk;
            if s < need*(1+1e-9), continue, end
            u = 10^floor(log10(s));
            while abs(s/u - round(s/u)) > 1e-9, u = u/2;  end
            c = ceil((hi - half*s)/u - 1e-9) * u;
            if c - half*s <= lo + 1e-9
                lim = [c - half*s, c + half*s];   tk = c + (-(N-1)/2:(N-1)/2)*s;   return
            end
        end
    end
    error('ylim_fit_:none', 'no nice frame for [%g %g]', lo, hi);
end
