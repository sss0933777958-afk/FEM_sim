% sph_shell  |grad(b.b)| on three nested shells, P1 excited, actuator frame.
%
%   The solid-harmonic model (L = 7, fitted on the six FEM excitations inside
%   r <= 150 um) is a continuous function, so the gradient can be asked for at any
%   point.  Here it is evaluated on three spheres, r = 150 / 100 / 50 um: the shell
%   is swept in (theta, phi), converted to Cartesian, and those coordinates are
%   already the actuator frame the model was fitted in -- no further rotation.
%
%   Colour = | d(b.b)/dx_j | (the three components combined), in mT^2/um, on ONE
%   shared scale for all three shells, so the shells can be compared directly.
%   The shells are drawn as coloured wireframes: an opaque outer surface would hide
%   the inner ones, and a translucent one would blend the colours and break the
%   colour bar.
%
%   Reads only (plot-scripts-pure); sph_shell.mat is written by the fitting run.
%   View 65/25 and the manual box edges follow figures/paper_fig_plot/plot/
%   plot_sphere_lattice_3d.m, the project's 3-D convention.
%
% Output: figures/long2016_hexapole_halfcut/sph_shell.png

HERE = fileparts(mfilename('fullpath'));                      % .../plot/<model>
FMX  = fileparts(fileparts(HERE));                            % .../Force/Maxwell
DAT  = fullfile(FMX,'utils','data');
FIG  = fullfile(FMX,'figures','long2016_hexapole_halfcut','current');

S   = load(fullfile(DAT, 'sph_shell.mat'));
SH  = S.SH;   nS = numel(SH);
allG = cell2mat(arrayfun(@(s) s.Gm(:), SH, 'UniformOutput', false).');
fprintf(['shells r = %s um | |grad(b.b)| %.4f ... %.4f mT^2/um' newline], ...
        mat2str([SH.R]), min(allG), max(allG));

FS = 60;  FSLAB = 36;  CANV = 14.5;  LWBOX = 5.0;
BH = 170;                                                     % box half-width [um]
ST = 4;                                                       % keep every ST-th mesh line
LWS = [2.0 2.6 3.2];                                          % outer -> inner, thicker inside
CUTPH = [110 380];                                            % outer shells: drop the quadrant
                                                              % nearest the camera (az = 65 deg)

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes('Parent',fig,'Position',[0.15 0.11 0.63 0.80]);  hold(ax,'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end

phd = linspace(0, 360, size(SH(1).X,2));                      % phi of the stored grid [deg]
for q = 1:nS                                                  % outermost first
    it = 1:ST:size(SH(q).X,1);   ip = 1:ST:size(SH(q).X,2);
    if q < nS                                                 % cut-away so the inner shells show
        ip = ip(phd(ip) >= CUTPH(1) & phd(ip) <= CUTPH(2));
    end
    surf(ax, SH(q).X(it,ip), SH(q).Y(it,ip), SH(q).Z(it,ip), SH(q).Gm(it,ip), ...
         'FaceColor','none', 'EdgeColor','interp', 'LineWidth',LWS(q));
end
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 16, 'LineWidth', 2.5);  % origin

colormap(ax, turbo);
[CL, ctk] = cticks_(min(allG), max(allG));                     % three inner ticks + both ends
caxis(ax, CL);                                                 %#ok<CAXIS> (clim() shadows the variable)

grid(ax,'off');  box(ax,'off');  daspect(ax,[1 1 1]);
xlim(ax,[-BH BH]);  ylim(ax,[-BH BH]);  zlim(ax,[-BH BH]);
view(ax, 65, 25);                                              % same view as sphere_lattice_3d
set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[0.025 0.025], ...
       'XTickLabelRotation',0,'YTickLabelRotation',0,'ZTickLabelRotation',0);
set(ax,'XTick',-100:100:100,'YTick',-100:100:100,'ZTick',-100:100:100);   % three ticks, no ends
draw_box_(ax, BH, LWBOX);
ax.Clipping = 'off';
xlabel(ax, '$\mathbf{x_a\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, '$\mathbf{y_a\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
zlabel(ax, '$\mathbf{z_a\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);

cb = colorbar(ax, 'Position',[0.815 0.22 0.032 0.58]);
cb.Ticks = ctk(2:end-1);                                       % inner ticks only; the two ends
cb.FontSize = FS;  cb.FontWeight = 'bold';  cb.TickDirection = 'out';   % get numbers without a
                                                               % tick mark, added below
for e = [1 numel(ctk)]                                         % end labels, drawn as text
    yy = cb.Position(2) + (ctk(e)-CL(1))/diff(CL)*cb.Position(4);
    annotation(fig, 'textbox', [cb.Position(1)+cb.Position(3)+0.012, yy-0.035, 0.12, 0.07], ...
               'String', sprintf('%g', ctk(e)), 'EdgeColor','none', 'FontSize',FS, ...
               'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle');
end
cb.Label.Interpreter = 'latex';                                % set BEFORE the string:
cb.Label.String = '$\mathbf{|\nabla(b\cdot b)|\;(mT^2/\mu m)}$';   % tex cannot parse \mathbf
cb.Label.FontSize = FSLAB;
hold(ax,'off');

out = fullfile(FIG, 'sph_shell.png');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 CANV CANV],'PaperSize',[CANV CANV]);
print(fig, out, '-dpng', '-r200');   close(fig);
fprintf(['colour bar %s (limits [%g %g])' newline 'wrote %s' newline], mat2str(ctk), CL, out);

% ---- local helpers ------------------------------------------------------------
function [lim, tk] = cticks_(lo, hi)
% Colour bar: three inner ticks equally spaced, and the two ends carry numbers as
% well (user, 2026-09-15).  The limits are widened to the nearest nice step so all
% five labels are round.
    s = (hi - lo)/4;
    u = 10^floor(log10(s));
    for m = [1 1.5 2 2.5 3 4 5 6 8 10]
        ss = m*u;
        lo2 = floor(lo/ss)*ss;
        if lo2 + 4*ss >= hi - 1e-12, lim = [lo2, lo2 + 4*ss];  tk = lo2 + (0:4)*ss;  return, end
    end
    error('cticks_:none', 'no nice colour scale for [%g %g]', lo, hi);
end

function draw_box_(ax, bh, lw)
% box off + the 9 cube edges drawn by hand (the 3 edges meeting the corner nearest
% the camera are dropped), copied from plot_sphere_lattice_3d.m.
    s = [-bh bh];  [Xc,Yc,Zc] = ndgrid(s,s,s);  C = [Xc(:) Yc(:) Zc(:)];
    E = [];
    for i = 1:8
        for j = i+1:8
            if nnz(abs(C(i,:)-C(j,:)) > 0) == 1, E = [E; i j];  end %#ok<AGROW>
        end
    end
    cp = campos(ax);  [~, near] = min(sum((C - cp).^2, 2));
    E = E(~(E(:,1)==near | E(:,2)==near), :);
    for k = 1:size(E,1)
        p1 = C(E(k,1),:);  p2 = C(E(k,2),:);
        plot3(ax, [p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k-', 'LineWidth', lw);
    end
end
