function plot_field_xy(sign_vec, tag, X, Y, BX, BY, BZ, scoreval)
% PLOT_FIELD_XY  xy-plane quiver of the superposed field at the WP slice for one
%   source/sink sign combination.  REAL FEM nodes (no interpolation), grid-cell
%   downsampled.  Arrows = in-plane (Bx,By) unit direction, color = |B| (full 3D).
%   Marks WP center (+) and the six pole-tip azimuths with their +/- sign.
%   Inputs (passed from sweep_field_cancellation.m, all over the WP-slice disc):
%     sign_vec 1x6 (+/-1 per pole P1..P6), tag (filename/title), X,Y node xy [m],
%     BX,BY,BZ  Nd x 6 single-coil fields [T], scoreval = that combo's mean|B| over
%     the WP-centre xy slice [T] (shown in the title).
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants
    cnst = mt_constants();
    out_dir = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\field_cancellation\figures'];
    if ~exist(out_dir,'dir'); mkdir(out_dir); end

    s  = sign_vec(:).';
    bx = BX*s.';  by = BY*s.';  bz = BZ*s.';           % superposed field for this combo
    bmag = sqrt(bx.^2 + by.^2 + bz.^2);
    Xmm = X*1e3;  Ymm = Y*1e3;  R = cnst.R_norm*1e3;   % mm

    %% grid-cell downsample (one node per cell, nearest to cell center)
    ng = 26;  xe = linspace(-R,R,ng+1);  ye = linspace(-R,R,ng+1);
    ix = discretize(Xmm,xe);  iy = discretize(Ymm,ye);
    xc = (xe(1:end-1)+xe(2:end))/2;  yc = (ye(1:end-1)+ye(2:end))/2;
    sel = zeros(0,1);
    for a = 1:ng, for b = 1:ng
        c = find(ix==a & iy==b);
        if ~isempty(c)
            [~,jj] = min((Xmm(c)-xc(a)).^2 + (Ymm(c)-yc(b)).^2);
            sel(end+1,1) = c(jj); %#ok<AGROW>
        end
    end, end
    Xs=Xmm(sel); Ys=Ymm(sel); BXs=bx(sel); BYs=by(sel); Bm=bmag(sel);

    fig = figure('Position',[80 80 860 780],'Color','w');
    set(fig,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',12); hold on;
    th = linspace(0,2*pi,200);  plot(R*cos(th),R*sin(th),'-','Color',[.55 .55 .55],'LineWidth',1);  % disc edge

    %% quiver: unit in-plane direction, colored by |B|
    bmax = max(quantile(Bm,0.90), eps);  arrow = 0.85*(2*R)/ng;
    nrm = hypot(BXs,BYs) + eps;  us = BXs./nrm*arrow;  vs = BYs./nrm*arrow;
    nb = 24;  eb = linspace(0,bmax,nb+1);  cm = turbo(nb);
    for q = 1:nb
        if q<nb, m = Bm>=eb(q)&Bm<eb(q+1); else, m = Bm>=eb(q); end
        if any(m), quiver(Xs(m),Ys(m),us(m),vs(m),0,'Color',cm(q,:),'LineWidth',1.0,'MaxHeadSize',0.45); end
    end
    colormap(turbo); clim([0 bmax]); cb = colorbar; ylabel(cb,'|B|  [T]','FontSize',11);

    %% WP center + pole-tip azimuths with sign
    plot(0,0,'k+','MarkerSize',14,'LineWidth',2);
    rt = cnst.R_norm_xy*1e3;
    for i = 1:6
        ti = cnst.pole_angles(i)*pi/180;  px = rt*cos(ti);  py = rt*sin(ti);
        if s(i) > 0, plot(px,py,'ko','MarkerSize',8,'MarkerFaceColor','k');
        else,        plot(px,py,'ko','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1.2); end
        sc = '+'; if s(i) < 0, sc = '-'; end
        text(px*1.16,py*1.16,sprintf('%s%s',cnst.pole_labels{i},sc),'FontSize',11,'FontWeight','bold', ...
             'HorizontalAlignment','center');
    end
    hold off;  axis equal;  xlim([-R*1.3 R*1.3]);  ylim([-R*1.3 R*1.3]);
    grid on; set(gca,'Layer','top','GridAlpha',0.18);  xlabel('x [mm]'); ylabel('y [mm]');
    title(sprintf('WP-slice xy field  (%s)   s=[%+d%+d%+d%+d%+d%+d]   mean|B|_{R50}=%.3f mT', ...
          tag, s, scoreval*1e3));
    out = fullfile(out_dir, sprintf('field_xy_%s.png', tag));
    exportgraphics(fig, out, 'Resolution',300);
    fprintf('saved %s\n', out);
end
