function plot_pole_circuit_charge_3d(p, showArrows, Rum)
% PLOT_POLE_CIRCUIT_CHARGE_3D(p, showArrows)
%   3D view of paper pole p: cone surface + WP + tip + fitted equivalent charge
%   (ell_hat*dhat_p, R=50um fit).  If showArrows (default true), overlay the
%   magnetic-circuit-in-iron quiver (REAL FEM nodes, source convention B=-B_FEM,
%   3D-grid-decimated, iron nodes only).  If false, charge-point only (no arrows).
%   p = 1 (P1, coil1) or 2 (P2, coil5).

    if nargin < 1, p = 1; end
    if nargin < 2, showArrows = true; end
    if nargin < 3, Rum = 50; end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();
    apdl_to_paper = [1,3,6,5,2,4];
    k = find(apdl_to_paper == p);
    out_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
    matf    = sprintf('G:\\my_workspace\\code\\FEM_sim\\kuo\\data\\long2016_hexapole_halfcut\\fit_KI_ball\\fit_KI_R%03d.mat', Rum);

    %% geometry + fitted charge (WP frame, mm)
    S = load(matf,'ell');  ellmm = S.ell*1e3;
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    tp   = tip(:,p);  dhat = tp/norm(tp);
    ax   = cnst.pole_axis(:,p); ax = ax/norm(ax);
    apex = tp*1e3;  chg = S.ell*dhat*1e3;
    L = cnst.POLE_CONE_LEN*1e3;  R = cnst.POLE_R*1e3;  SPHmm = cnst.SPH_OFST*1e3;
    bh = 3.0;
    fprintf('P%d: ell=%.3f mm, charge=(%.2f,%.2f,%.2f) mm, arrows=%d\n', p, ellmm, chg, showArrows);

    fig = figure('Position',[60 60 980 820],'Color','w'); hold on;

    %% optional magnetic-circuit-in-iron quiver
    if showArrows
        res_dir = sprintf('G:\\my_workspace\\code\\FEM_sim\\kuo\\results\\long2016_hexapole_halfcut\\coil%d', k);
        d  = import_ansys_data(res_dir, 'all', sprintf('coil%d',k));
        xw = d.x*1e3; yw = d.y*1e3; zw = d.z*1e3 - SPHmm;
        bx = -d.bx; by = -d.by; bz = -d.bz;
        air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize', false));
        iron = ~air;  nc = 20;
        m = iron & abs(xw)<bh & abs(yw)<bh & abs(zw)<bh;
        xw=xw(m); yw=yw(m); zw=zw(m); bx=bx(m); by=by(m); bz=bz(m);
        ed = linspace(-bh,bh,nc+1); cc = (ed(1:end-1)+ed(2:end))/2;
        ix=discretize(xw,ed); iy=discretize(yw,ed); iz=discretize(zw,ed);
        cid = ix + (iy-1)*nc + (iz-1)*nc*nc;
        dist = (xw-cc(ix).').^2 + (yw-cc(iy).').^2 + (zw-cc(iz).').^2;
        [~,order] = sort(dist);  [~,ia] = unique(cid(order),'stable');  sel = order(ia);
        X=xw(sel);Y=yw(sel);Z=zw(sel);BX=bx(sel);BY=by(sel);BZ=bz(sel);
        bmag=sqrt(BX.^2+BY.^2+BZ.^2);
        b_clip=quantile(bmag,0.90); arrow=0.9*(2*bh/nc); s=arrow./max(bmag,eps);
        nb=20; eb=linspace(0,b_clip,nb+1); cm=turbo(nb);
        for q=1:nb
            if q<nb, in=bmag>=eb(q)&bmag<eb(q+1); else, in=bmag>=eb(q); end
            if any(in)
                quiver3(X(in),Y(in),Z(in),BX(in).*s(in),BY(in).*s(in),BZ(in).*s(in),0, ...
                        'Color',cm(q,:),'LineWidth',1.0,'MaxHeadSize',0.5);
            end
        end
        colormap(turbo); clim([0 b_clip]); cb=colorbar; ylabel(cb,'|B|  [T]');
        coneAlpha = 0.12;
    else
        coneAlpha = 0.30;
    end

    %% cone surface  (LOWER poles are halfcut -> half cone + flat cut face;
    %                 UPPER poles -> full cone)
    isLower = logical(cnst.pole_is_lower(p));
    tmax = min(1, 2.6/L);  tt = linspace(0,tmax,12);
    if isLower
        v1 = cross(ax,[0;0;1]); v1 = v1/norm(v1);  v2 = [0;0;1];   % keep z <= apex_z
        ph = linspace(pi, 2*pi, 22);                                % lower half only
    else
        vN = null(ax.'); v1 = vN(:,1); v2 = vN(:,2);  ph = linspace(0,2*pi,40);
    end
    [TT,PP]=meshgrid(tt,ph);
    CX = apex(1)+L*TT*ax(1)+(R*TT).*(cos(PP)*v1(1)+sin(PP)*v2(1));
    CY = apex(2)+L*TT*ax(2)+(R*TT).*(cos(PP)*v1(2)+sin(PP)*v2(2));
    CZ = apex(3)+L*TT*ax(3)+(R*TT).*(cos(PP)*v1(3)+sin(PP)*v2(3));
    surf(CX,CY,CZ,'FaceColor',[0.55 0.55 0.55],'FaceAlpha',coneAlpha,'EdgeColor','none');
    if isLower   % flat halfcut face (triangle in the z=apex plane)
        bc = apex + L*tmax*ax;  e1 = bc + R*tmax*v1;  e2 = bc - R*tmax*v1;
        fill3([apex(1) e1(1) e2(1)],[apex(2) e1(2) e2(2)],[apex(3) e1(3) e2(3)], ...
              [0.55 0.55 0.55],'FaceAlpha',coneAlpha,'EdgeColor','none');
    end

    %% markers
    plot3(0,0,0,'k+','MarkerSize',14,'LineWidth',2.2); text(0,0,0.16,'WP','FontSize',12,'FontWeight','bold');
    plot3(apex(1),apex(2),apex(3),'ko','MarkerSize',7,'MarkerFaceColor','k');
    text(apex(1),apex(2),apex(3)+0.16,'tip','FontSize',11);
    plot3(chg(1),chg(2),chg(3),'o','MarkerSize',13,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor','k','LineWidth',1.4);
    text(chg(1)+0.1,chg(2),chg(3)-0.16,sprintf('q_{P%d}',p),'FontSize',13,'FontWeight','bold','Color',[0.6 0 0.6]);

    hold off; grid on; box on; axis equal;
    xlabel('x_{wp} [mm]'); ylabel('y_{wp} [mm]'); zlabel('z_{wp} [mm]');
    xlim([-bh bh]); ylim([-bh bh]); zlim([-bh bh]); view(40,22); set(gca,'FontSize',11);
    if showArrows
        title(sprintf('P%d magnetic circuit in iron + charge  (3D, l = %.3f mm)', p, ellmm));
        out = sprintf('P%d_circuit_charge_R%d_3D.png', p, Rum);
    else
        title(sprintf('P%d equivalent charge position  (3D, l = %.3f mm)', p, ellmm));
        out = sprintf('P%d_charge_only_R%d_3D.png', p, Rum);
    end
    if ~exist(out_dir,'dir'); mkdir(out_dir); end
    exportgraphics(fig, fullfile(out_dir,out), 'Resolution',220);
    fprintf('  saved %s\n', out);
end
