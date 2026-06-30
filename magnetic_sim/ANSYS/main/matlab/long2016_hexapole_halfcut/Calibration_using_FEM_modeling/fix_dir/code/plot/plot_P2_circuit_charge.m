function plot_P2_circuit_charge(showArrows, Rum)
% PLOT_P2_CIRCUIT_CHARGE(showArrows)
%   P2 side view (xz, y~0): full cone outline + WP + fitted equivalent charge
%   (ell_hat*dhat_2, R=50um fit).  P2 is an UPPER pole (NOT halfcut), at -x/+z,
%   excited by coil5.  showArrows (default true): overlay magnetic-circuit quiver
%   (real nodes, source convention B=-B_FEM, nearest-y=0 node per grid cell).

    if nargin < 1, showArrows = true; end
    if nargin < 2, Rum = 50; end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();
    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil5\standard';   % P2
    out_dir = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\fix_dir\figures'];
    matf    = sprintf(['G:\\my_workspace\\code\\FEM_sim\\magnetic_sim\\ANSYS\\main\\matlab\\' ...
                       'long2016_hexapole_halfcut\\Calibration_using_FEM_modeling\\fix_dir\\data\\fit_KI_R%03d.mat'], Rum);
    if ~exist(out_dir,'dir'); mkdir(out_dir); end

    SPH = cnst.SPH_OFST*1e3;  POLE_R = cnst.POLE_R*1e3;  POLE_LEN = cnst.POLE_CONE_LEN*1e3;

    %% P2 geometry + fitted charge (WP frame -> plot frame, mm)
    S = load(matf,'ell');  ellmm = S.ell*1e3;
    tip = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    tp2 = tip(:,2);  d2 = tp2/norm(tp2);  cwp = S.ell*d2;
    Tx = tp2(1)*1e3;  Tz = tp2(3)*1e3 + SPH;            % P2 tip (plot coords)
    cx = cwp(1)*1e3;  cz = cwp(3)*1e3 + SPH;            % charge
    ax2 = cnst.pole_axis(:,2);  ux = ax2(1); uz = ax2(3);  un = hypot(ux,uz); ux=ux/un; uz=uz/un;

    %% view limits (zoom on P2 tip + WP)
    xlim_v=[-4.5,1.5]; zlim_v=[-13.2,-9.8]; ytol=0.20; gx=46; gz=30;

    fig=figure('Position',[60 60 1400 760],'Color','w');
    set(fig,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',12); hold on;

    %% optional quiver (real nodes)
    if showArrows
        d = import_ansys_data(res_dir,'all','coil5');
        % all-source: P2 is an UPPER pole (raw FEM already source), so KEEP sign
        % (s=+1 per [-1,+1,-1,+1,+1,-1]); arrows radiate OUT of the P2 tip.
        xmm=d.x*1e3; ymm=d.y*1e3; zmm=d.z*1e3; bx=d.bx; by=d.by; bz=d.bz;
        m = abs(ymm)<ytol & xmm>=xlim_v(1)&xmm<=xlim_v(2)&zmm>=zlim_v(1)&zmm<=zlim_v(2);
        xm=xmm(m);ym=ymm(m);zm=zmm(m);bxm=bx(m);bym=by(m);bzm=bz(m);
        xe=linspace(xlim_v(1),xlim_v(2),gx+1); ze=linspace(zlim_v(1),zlim_v(2),gz+1);
        ix=discretize(xm,xe); iz=discretize(zm,ze); sel=zeros(0,1);
        for a=1:gx, for b=1:gz
            c=find(ix==a&iz==b); if ~isempty(c), [~,j]=min(abs(ym(c))); sel(end+1,1)=c(j); end
        end, end
        X=xm(sel);Z=zm(sel);BX=bxm(sel);BZ=bzm(sel);
        bmag=sqrt(bxm(sel).^2+bym(sel).^2+bzm(sel).^2);
        b_clip=quantile(bmag,0.90); arrow=0.85*(xlim_v(2)-xlim_v(1))/gx; s=arrow./max(bmag,eps);
        nb=24; eb=linspace(0,b_clip,nb+1); cm=turbo(nb); lwr=[0.55 1.9];
        for q=1:nb
            if q<nb, in=bmag>=eb(q)&bmag<eb(q+1); else, in=bmag>=eb(q); end
            if any(in)
                lw=lwr(1)+(q-1)/(nb-1)*(lwr(2)-lwr(1));
                quiver(X(in),Z(in),BX(in).*s(in),BZ(in).*s(in),0,'Color',cm(q,:),'LineWidth',lw,'MaxHeadSize',0.6);
            end
        end
        colormap(turbo); clim([0 b_clip]); cb=colorbar; ylabel(cb,'|B|  [T]','FontSize',11);
    end

    %% P2 full cone outline (two generators about the axis + base line)
    ha = atan(POLE_R/POLE_LEN); ca=cos(ha); sa=sin(ha);
    g1 = [ux*ca-uz*sa,  ux*sa+uz*ca];     % axis rotated +ha
    g2 = [ux*ca+uz*sa, -ux*sa+uz*ca];     % axis rotated -ha
    slant = POLE_LEN/ca;
    B1 = [Tx Tz] + slant*g1;  B2 = [Tx Tz] + slant*g2;
    plot([Tx B1(1)],[Tz B1(2)],'k-','LineWidth',2);
    plot([Tx B2(1)],[Tz B2(2)],'k-','LineWidth',2);
    plot([B1(1) B2(1)],[B1(2) B2(2)],'k-','LineWidth',2);

    %% WP + charge
    plot(0,SPH,'k+','MarkerSize',15,'LineWidth',2.2);
    text(0+0.15,SPH-0.18,'WP','FontSize',12,'FontWeight','bold');
    plot(cx,cz,'o','MarkerSize',12,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor','k','LineWidth',1.4);
    text(cx+0.18,cz+0.18,'q_{P2}','FontSize',13,'FontWeight','bold','Color',[0.6 0 0.6]);
    hold off;

    axis equal; grid on; xlim(xlim_v); ylim(zlim_v);
    set(gca,'GridAlpha',0.18,'Layer','top'); xlabel('x [mm]'); ylabel('z [mm]');
    if showArrows
        title(sprintf('P2 magnetic circuit + equivalent charge  (l = %.3f mm)', ellmm));
        out=sprintf('P2_circuit_charge_R%d_zoom.png', Rum);
    else
        title(sprintf('P2 equivalent charge position  (l = %.3f mm)', ellmm));
        out=sprintf('P2_charge_only_R%d_zoom.png', Rum);
    end
    exportgraphics(fig, fullfile(out_dir,out), 'Resolution',300);
    fprintf('P2: charge=(%.2f,%.2f) tip=(%.2f,%.2f)  saved %s\n', cx,cz,Tx,Tz, out);
end
