function plot_P1_circuit_charge(mode, showArrows, Rum)
% PLOT_P1_CIRCUIT_CHARGE(mode, showArrows)
%   P1 side view (xz, y~0): pole outline + WP + fitted equivalent charge
%   (ell_hat*dhat_1, R=50um fit).  If showArrows (default true), overlay the
%   magnetic-circuit quiver (REAL FEM nodes, source convention B=-B_FEM,
%   nearest-y=0 node per grid cell).  If false, charge-point only (no arrows).
%   mode = 'full' | 'zoom'.

    if nargin < 1, mode = 'full'; end
    if nargin < 2, showArrows = true; end
    if nargin < 3, Rum = 50; end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();
    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil1\standard';
    out_dir = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\fix_dir\figures'];
    matf    = sprintf(['G:\\my_workspace\\code\\FEM_sim\\magnetic_sim\\ANSYS\\main\\MATLAB_data\\' ...
                       'long2016_hexapole_halfcut\\charge_fit\\fit_KI_ball\\fit_KI_R%03d.mat'], Rum);
    if ~exist(out_dir,'dir'); mkdir(out_dir); end

    %% geometry (mm)
    R_norm_xy = cnst.R_norm_xy*1e3;  R_norm_z = cnst.R_norm_z*1e3;
    SPH_OFST  = cnst.SPH_OFST*1e3;  POLE_R = cnst.POLE_R*1e3;  POLE_LEN = cnst.POLE_CONE_LEN*1e3;
    Z_CONE = SPH_OFST - R_norm_z;
    x_tip = R_norm_xy;  x_base = x_tip + POLE_LEN;  x_end = (42+53)/2 - 10;

    %% fitted ell + P1 charge (WP frame -> plot frame)
    S = load(matf,'ell');  ellmm = S.ell*1e3;
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    d1   = tip(:,1)/norm(tip(:,1));  cwp = S.ell*d1;
    cx   = cwp(1)*1e3;  cz = cwp(3)*1e3 + SPH_OFST;

    %% view limits
    switch mode
        case 'zoom', xlim_v=[-2,2]; zlim_v=[-14,-11.3]; ytol=0.20; gx=26; gz=26; tag='zoom';
        otherwise,   xlim_v=[-3,15];  zlim_v=[-16,-10];   ytol=0.30; gx=58; gz=26; tag='';
    end

    fig=figure('Position',[60 60 1400 640],'Color','w');
    set(fig,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',12); hold on;

    %% optional magnetic-circuit quiver (real nodes)
    if showArrows
        d  = import_ansys_data(res_dir, 'all', 'coil1');
        xmm=d.x*1e3; ymm=d.y*1e3; zmm=d.z*1e3; bx=-d.bx; by=-d.by; bz=-d.bz;
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

    %% P1 pole outline + WP + charge
    plot([x_tip x_base],[Z_CONE Z_CONE-POLE_R],'k-','LineWidth',2);
    plot([x_base x_end],[Z_CONE-POLE_R Z_CONE-POLE_R],'k-','LineWidth',2);
    plot([x_end x_end],[Z_CONE-POLE_R Z_CONE],'k-','LineWidth',2);
    plot([x_end x_tip],[Z_CONE Z_CONE],'k-','LineWidth',2);
    plot(0,SPH_OFST,'k+','MarkerSize',15,'LineWidth',2.2);
    text(0+0.3,SPH_OFST+0.35,'WP','FontSize',12,'FontWeight','bold');
    plot(cx,cz,'o','MarkerSize',12,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor','k','LineWidth',1.4);
    text(cx+0.35,cz-0.15,'q_{P1}','FontSize',13,'FontWeight','bold','Color',[0.6 0 0.6]);
    hold off;

    axis equal; grid on; xlim(xlim_v); ylim(zlim_v);
    set(gca,'GridAlpha',0.18,'Layer','top'); xlabel('x [mm]'); ylabel('z [mm]');
    if showArrows
        title(sprintf('P1 magnetic circuit + equivalent charge  (l = %.3f mm)', ellmm));
        out = sprintf('P1_circuit_charge_R%d%s.png', Rum, ['_' tag]);  out=strrep(out,'_.','.');
    else
        title(sprintf('P1 equivalent charge position  (l = %.3f mm)', ellmm));
        out = sprintf('P1_charge_only_R%d%s.png', Rum, ['_' tag]);  out=strrep(out,'_.','.');
    end
    exportgraphics(fig, fullfile(out_dir,out), 'Resolution',300);
    fprintf('saved %s\n', out);
end
