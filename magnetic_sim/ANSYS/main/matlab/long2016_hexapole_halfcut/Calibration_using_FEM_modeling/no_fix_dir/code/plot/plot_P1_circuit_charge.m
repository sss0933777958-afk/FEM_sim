function plot_P1_circuit_charge(mode, showArrows, Rum)
% PLOT_P1_CIRCUIT_CHARGE(mode, showArrows, Rum)   [no_fix_dir / 18-param bias]
%   P1 side view (xz, y~0): pole outline + WP + fitted equivalent charge.
%   Charge is the 18-PARAMETER BIAS solution (OFF-AXIS): c = ell_hat*(R'*Pc_18(:,1)),
%   loaded from calib_bias.mat (fit at R_select=150um, actuator frame -> measure frame).
%   This differs from fix_dir (on-axis ell*dhat).  Field arrows + geometry are identical
%   to fix_dir (same coil1 FEM, all-source B=-B_FEM for the lower pole -> tip radiates).
%   mode = 'full' | 'zoom';  showArrows default true.

    if nargin < 1, mode = 'zoom'; end
    if nargin < 2, showArrows = true; end
    if nargin < 3, Rum = 150; end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    cnst = mt_constants();
    res_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\coil1\standard';
    out_dir = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\no_fix_dir\figures'];
    calib_bias = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\' ...
                  'long2016_hexapole_halfcut\charge_fit\calib_bias.mat'];
    if ~exist(out_dir,'dir'); mkdir(out_dir); end

    %% geometry (mm)
    R_norm_xy = cnst.R_norm_xy*1e3;  R_norm_z = cnst.R_norm_z*1e3;
    SPH_OFST  = cnst.SPH_OFST*1e3;  POLE_R = cnst.POLE_R*1e3;  POLE_LEN = cnst.POLE_CONE_LEN*1e3;
    Z_CONE = SPH_OFST - R_norm_z;
    x_tip = R_norm_xy;  x_base = x_tip + POLE_LEN;  x_end = (42+53)/2 - 10;

    %% 18-param bias: ell_hat + OFF-AXIS P1 charge (actuator -> measure -> plot frame)
    CB = load(calib_bias, 'R','Pc_18','ell_hat');  ellmm = CB.ell_hat*1e3;
    c_meas = CB.ell_hat * (CB.R.' * CB.Pc_18(:,1));        % measure/WP frame [m] (k=1 -> P1)
    cx = c_meas(1)*1e3;  cz = c_meas(3)*1e3 + SPH_OFST;  cy = c_meas(2)*1e3;
    fprintf('P1 bias charge (measure): x=%.3f z=%.3f mm, off-plane dy=%.4f mm\n', cx, cz, cy);

    %% view limits
    switch mode
        case 'zoom', xlim_v=[-2,2]; zlim_v=[-14,-11.3]; ytol=0.20; gx=26; gz=26; tag='zoom';
        otherwise,   xlim_v=[-3,15];  zlim_v=[-16,-10];   ytol=0.30; gx=58; gz=26; tag='';
    end

    fig=figure('Position',[60 60 1400 640],'Color','w');
    set(fig,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',12); hold on;

    %% optional magnetic-circuit quiver (real nodes, all-source B=-B_FEM for lower pole)
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
        title(sprintf('P1 magnetic circuit + equivalent charge  (18-param bias, l = %.3f mm)', ellmm));
        out = sprintf('P1_circuit_charge_R%d%s.png', Rum, ['_' tag]);  out=strrep(out,'_.','.');
    else
        title(sprintf('P1 equivalent charge position  (18-param bias, l = %.3f mm)', ellmm));
        out = sprintf('P1_charge_only_R%d%s.png', Rum, ['_' tag]);  out=strrep(out,'_.','.');
    end
    exportgraphics(fig, fullfile(out_dir,out), 'Resolution',300);
    fprintf('saved %s\n', out);
end
