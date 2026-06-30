function plot_Bvector_lowerfilled_P1(HW, out_path, METHOD, ARROWS, STYLE)
% PLOT_BVECTOR_LOWERFILLED_P1
%   Lower-filled hexapole, coil1 (P1) excited, y=0 side view (x-z, APDL frame).
%   Data = data/coil1/lower_filled (coil1=-1 -> P1 SOURCE, B out of tip), no flip.
%   Field sampled by 3D INTERPOLATION onto a y=0 grid (-> shows the WP-front fan).
%   Arrows: ARROWS='uniform' -> all same length (dir only, colour = |B|; matches the
%   image-23 reference, weak back/fan field stays visible); 'var' -> length ~ |B|.
%   Restyled: big bold fonts/ticks, box on, no grid, halved ticks, units (mm)/(T).
%
%   HW (optional): [xlo xhi zlo zhi] mm. Default full [-2 12 -16 -10].
%   out_path (optional): PNG path.
%   METHOD (optional): 'interp' (default) or 'raw' (real nodes).
%   ARROWS (optional): 'uniform' (default, like image 23) or 'var' (length ~ |B|).

    if nargin < 1 || isempty(HW), HW = [-2 12 -16 -10]; end
    if nargin < 2 || isempty(out_path)
        out_path = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
                    'long2016_hexapole_halfcut\field_viz\figures\lowerfilled_P1_field.png'];
    end
    if nargin < 3 || isempty(METHOD), METHOD = 'interp'; end   % 'interp' or 'raw' (real nodes)
    if nargin < 4 || isempty(ARROWS), ARROWS = 'uniform'; end  % 'uniform' (img-23) or 'var'
    if nargin < 5 || isempty(STYLE),  STYLE  = 'new'; end       % 'new'(box/no-grid/(mm)/bold) or 'old'(grid/[mm]/normal)

    FS=16; LWAX=2.0; LW_STEEL=2.2; DPI=150;
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
    cnst = mt_constants();
    beta=atan2(3.0,15.0);
    R_norm_xy=cnst.R_norm_xy*1e3; R_norm_z=cnst.R_norm_z*1e3; SPH_OFST=cnst.SPH_OFST*1e3;
    POLE_R=cnst.POLE_R*1e3; POLE_LEN=cnst.POLE_CONE_LEN*1e3; Z_CONE=SPH_OFST-R_norm_z;

    %% --- load coil1/lower_filled (already SOURCE; no flip) ---
    res_dir = ansys_path('long2016_hexapole_halfcut','data','coil1','lower_filled');
    fprintf('Loading coil1/lower_filled (P1 source)...\n');
    d = import_ansys_data(res_dir,'all','coil1');
    fprintf('  Matched %d nodes, |B| max %.4g T\n', numel(d.x), max(sqrt(d.bx.^2+d.by.^2+d.bz.^2)));

    xlim_v=[HW(1) HW(2)]; zlim_v=[HW(3) HW(4)];

    gx=60;                                              % fixed arrow density (window-independent)
    gz=round(gx*(zlim_v(2)-zlim_v(1))/(xlim_v(2)-xlim_v(1)));
    if strcmpi(METHOD,'raw')
        %% --- y=0 REAL-NODE grid sampling (no interpolation) ---
        xmm=d.x*1e3; ymm=d.y*1e3; zmm=d.z*1e3;
        inwin=xmm>=xlim_v(1)&xmm<=xlim_v(2)&zmm>=zlim_v(1)&zmm<=zlim_v(2)&abs(ymm)<0.8;
        xs=xmm(inwin); ys=ymm(inwin); zs=zmm(inwin);
        bxs=d.bx(inwin); bys=d.by(inwin); bzs=d.bz(inwin);
        xe=linspace(xlim_v(1),xlim_v(2),gx+1); ze=linspace(zlim_v(1),zlim_v(2),gz+1);
        ix=discretize(xs,xe); iz=discretize(zs,ze); ok=~isnan(ix)&~isnan(iz);
        cellid=(iz(ok)-1)*gx+ix(ok);
        xv=xs(ok);zv=zs(ok);bxv=bxs(ok);byv=bys(ok);bzv=bzs(ok);ayv=abs(ys(ok));
        [gc,~]=findgroups(cellid); gi=(1:numel(cellid))';
        pick=accumarray(gc,gi,[],@(r) r(find(ayv(r)==min(ayv(r)),1)));
        Xg=xv(pick);Zg=zv(pick);Bx_g=bxv(pick);By_g=byv(pick);Bz_g=bzv(pick);
        bsum_g=sqrt(Bx_g.^2+By_g.^2+Bz_g.^2);
        fprintf('  raw real-node arrows: %d\n', numel(Xg));
    else
        %% --- y=0 field via 3D INTERPOLATION onto a regular grid (smooth -> WP-front fan) ---
        F_bx=scatteredInterpolant(d.x,d.y,d.z,d.bx,'linear','none');
        F_by=scatteredInterpolant(d.x,d.y,d.z,d.by,'linear','none');
        F_bz=scatteredInterpolant(d.x,d.y,d.z,d.bz,'linear','none');
        xc=linspace(xlim_v(1)+0.1,xlim_v(2)-0.1,gx);
        zc=linspace(zlim_v(1)+0.1,zlim_v(2)-0.1,gz);
        [Xg,Zg]=meshgrid(xc,zc); Xg=Xg(:); Zg=Zg(:);
        Bx_g=F_bx(Xg*1e-3,zeros(size(Xg)),Zg*1e-3);
        By_g=F_by(Xg*1e-3,zeros(size(Xg)),Zg*1e-3);
        Bz_g=F_bz(Xg*1e-3,zeros(size(Xg)),Zg*1e-3);
        keep=~isnan(Bx_g)&~isnan(Bz_g);
        Xg=Xg(keep);Zg=Zg(keep);Bx_g=Bx_g(keep);By_g=By_g(keep);Bz_g=Bz_g(keep);
        bsum_g=sqrt(Bx_g.^2+By_g.^2+Bz_g.^2);
        fprintf('  interpolated arrows: %d\n', numel(Xg));
    end

    %% --- arrows: colour = |B| (cap); length uniform (img-23) or ~ |B| ---
    CAP=0.045;
    mag_c=min(bsum_g,CAP);
    if strcmpi(ARROWS,'var')
        clampf=min(bsum_g,CAP)./max(bsum_g,eps);
        arrow_max=0.55*(xlim_v(2)-xlim_v(1))/14;
        scale=arrow_max/CAP;
        bx_q=Bx_g.*clampf*scale; bz_q=Bz_g.*clampf*scale;
    else                                                % 'uniform' (default, like image 23)
        arrow_len=0.9*(xlim_v(2)-xlim_v(1))/gx;         % fixed length (mm)
        Bip=hypot(Bx_g,Bz_g); Bip(Bip<eps)=eps;         % in-plane magnitude
        bx_q=Bx_g./Bip*arrow_len; bz_q=Bz_g./Bip*arrow_len;
    end

    %% --- plot ---
    fig=figure('Position',[60 60 1180 600],'Color','w');
    ax=axes(fig); hold(ax,'on');
    nb=28; ed=linspace(0,CAP,nb+1); cmap=turbo(nb);
    for k=1:nb
        in=mag_c>=ed(k)&mag_c<ed(k+1); if k==nb, in=in|(mag_c>=ed(end)); end
        if any(in)
            quiver(Xg(in),Zg(in),bx_q(in),bz_q(in),0,'Color',cmap(k,:),'LineWidth',1.1,'MaxHeadSize',0.6);
        end
    end
    % full-cone P1 outline + WP
    T=[R_norm_xy,Z_CONE]; Lsl=sqrt(POLE_LEN^2+POLE_R^2);
    Cu=T+Lsl*[cos(beta),sin(beta)]; Cl=T+Lsl*[cos(beta),-sin(beta)];
    plot([T(1) Cu(1)],[T(2) Cu(2)],'k-','LineWidth',LW_STEEL);
    plot([T(1) Cl(1)],[T(2) Cl(2)],'k-','LineWidth',LW_STEEL);
    plot([Cu(1) Cl(1)],[Cu(2) Cl(2)],'k-','LineWidth',LW_STEEL);
    plot(0,SPH_OFST,'k+','MarkerSize',18,'LineWidth',3);
    OLD = strcmpi(STYLE,'old');
    if OLD, wpw='normal'; else, wpw='bold'; end
    text(0+0.25,SPH_OFST+0.4,'WP','FontSize',FS-(OLD*3),'FontWeight',wpw);
    hold(ax,'off');

    if OLD
        %% --- OLD style: grid on, [mm]/[T], normal weight, full ticks ---
        colormap(turbo); clim([0 CAP]);
        cb=colorbar; cb.Ticks=0:0.005:0.045; ylabel(cb,'|B| [T]');
        axis(ax,'equal'); box(ax,'on'); grid(ax,'on');
        xlim(ax,xlim_v); ylim(ax,zlim_v);
        xlabel(ax,'x [mm]'); ylabel(ax,'z [mm]');
    else
        %% --- NEW style: box, no grid, (mm)/(T), bold, halved ticks ---
        colormap(turbo); clim([0 CAP]);
        cb=colorbar; cb.Ticks=0:0.01:0.04;
        ylabel(cb,'|B| (T)','FontSize',FS,'FontWeight','bold');
        set(cb,'FontSize',FS,'FontWeight','bold','LineWidth',LWAX);
        axis(ax,'equal'); box(ax,'on'); grid(ax,'off');
        xlim(ax,xlim_v); ylim(ax,zlim_v);
        set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWAX,'TickLength',[0.018 0.018]);
        xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));     % halve x ticks
        yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));     % halve z ticks
        xlabel(ax,'x (mm)','FontSize',FS,'FontWeight','bold');
        ylabel(ax,'z (mm)','FontSize',FS,'FontWeight','bold');
    end

    exportgraphics(fig,out_path,'Resolution',DPI);
    fprintf('Saved: %s\n', out_path);
end
