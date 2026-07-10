function plot_P2pole_circuit_2d(PREVIEW)
% PLOT_P2POLE_CIRCUIT_2D  P2 整根磁極 2D 磁路箭頭圖（WP → 支撐座）。
%   graded 密網格、P1 激發（CURR=-1 已 all-source）、真實 FEM 節點、x-z 剖面（y≈0 薄片）。
%   箭頭=單位方向、顏色=|B|(log)。
    if nargin < 1 || isempty(PREVIEW), PREVIEW = false; end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\code\function']);
    cnst = mt_constants();
    rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

    d = import_ansys_data(fullfile(rr,'coil1','graded_p2'),'p2reg_full','coil1');   % 已 all-source（不翻）；_full=Z 擴到 holder 頂
    X = [d.x, d.y, d.z-cnst.SPH_OFST]*1e3;  B = [d.bx, d.by, d.bz];

    % y≈0 薄片 + 視窗（視角拉高，含整個支撐座）
    BX=[-58 3]; BZ=[-3 34];
    m = abs(X(:,2))<1.5 & X(:,1)>BX(1) & X(:,1)<BX(2) & X(:,3)>BZ(1) & X(:,3)<BZ(2);
    Px=X(m,1); Pz=X(m,3); Bx=B(m,1); Bz=B(m,3); Bmag=vecnorm(B(m,:),2,2);

    % 格點抽樣降密度（真實節點）
    CELL=1.0; [~,iu]=unique(round([Px Pz]/CELL),'rows','stable');
    Px=Px(iu); Pz=Pz(iu); Bx=Bx(iu); Bz=Bz(iu); Bmag=Bmag(iu);
    Bip=hypot(Bx,Bz); Bip(Bip<eps)=eps; ARROW=0.95;
    ux=Bx./Bip*ARROW; uz=Bz./Bip*ARROW;

    fig=figure('Position',[60 60 1200 560],'Color','w'); ax=axes(fig); hold(ax,'on');
    Cv=log10(max(Bmag,1e-9)); nb=20; ed=linspace(min(Cv),max(Cv),nb+1); cmap=turbo(nb);
    for k=1:nb
        in=Cv>=ed(k)&Cv<ed(k+1); if k==nb, in=in|(Cv>=ed(end)); end
        if any(in), quiver(Px(in),Pz(in),ux(in),uz(in),0,'Color',cmap(k,:),'LineWidth',1.1,'MaxHeadSize',1.2); end
    end
    colormap(turbo); clim([min(Cv) max(Cv)]); cb=colorbar;
    lo=floor(min(Cv)); hi=ceil(max(Cv)); tk=lo:hi;
    cb.Ticks=tk; cb.TickLabels=arrayfun(@(t)sprintf('10^{%d}',t),tk,'UniformOutput',false); ylabel(cb,'|B| [T]');

    % ---- 幾何輪廓（鋼件外框：P2 錐 + protrusion + yoke）----
    inc=cnst.upper_incline; beta=atan2(cnst.POLE_R, cnst.POLE_CONE_LEN);
    T=[-cnst.R_norm_xy, cnst.R_norm_z]*1e3;                       % P2 尖端 (x,z) [mm]
    adir=[-cos(inc);sin(inc)]; pp=[sin(inc);cos(inc)];           % 極軸 / 垂向（x-z）
    Lc=cnst.POLE_CONE_LEN*1e3/cos(beta); Cend=T(:)+Lc*adir;      % 錐末（到 POLE_R 寬）
    R=cnst.POLE_R*1e3; Cu=Cend+R*pp; Cl=Cend-R*pp; Lb=31;        % 桿身平行延伸到 holder
    Bu=Cu+Lb*adir; Bl=Cl+Lb*adir;
    plot([T(1) Cu(1) Bu(1)],[T(2) Cu(2) Bu(2)],'k-','LineWidth',1.6);   % 上 flank
    plot([T(1) Cl(1) Bl(1)],[T(2) Cl(2) Bl(2)],'k-','LineWidth',1.6);   % 下 flank
    z0=@(zans)(zans-cnst.SPH_OFST)*1e3;                          % ANSYS z → WP z [mm]
    % 支撐座（剖面，由下而上）：yoke 環 → protrusion → holder 塊
    rectangle('Position',[-53,   z0(0),    11, z0(2e-3)-z0(0)],    'EdgeColor','k','LineWidth',1.6);  % yoke 環
    rectangle('Position',[-52.5, z0(2e-3), 10, z0(9e-3)-z0(2e-3)], 'EdgeColor','k','LineWidth',1.6);  % protrusion
    rectangle('Position',[-56,   z0(9e-3), 25, z0(19e-3)-z0(9e-3)],'EdgeColor','k','LineWidth',1.6);  % holder（pole 接於此）
    plot(0,0,'p','MarkerSize',16,'MarkerFaceColor',[1 .84 0],'MarkerEdgeColor','k'); text(0.8,-0.6,'WP','FontWeight','bold','FontSize',13);

    hold(ax,'off'); axis(ax,'equal'); grid(ax,'on'); xlim(BX); ylim(BZ);
    xlabel('x [mm]'); ylabel('z [mm]');
    title('P2 整根磁路（WP → 支撐座）｜ P1 激發 ｜ 真實節點');
    ax.Toolbar.Visible='off';

    if PREVIEW, out=fullfile(tempdir,'p2pole2d.png');
    else, out=['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\figures\P2pole_circuit_2d.png']; end
    exportgraphics(fig,out,'Resolution',200); fprintf('saved: %s\n', out);
end
