function plot_P2sensor_Braw_P1exc()
% PLOT_P2SENSOR_BRAW_P1EXC
%   P2-sensor 局部磁路箭頭圖。coil1 (P1 excited) RAW FEM（不翻號、不變號），
%   真實 FEM 節點（不內插），y=0 薄 slab 側視，聚焦 P2 cone 面 sensor。
%   目的：看「P1 激發時，對極 P2 的 sensor 處 raw B 方向 vs n+」——
%   為何 B·n+ 近零且與自激同號。
%   箭頭=單位方向（看角度），顏色=|B|(3D, Tesla)。
%
%   PREVIEW 模式：輸出到暫存 png（定案前不落 figures/）。

    %% ---- config（preview 時原地調這裡）----
    SLAB_MM   = 0.30;     % |y| < SLAB_MM 的真實節點切片
    HW_MM     = 2.0;      % zoom 半寬（x,z 各 ±HW_MM，以 P2 sensor 為中心）
    ARROW_MM  = 0.18;     % 單位方向箭頭長度 [mm]
    PREVIEW   = false;    % true=暫存預覽；false=存 figures/（草稿亦放此，方便檢視）
    DPI       = 200;

    %% ---- paths ----
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_no_fix_dir\code\function']);
    cnst = mt_constants();
    rr   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

    %% ---- P2 sensor 幾何（i=2）----
    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);   % WP 框 [m]
    sx = sensor_pos(1,2)*1e3;  sz = sensor_pos(3,2)*1e3;     % P2 sensor (x,z) [mm]
    nx = sensor_n(1,2);        nz = sensor_n(3,2);           % n+ (x,z 分量)
    nn = hypot(nx,nz); nx=nx/nn; nz=nz/nn;                   % 投到 x-z 後正規化
    fprintf('P2 sensor (x,z) = (%.3f, %.3f) mm ; n+(x,z)=(%.3f,%.3f)\n', sx,sz,nx,nz);

    %% ---- 載 coil1 (P1 excited)，套 all-source（P1 下極→negate 使尖端射出/source）----
    %  per charge-model-source-convention：每顆激發極一律當 source。負號 ≡ 反向繞線重跑，
    %  仍是真實 FEM 場（節點原值、不內插）。
    d = import_ansys_data(fullfile(rr,'coil1','standard'),'all','coil1');
    d.bx = -d.bx;  d.by = -d.by;  d.bz = -d.bz;                     % all-source flip (P1 lower)
    xw = d.x*1e3;  yw = d.y*1e3;  zw = (d.z - cnst.SPH_OFST)*1e3;   % WP 框 [mm]

    %% ---- 真實節點切片：y 薄片 + zoom 盒 ----
    m = abs(yw) < SLAB_MM & abs(xw - sx) < HW_MM & abs(zw - sz) < HW_MM;
    fprintf('zoom 盒內 |y|<%.2fmm 真實節點數 = %d\n', SLAB_MM, nnz(m));
    Xx = xw(m); Zz = zw(m);
    Bx = d.bx(m); By = d.by(m); Bz = d.bz(m);
    Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);              % |B| 3D [T]
    Bip  = hypot(Bx, Bz);  Bip(Bip<eps)=eps;         % in-plane mag
    ux = Bx./Bip * ARROW_MM;  uz = Bz./Bip * ARROW_MM;   % 單位方向 × 箭頭長
    % 格點抽樣降密度（每格留一個真實節點，仍是節點原值、非內插）
    CELL = 0.16;
    [~,iu] = unique(round([Xx Zz]/CELL),'rows','stable');
    Xx=Xx(iu); Zz=Zz(iu); ux=ux(iu); uz=uz(iu); Bmag=Bmag(iu);

    %% ---- 圖 ----
    fig = figure('Position',[60 60 1000 950],'Color','w');
    hold on;
    % 依 |B| 上色（turbo），單位方向箭頭
    nb=24; ed=linspace(min(Bmag),max(Bmag),nb+1); cmap=turbo(nb);
    for k=1:nb
        in = Bmag>=ed(k) & Bmag<ed(k+1); if k==nb, in=in|(Bmag>=ed(end)); end
        if any(in)
            quiver(Xx(in),Zz(in),ux(in),uz(in),0,'Color',cmap(k,:), ...
                   'LineWidth',1.3,'MaxHeadSize',2.0);
        end
    end
    colormap(turbo); clim([min(Bmag) max(Bmag)]);
    cb=colorbar; ylabel(cb,'|B| (3D) [Tesla]');

    %% ---- P2 cone 局部外框（WP 框 x-z，非鏡像）----
    inc=cnst.upper_incline; beta=atan2(3.0,15.0);
    T=[cnst.R_norm_xy*cosd(180), cnst.R_norm_z]*1e3;             % P2 tip (x,z) [mm]
    Lsl=sqrt((cnst.POLE_CONE_LEN*1e3)^2+(cnst.POLE_R*1e3)^2);
    ax=[-cos(inc); sin(inc)];                                    % P2 pole axis 在 x-z（朝 -x+z）
    rot=@(v,a)[cos(a)*v(1)-sin(a)*v(2); sin(a)*v(1)+cos(a)*v(2)];
    Ct=T(:)+Lsl*rot(ax,+beta); Cb=T(:)+Lsl*rot(ax,-beta);
    plot([T(1) Ct(1)],[T(2) Ct(2)],'k-','LineWidth',2);
    plot([T(1) Cb(1)],[T(2) Cb(2)],'k-','LineWidth',2);

    %% ---- sensor + n+ + disc ----
    plot(sx,sz,'o','MarkerSize',8,'MarkerFaceColor',[.1 .7 .25],'MarkerEdgeColor','k');
    nt=[sx;sz]+0.85*[nx;nz];                                     % n+ 箭頭 0.85mm
    plot([sx nt(1)],[sz nt(2)],'-','Color','w','LineWidth',5.5);             % 白邊襯底
    plot([sx nt(1)],[sz nt(2)],'-','Color',[.92 .15 .15],'LineWidth',3.6);
    plot(nt(1),nt(2),'^','MarkerSize',13,'MarkerFaceColor',[.92 .15 .15],'MarkerEdgeColor','w','LineWidth',1.0);
    text(nt(1)+0.08,nt(2),'n_+','Color',[.92 .15 .15],'FontSize',15,'FontWeight','bold','Interpreter','tex');
    % disc 平面（⊥ n+）短線
    tdir=[-nz; nx]; dd=[sx;sz]-0.15*tdir; de=[sx;sz]+0.15*tdir;
    plot([dd(1) de(1)],[dd(2) de(2)],'-','Color',[.1 .5 .15],'LineWidth',3);

    hold off; axis equal; grid on;
    xlim([sx-HW_MM sx+HW_MM]); ylim([sz-HW_MM sz+HW_MM]);
    xlabel('x [mm] (WP frame)'); ylabel('z [mm] (WP frame)');
    title('P2 sensor ｜ P1 激發 (all-source)', 'Interpreter','none');

    %% ---- 輸出 ----
    if PREVIEW
        out = fullfile(tempdir,'p2sensor_preview.png');
    else
        out = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_no_fix_dir\figures\P2sensor_Braw_P1exc.png'];
    end
    exportgraphics(fig,out,'Resolution',DPI);
    fprintf('saved: %s\n', out);
end
