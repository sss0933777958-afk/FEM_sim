function plot_P2sensor_tets_3d(PREVIEW)
% PLOT_P2SENSOR_TETS_3D  P2 sensor 取樣範圍被「真實 FEM 四面體」包住的 3D 圖。
% -------------------------------------------------------------------------
%   真實情況：standard 粗網格在 P2 sensor 取樣圓柱(Ø0.3×0.1mm)內 0 節點，
%   168 個取樣點(= 加密版節點位置)落在少數幾顆「大」原始 tet 內。本圖畫出
%   **真實連接性**(mesh_baseline.db CDWRITE 的 sensor_local_*.csv)那幾顆 tet
%   (半透明)+ 168 取樣點(紅)+ sensor 中心 + n+，一眼看出「粗 tet vs 細取樣」。
%   tet 來自真實 FEM 網格，非 Delaunay 估算、非內插。
%   PREVIEW=true → 暫存 png；false → 存 figures/。

    if nargin < 1 || isempty(PREVIEW), PREVIEW = true; end
    DPI = 200;

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
             'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\code\function']);
    cnst = mt_constants();
    [sp, sn] = build_sensor_geometry(cnst);
    rr   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
    MESH = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\mesh\standard\csv';
    SENSOR_R = 0.15e-3; AXIAL_TOL = 0.10e-3;

    % ---- 真實局部網格 → triangulation ----
    N = readmatrix(fullfile(MESH,'sensor_local_nodes.csv'));
    E = readmatrix(fullfile(MESH,'sensor_local_elems.csv'));
    nid = N(:,1); P = N(:,2:4); mxid = max(nid); g2l = zeros(mxid,1); g2l(nid) = 1:numel(nid);
    S = E(:,2:9); tets = zeros(size(S,1),4); kk = 0;
    for r = 1:size(S,1)
        u = unique(S(r,:),'stable'); if numel(u)==4, kk=kk+1; tets(kk,:)=g2l(u); end
    end
    tets = tets(1:kk,:);
    v1=P(tets(:,2),:)-P(tets(:,1),:); v2=P(tets(:,3),:)-P(tets(:,1),:); v3=P(tets(:,4),:)-P(tets(:,1),:);
    vol=dot(v1,cross(v2,v3,2),2); bad=vol<0; tets(bad,[3 4])=tets(bad,[4 3]);
    TR = triangulation(tets, P);

    % ---- P2 取樣點(加密版圓柱內 168 點) → 找包含 tet ----
    dref = import_ansys_data(fullfile(rr,'coil1','sensor_spheres'),'all','coil1');
    c = sp(:,2); ni = sn(:,2);
    r3 = [dref.x,dref.y,dref.z-cnst.SPH_OFST]-c.'; ax=r3*ni; rho=vecnorm(r3-ax*ni.',2,2);
    m = (rho<=SENSOR_R)&(ax>=0)&(ax<=AXIAL_TOL);
    pts = [dref.x(m), dref.y(m), dref.z(m)];                 % ANSYS 框
    ti = pointLocation(TR, pts); tu = unique(ti(~isnan(ti))); % 包含 tet（相異）
    fprintf('P2 取樣點 %d 個，落入 %d 顆真實 FEM tet\n', size(pts,1), numel(tu));

    % ---- 轉 WP 框 [mm]（z 扣 SPH_OFST；sp 已是 WP 框）----
    toWP = @(A) [A(:,1), A(:,2), A(:,3)-cnst.SPH_OFST]*1e3;
    Pw   = toWP(P);  ptsw = toWP(pts);
    scmm = sp(:,2).'*1e3;   nmm = ni.';

    % ---- 3D 圖 ----
    fig = figure('Position',[80 80 1000 880],'Color','w');
    ax = axes(fig,'Position',[0.12 0.10 0.80 0.80]); hold(ax,'on');  % 留白給標題/軸標
    tetramesh(TR.ConnectivityList(tu,:), Pw, 'FaceAlpha',0.16, ...   % 真實 tet 半透明（淡，露出內部點）
              'EdgeColor',[.2 .2 .2],'LineWidth',1.0);
    scatter3(ptsw(:,1),ptsw(:,2),ptsw(:,3), 16,[.85 .1 .1],'filled');  % 168 取樣點
    plot3(scmm(1),scmm(2),scmm(3),'p','MarkerSize',20, ...
          'MarkerFaceColor',[.1 .8 .25],'MarkerEdgeColor','k','LineWidth',1.2);  % sensor 中心
    quiver3(scmm(1),scmm(2),scmm(3), nmm(1),nmm(2),nmm(3), 0.35, ...   % n+
            'Color',[.1 .1 .9],'LineWidth',2.5,'MaxHeadSize',0.6);
    text(scmm(1)+0.35*nmm(1), scmm(2)+0.35*nmm(2), scmm(3)+0.35*nmm(3), ' n_+', ...
         'Color',[.1 .1 .9],'FontSize',13,'FontWeight','bold');

    % 視窗縮到包含 tet 的範圍 + 餘裕
    bb = Pw(unique(TR.ConnectivityList(tu,:)),:); mrg=0.08;
    xlim([min(bb(:,1))-mrg max(bb(:,1))+mrg]);
    ylim([min(bb(:,2))-mrg max(bb(:,2))+mrg]);
    zlim([min(bb(:,3))-mrg max(bb(:,3))+mrg]);
    grid(ax,'on'); box(ax,'on'); view(35,20); daspect(ax,[1 1 1]);
    xlabel('x [mm] (WP frame)'); ylabel('y [mm] (WP frame)'); zlabel('z [mm] (WP frame)');
    title(sprintf('P2 sensor 取樣範圍：%d 顆真實 FEM tet 包住 %d 取樣點（standard 粗網格）', ...
                  numel(tu), size(pts,1)), 'Interpreter','none');
    colormap(ax, lines(numel(tu)));
    ax.Toolbar.Visible = 'off';

    if PREVIEW
        out = fullfile(tempdir,'p2sensor_tets3d_preview.png');
    else
        out = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\figures\P2sensor_tets_3d.png'];
    end
    exportgraphics(fig, out, 'Resolution', DPI);
    fprintf('saved: %s\n', out);
end
