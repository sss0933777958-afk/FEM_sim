function plot_sensor_quadrant_B_3d(REGION)
% plot_sensor_quadrant_B_3d — P1 sensor 圓柱內的 B（大小 + 方向）3D 圖。（檔名 legacy：已不限象限）
%   coil1 激發。取 sensor 圓柱（R0.15×H0.10 @ (13.15,-0.99,9.46)）內的內插取樣點，
%   畫各點的**完整 3D B 向量**（長度 ∝ |B|、顏色 = |B| 色階）。**不做 n+ 投影**（使用者指定）。
%   含襯套 / 無襯套各一張，**共用色階與長度尺度**。
%
%   [ADDED] REGION：取哪些點
%     'all'（預設）→ **整個圓柱的 1000 點**（抽稀到 NSHOW 支箭頭）→ sensor_area_B_3d_*.png
%     'quadrants'  → 只取 +x+y 與 -x-y 兩個對角象限（舊行為）    → sensor_quadrant_B_3d_*.png
%
%   ⚠ 此圖的場為**重心內插值**（非 FEM 節點原值）—— 使用者指定要看 build_V_matrix 內插出來的場。
%     內插法與 build_V_matrix 完全相同：中心 0.30mm 內真實節點建 Delaunay → 圓柱內面積均勻撒點
%     → pointLocation + cartesianToBarycentric → B(p)=Σλ_i·B_i。
%
%   出 2 檔 → figures/shared/sensor_quadrant_B_3d_{sleeve,nosleeve}.png

    % ======================== 參數 ========================
    NSAMP  = 1000;          % 圓柱內取樣點數（同 build_V_matrix 的 n_uniform）
    NSHOW  = 1000;          % [MODIFIED] 畫出全部 1000 點（使用者指定）
    LMAX_F = 0.42;          % 最大 |B| 對應的箭頭長度 = LMAX_F × sensor 半徑 [mm]
    VIEW   = [-32 12];      % [MODIFIED] 仰角壓低 → 看得出箭頭在圓柱「裡面」
    % =====================================================

    here   = fileparts(mfilename('fullpath'));
    vbase  = fileparts(fileparts(here));
    figdir = fullfile(vbase,'figures','shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath(fullfile(vbase,'code','function'));
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live 樹。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    % ⚠ addpath 預設「前置」→ 最後加的贏。NTU mt_constants 必須排在 backup 之後才蓋得過 long2016 同名檔。
    addpath(fullfile(fileparts(vbase),'current_base','code','main_function'));   % NTU mt_constants
    droot = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';

    cnst = model_config('NTU_hexapole');
    c  = cnst.sensor_pos(:,1)*1e3;          % P1 sensor 中心 [mm]
    ni = cnst.sensor_n(:,1);                % n+ = +z
    SR = cnst.SENSOR_R*1e3;  SH = cnst.SENSOR_H*1e3;

    if nargin < 1 || isempty(REGION), REGION = 'all'; end     % [ADDED] 預設整個圓柱
    switch lower(REGION)
        case 'all',       base = 'sensor_area_B_3d';
        case 'quadrants', base = 'sensor_quadrant_B_3d';
        otherwise, error('REGION must be ''all'' or ''quadrants''.');
    end
    variants = {'full_assembly_sleeve','full_assembly'};
    tags     = {'sleeve','nosleeve'};

    % ---- Pass 1：兩 variant 內插 + 取兩象限，算共用色階/長度尺度 ----
    S = cell(1,2);  allB = [];
    for v = 1:2
        S{v} = interp_points(fullfile(droot,variants{v},'coil1'), 'coil1', c, ni, SR, SH, NSAMP, NSHOW, REGION);
        allB = [allB; S{v}.Bmag];                                          % [MODIFIED] 回到 |B|（不投影）
        fprintf('%-22s: %4d 支箭頭 | mean |B| = %.4f mT | mean B.n+ = %+.4f mT\n', ...
                variants{v}, numel(S{v}.Bmag), mean(S{v}.Bmag), mean(S{v}.B*ni));
    end
    CLIM  = [min(allB) max(allB)];          % 共用色階（|B|）
    SCALE = (LMAX_F*SR) / CLIM(2);          % 共用：最大 |B| → LMAX_F×R 長
    fprintf('共用 clim(|B|) = [%.3f %.3f] mT | SCALE = %.4g mm/mT\n', CLIM, SCALE);

    for v = 1:2
        render(S{v}, c, ni, SR, SH, CLIM, SCALE, VIEW, ...
               fullfile(figdir, sprintf('%s_%s.png', base, tags{v})));
    end
end

% ============================================================================
function s = interp_points(dir_in, prefix, c, ni, SR, SH, NSAMP, NSHOW, REGION)
% 內插法逐行比照 code/main_function/build_V_matrix.m（LOC 0.30mm / 面積均勻 / 重心內插）
    LOC = 0.30;  OVER = 2*NSAMP;  rng(1);
    d = import_ansys_data(dir_in, 'all', prefix);
    X = [d.x d.y d.z]*1e3;  B3 = [d.bx d.by d.bz]*1e3;          % mm / mT
    loc = vecnorm(X - c.', 2, 2) < LOC;
    P = X(loc,:);  BP = B3(loc,:);
    TR = delaunayTriangulation(P);
    Tb = null(ni.');  t1 = Tb(:,1);  t2 = Tb(:,2);
    a  = SH*rand(OVER,1);  r = SR*sqrt(rand(OVER,1));  th = 2*pi*rand(OVER,1);
    pts = c.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';
    ti = pointLocation(TR, pts);  good = find(~isnan(ti));
    good = good(1:min(NSAMP, numel(good)));
    bc = cartesianToBarycentric(TR, ti(good), pts(good,:));
    conn = TR.ConnectivityList(ti(good),:);
    Bint = zeros(numel(good),3);
    for cc = 1:4, Bint = Bint + bc(:,cc).*BP(conn(:,cc),:); end
    Pg = pts(good,:);

    % [MODIFIED] REGION='all' → 整個圓柱；'quadrants' → 只取 +x+y 與 -x-y 兩個對角象限
    dx = Pg(:,1)-c(1);  dy = Pg(:,2)-c(2);                       % 以 sensor 中心為原點（world xy）
    if strcmpi(REGION,'quadrants')
        keep = [];
        for q = 1:2
            if q==1, idx = find(dx>0 & dy>0); else, idx = find(dx<0 & dy<0); end
            idx = idx(round(linspace(1, numel(idx), min(NSHOW, numel(idx)))));
            keep = [keep; idx];                                              %#ok<AGROW>
        end
    else
        idx  = (1:size(Pg,1)).';
        keep = idx(round(linspace(1, numel(idx), min(NSHOW, numel(idx)))));  % 全圓柱均勻抽 NSHOW 支
    end
    s.P = Pg(keep,:);  s.B = Bint(keep,:);
    s.Bmag = vecnorm(s.B,2,2);
end

% ============================================================================
function render(s, c, ni, SR, SH, CLIM, SCALE, VIEW, outpng)
    colS = [0.85 0.10 0.10];
    f = figure('Color','w','Position',[100 100 1250 1000]); hold on;

    % ---- sensor 圓柱輪廓（P1_full 風格：半透明面 + col*0.6 邊線；alpha 壓低才看得到裡面）----
    % [MODIFIED] 圓柱畫在**取樣實際所在的位置**：build_V_matrix 的 a=axial_tol*rand → pts = c + a*n̂
    %   → 圓柱從 sensor 中心往 +n 長（z ∈ [c, c+H]），不是以 c 為中心。原本畫成置中 → 箭頭浮在柱外。
    [Xc,Yc,Zc] = cylinder(SR,48);  Xc=Xc+c(1); Yc=Yc+c(2);  Zc = Zc*SH + c(3);
    surf(Xc,Yc,Zc,'FaceColor',colS,'FaceAlpha',0.10,'EdgeColor','none');
    th = linspace(0,2*pi,49);
    for zc = [c(3), c(3)+SH]
        patch(c(1)+SR*cos(th), c(2)+SR*sin(th), zc*ones(size(th)), colS, ...
              'FaceAlpha',0.10,'EdgeColor','none');
        plot3(c(1)+SR*cos(th), c(2)+SR*sin(th), zc*ones(size(th)), '-', ...
              'Color',colS*0.6,'LineWidth',1.4);
    end

    % ---- 箭頭 = **完整 3D B 向量**（不投影）：長度 ∝ |B|、顏色 = |B|（共用 SCALE 與 clim）----
    V = s.B * SCALE;
    cmap = turbo(64);
    idx = round( (s.Bmag - CLIM(1)) / max(CLIM(2)-CLIM(1), eps) * 63 ) + 1;
    idx = min(max(idx,1),64);
    for k = 1:64
        sel = idx==k;
        if ~any(sel); continue; end
        quiver3(s.P(sel,1),s.P(sel,2),s.P(sel,3), V(sel,1),V(sel,2),V(sel,3), 0, ...
                'Color',cmap(k,:),'LineWidth',1.2,'MaxHeadSize',0.6);
    end
    colormap(gca,turbo(64));  clim(CLIM);
    cb = colorbar; cb.FontSize=14; cb.FontWeight='bold';
    cb.Label.String='|B| (mT)'; cb.Label.FontWeight='bold'; cb.Label.FontSize=15;

    ax = gca;
    m = 1.45*SR;
    xlim([c(1)-m c(1)+m]); ylim([c(2)-m c(2)+m]); zlim([c(3)-1.2*SH c(3)+2.4*SH]);   % [MODIFIED] 配合圓柱往上長
    grid(ax,'off'); box(ax,'on'); daspect(ax,[1 1 1]); view(ax,VIEW(1),VIEW(2));
    set(ax,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    xlabel(ax,'x (mm)','FontWeight','bold'); ylabel(ax,'y (mm)','FontWeight','bold');
    zlabel(ax,'z (mm)','FontWeight','bold');
    xt=get(ax,'XTick'); if numel(xt)>3; set(ax,'XTick',xt(1:2:end)); end
    yt=get(ax,'YTick'); if numel(yt)>3; set(ax,'YTick',yt(1:2:end)); end
    ax.Toolbar.Visible='off';
    exportgraphics(f, outpng, 'Resolution', 150);  close(f);
    fprintf('saved %s\n', outpng);
end
