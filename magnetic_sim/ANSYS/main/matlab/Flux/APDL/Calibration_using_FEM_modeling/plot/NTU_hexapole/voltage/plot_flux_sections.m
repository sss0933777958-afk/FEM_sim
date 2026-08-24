function plot_flux_sections()
% plot_flux_sections — 通量積分用的兩個截面 + 微小面積分割（檢查用，純幾何）。
%   ① 導柱截面 @ z = 15.00 mm：圓盤 r2.5、圓心 (20.5, 0)、法線 n = +z → Φ = ∫Bz dA
%   ② 磁極板截面 @ x = 13.15 mm：矩形 y ∈ [-2.9942, 2.9942] × z ∈ [8.80, 9.05]、n = +x → Φ = ∫Bx dA
%   分割：邊長 H 的正方形微小面積；圓盤取「格心落在 r<2.5 內」的格；矩形**強制整除鋪滿**
%   （N = round(L/H)、實際格距 L/N）—— 否則格子鋪不滿會直接變成面積誤差（板厚 0.25mm 用
%   H=0.10 會得到 3 排 = 0.30mm、面積多算 +20%，踩過）。
%   z=15 處只有導柱（襯套 flange 從 15.60 才開始）→ 兩個 variant 幾何相同。
%   出 2 檔 → figures/shared/flux_section_{post,pole}.png

    % ======================== 參數 ========================
    H     = 0.05;                 % 微小面積邊長 [mm]（收斂測試：0.05→0.01 通量只變 0.3~0.8%）
    XCUT  = 13.15;                % 磁極板切面 x [mm]
    ZCUT  = 15.00;                % 導柱切面 z [mm]
    PC    = [20.5, 0];  PR = 2.5; % 導柱圓心 / 半徑 [mm]
    ZB    = 8.80;  ZT = 9.05;     % 板底 / 板頂 [mm]
    % =====================================================

    here   = fileparts(mfilename('fullpath'));
    vbase  = fileparts(fileparts(here));
    figdir = fullfile(vbase,'figures','shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath(fullfile(vbase,'code','function'));       % ntu_pole_profile

    % ---- 板半寬 w(XCUT)：由 9 段輪廓內插（不寫死）----
    prof = ntu_pole_profile();
    lo = prof(prof(:,2)<0,:);  [~,i] = sort(lo(:,1));  lo = lo(i,:);
    W  = abs(interp1(lo(:,1), lo(:,2), XCUT));

    % ================= ① 導柱圓盤 @ z = ZCUT =================
    g = -PR+H/2 : H : PR;   [GX,GY] = ndgrid(g,g);
    in = GX.^2 + GY.^2 <= PR^2;                       % 格心在圓內
    cx = PC(1)+GX(in);  cy = PC(2)+GY(in);
    f = figure('Color','w','Position',[100 100 1000 980]);  ax = axes(f);  hold(ax,'on');
    draw_cells(cx, cy, H, H, [0.55 0.75 0.95]);
    th = linspace(0,2*pi,400);
    plot(ax, PC(1)+PR*cos(th), PC(2)+PR*sin(th), '-', 'Color',[0.15 0.25 0.45], 'LineWidth',2.5);
    style2d(ax, 'x (mm)', 'y (mm)', true);   % [MODIFIED] 橫軸畫的是 cx → 標 x（原本標反）
    xlim(ax, PC(1)+[-PR-0.35 PR+0.35]);  ylim(ax, PC(2)+[-PR-0.35 PR+0.35]);
    title(ax, sprintf('post section @ z = %.2f mm   |   %d cells of %.2f \\times %.2f mm   |   A = %.3f mm^2', ...
          ZCUT, nnz(in), H, H, nnz(in)*H^2), 'FontSize',15,'FontWeight','bold');
    exportgraphics(f, fullfile(figdir,'flux_section_post.png'), 'Resolution',150);  close(f);
    fprintf('post @ z=%.2f : %d 格 × %.2fmm → A=%.4f mm^2（理論 %.4f，誤差 %+.2f%%）\n', ...
            ZCUT, nnz(in), H, nnz(in)*H^2, pi*PR^2, 100*(nnz(in)*H^2-pi*PR^2)/(pi*PR^2));

    % ================= ② 磁極板矩形 @ x = XCUT =================
    Nz = round((ZT-ZB)/H);   dz = (ZT-ZB)/Nz;         % 強制整除鋪滿
    Ny = round(2*W/H);       dy = 2*W/Ny;
    gz = ZB+dz/2 : dz : ZT-dz/2;
    gy = -W+dy/2 : dy : W-dy/2;
    [QY,QZ] = ndgrid(gy,gz);
    f = figure('Color','w','Position',[100 100 1500 640]);  ax = axes(f);  hold(ax,'on');
    draw_cells(QY(:), QZ(:), dy, dz, [0.55 0.75 0.95]);
    plot(ax, [-W W W -W -W], [ZB ZB ZT ZT ZB], '-', 'Color',[0.15 0.25 0.45], 'LineWidth',2.5);
    style2d(ax, 'y (mm)', 'z (mm)', false);           % z 軸拉伸（長寬比 24:1，等比例會變一條線）
    xlim(ax, [-W-0.15 W+0.15]);  ylim(ax, [ZB-0.03 ZT+0.03]);
    title(ax, sprintf(['pole section @ x = %.2f mm   |   %d \\times %d = %d cells of %.4f \\times %.4f mm' ...
          '   |   A = %.3f mm^2   (z axis stretched)'], XCUT, Ny, Nz, Ny*Nz, dy, dz, Ny*Nz*dy*dz), ...
          'FontSize',14,'FontWeight','bold');
    exportgraphics(f, fullfile(figdir,'flux_section_pole.png'), 'Resolution',150);  close(f);
    fprintf('pole @ x=%.2f : %d×%d=%d 格 (%.4f×%.4f mm) → A=%.4f mm^2（理論 %.4f，誤差 %+.2f%%）\n', ...
            XCUT, Ny, Nz, Ny*Nz, dy, dz, Ny*Nz*dy*dz, 2*W*(ZT-ZB), 100*(Ny*Nz*dy*dz-2*W*(ZT-ZB))/(2*W*(ZT-ZB)));
end

% ============================================================================
function draw_cells(cx, cy, dx, dy, col)
% 一次畫 N 個格子（單一 patch、Faces/Vertices）：格心 (cx,cy)、尺寸 dx×dy
    n = numel(cx);
    V = [cx(:)-dx/2 cy(:)-dy/2; cx(:)+dx/2 cy(:)-dy/2; cx(:)+dx/2 cy(:)+dy/2; cx(:)-dx/2 cy(:)+dy/2];
    F = [(1:n).', (n+1:2*n).', (2*n+1:3*n).', (3*n+1:4*n).'];
    patch('Faces',F, 'Vertices',V, 'FaceColor',col, 'FaceAlpha',0.45, ...
          'EdgeColor',[0.35 0.35 0.35], 'LineWidth',0.25);
end

function style2d(ax, xl, yl, eq)
    if eq, daspect(ax,[1 1 1]); end
    box(ax,'on'); grid(ax,'off');
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    xlabel(ax, xl, 'FontWeight','bold');  ylabel(ax, yl, 'FontWeight','bold');
    xt=get(ax,'XTick'); if numel(xt)>3; set(ax,'XTick',xt(1:2:end)); end
    yt=get(ax,'YTick'); if numel(yt)>3; set(ax,'YTick',yt(1:2:end)); end
    ax.Toolbar.Visible='off';
end
