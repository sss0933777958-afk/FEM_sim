function plot_sensor_cyl_B_3d(face, SHARED)
% plot_sensor_cyl_B_3d -- Hall 感測圓柱內 500 個撒點的 B 向量 3D 圖(Maxwell)
% =========================================================================
%   三個貼附面各一張（每面用該極自己的自激發場）：
%     'p1flat' → sensor_cylinder_B_3d_P1_flat_maxwell.png （P1 半切平切上表面，n̂⁺=+z）
%     'p1cone' → sensor_cylinder_B_3d_P1_cone_maxwell.png （P1 半切件底錐面，n̂⁺ 朝下出鋼）
%     'p2cone' → sensor_cylinder_B_3d_P2_cone_maxwell.png （P2 上錐面，n̂⁺ 朝上出鋼）
%   圓柱 R0.15 × H0.10mm、離面氣隙 0.41mm、SOFF 4.572mm、每柱 500 點。
%   箭頭 = B 方向(等長)、依 ‖b‖ 分 28 個 turbo bin 上色；紅圈 = 圓柱上下緣。
%   座標系 = **WP 框、mm**。
%
%   ⚠ 這是**內插**(Maxwell 格距 0.1mm)，與校正管線 V_METHOD='scattered' 同一套。
%     已驗證：500 點的內插四面體**無一含鋼件頂點**（撒點離最近鋼格點 0.385~0.524mm）。
%   ★ 撒點與場向量**一律取自 plot_sensor_B_hist 的 per-pole 快取**
%     （`sensor_B_hist_P<k>_maxwell_soff<S>_n<N>.mat`），本檔不自建快取、也不讀 .fld
%     → 3D 圖與直方圖保證是同一組點。快取不存在時先跑 `plot_sensor_B_hist(true)`。
%
%   plot_sensor_cyl_B_3d()               % 三面都畫
%   plot_sensor_cyl_B_3d('p2cone')       % 只畫 P2 上錐面
%   plot_sensor_cyl_B_3d([], true)       % 三面共用色階（預設 false = 各自貼齊自己範圍）
%
%   風格①粗體框圖 3D 變體 A（同尺度立方：daspect[1 1 1] + box off + 手動框邊省最近角）。
% =========================================================================
    clc;
    if nargin < 1 || isempty(face),   face   = 'all';  end
    if nargin < 2 || isempty(SHARED), SHARED = false;  end

    here   = fileparts(fileparts(mfilename('fullpath')));       % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    ALL = {'p1flat','p1cone','p2cone'};
    OUT = struct('p1flat','P1_flat', 'p1cone','P1_cone', 'p2cone','P2_cone');
    if strcmpi(face,'all'), cs = ALL; else, cs = {lower(face)}; end
    if ~all(ismember(cs, ALL)), error('face 必為 %s 或 ''all''', strjoin(ALL,' | ')); end

    D = struct();
    for i = 1:numel(cs), D.(cs{i}) = load_case(cs{i}, here); end

    CL = [];
    if SHARED                                   % 同類比較圖共用色階(figure-style)
        lo = inf; hi = -inf;
        for i = 1:numel(cs)
            b = vecnorm(D.(cs{i}).B,2,2);  lo = min(lo,min(b));  hi = max(hi,max(b));
        end
        CL = [lo hi];
    end
    for i = 1:numel(cs)
        b = vecnorm(D.(cs{i}).B,2,2);
        fprintf('%-7s ‖b‖ %.3f~%.3f mT (mean %.4f)\n', cs{i}, min(b), max(b), mean(b));
        render(D.(cs{i}).P, D.(cs{i}).B, CL, figdir, OUT.(cs{i}));
    end
end

% ============================================================================
function C = load_case(name, here)
% 由 plot_sensor_B_hist 的 per-pole 快取取出該面的撒點與場向量。
    kpole = 1;  if strncmpi(name,'p2',2), kpole = 2; end
    f = fullfile(here,'data', sprintf('sensor_B_hist_P%d_maxwell_soff4.572_n500.mat', kpole));
    if ~exist(f,'file')
        error('load_case:noCache', ['缺快取 %s\n' ...
              '先跑 plot_sensor_B_hist([],[],[],[],''%s'') 產生（P2 需讀 2GB .fld）。'], ...
              f, ternary_str(kpole==2,'p1p2','p1'));
    end
    S = load(f);
    if ~isfield(S, ['P_' name])
        error('load_case:noFace', '快取 %s 內找不到面 ''%s''', f, name);
    end
    C.P = S.(['P_' name]);   C.B = S.(['B_' name]);
end

% ============================================================================
function s = ternary_str(c,a,b), if c, s=a; else, s=b; end, end

% ============================================================================
function render(P, B, CL, figdir, tag)
% P: 500×3 撒點 [mm, WP 框]；B: 500×3 [mT]
    FS = 36;  NB = 28;                                   % 字級 36、turbo 28 bin(figure-style)
    LARR = 0.035;                                        % 箭頭等長 [mm]（≈柱高 1/3，包在上下緣內）
    HALF = 0.18;                                         % 立方框半邊長 [mm]
    bm = vecnorm(B,2,2);
    % 色階範圍 = 資料範圍。
    %   ⚠ figure-style 的「clim([0 CLIM])、CLIM=ceil(max/50)*50」是為 0~數百 mT 的磁路場圖訂的；
    %     本圖 ‖b‖ 只在 10.4~11.1mT 的窄帶，套 [0,50] 會讓 500 支箭頭全變同色、柱內結構消失
    %     （已實測）。故此處保留資料範圍——**colorbar 的樣式**（平滑 turbo、粗體 36、‖b‖ 標題、
    %     兩端有刻度數字）仍完全照規則。
    if isempty(CL), CL = [min(bm) max(bm)]; end

    cmid = mean(P,1);                                    % 圓柱幾何中心
    cc   = round(cmid/0.1)*0.1;                          % 框心吸附到 0.1mm 格(刻度才是整數)
    XL = cc(1)+[-HALF HALF];  YL = cc(2)+[-HALF HALF];  ZL = cc(3)+[-HALF HALF];

    fig = figure('Color','w','Position',[100 100 1150 860]);  ax = axes(fig);  hold(ax,'on');
    % ⚠ colormap 要**平滑 turbo**（colorbar 才是連續漸層）；28 bin 只用在箭頭 truecolor 上。
    %   餵 turbo(NB) 會讓 colorbar 變成 28 塊粗色塊 —— 不符規則（source of truth: plot_circuit_side.m:302）。
    colormap(ax, turbo);  clim(ax, CL);

    % ---- 圓柱外框:上下緣紅圈 + 半透明側面(撒點的實際取樣體積)----
    draw_cylinder(ax, P);

    % ---- 箭頭:等長、依 ‖b‖ 分 bin 上色 ----
    cmap = turbo(NB);  edg = linspace(CL(1), CL(2), NB+1);
    U = B ./ bm * LARR;                                  % 單位方向 × 固定長度
    for k = 1:NB
        if k < NB, m = bm >= edg(k) & bm < edg(k+1); else, m = bm >= edg(k) & bm <= edg(k+1); end
        if ~any(m), continue; end
        quiver3(ax, P(m,1),P(m,2),P(m,3), U(m,1),U(m,2),U(m,3), 0, ...
                'Color',cmap(k,:), 'LineWidth',1.0, 'MaxHeadSize',0.6, 'AutoScale','off');
    end

    % ---- 風格①3D 變體 A:同尺度立方、手動框邊(省最近角 3 邊)----
    grid(ax,'off');  box(ax,'off');  daspect(ax,[1 1 1]);
    xlim(ax,XL); ylim(ax,YL); zlim(ax,ZL);
    view(ax, -37.5, 22);                                 % 先設 view(draw_box_edges3 用 campos)
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015]);
    set(ax,'XTick', cc(1)+(-1:1)*0.1, 'YTick', cc(2)+(-1:1)*0.1, 'ZTick', cc(3)+(-1:1)*0.1);
    draw_box_edges3(ax, XL, YL, ZL, 3.0);                % ⚠ 必在 hold off 之前
    % colorbar：3D 透視圖預設會把標題擠出畫布 → 手動縮軸 + 明確定位(figure-style)
    cb = colorbar(ax);  style_cbar(cb, FS);              % 刻度用 MATLAB 自動(等距、兩端有數字)
    ax.Toolbar.Visible = 'off';  hold(ax,'off');
    drawnow;
    % ⚠ 左邊界不可壓太小：z 刻度標籤畫在座標軸左外側，left<0.05 會把 "0.1" 切成 "1"
    %   （看起來像被抽了 ×10^-1，其實是字元被畫布裁掉）。
    set(ax, 'Position', [0.10 0.06 0.60 0.90]);
    set(cb, 'Position', [0.760 0.10 0.030 0.82]);
    % 明寫刻度文字，順便擋掉 MATLAB 對小數值自動抽指數的行為。
    set(ax, 'XTickLabel', compose('%g', get(ax,'XTick')), ...
            'YTickLabel', compose('%g', get(ax,'YTick')), ...
            'ZTickLabel', compose('%g', get(ax,'ZTick')));

    out = fullfile(figdir, sprintf('sensor_cylinder_B_3d_%s_maxwell.png', tag));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function draw_cylinder(ax, P)
% 由撒點反推圓柱的軸/半徑/高，畫上下緣紅圈 + 半透明側面。
    c0 = mean(P,1);
    [~,~,V] = svd(P - c0, 0);   n = V(:,3).';            % 最小變異方向 = 柱軸(高 0.10 << 直徑 0.30)
    s  = (P - c0) * n.';   n = n * sign(max(s)+min(s)+eps);   s = (P - c0)*n.';
    lo = c0 + min(s)*n;    hi = c0 + max(s)*n;           % 上下緣中心
    d  = P - c0 - s.*n;    R  = max(vecnorm(d,2,2));
    t1 = null(n);  t2 = t1(:,2).';  t1 = t1(:,1).';
    th = linspace(0,2*pi,181).';
    Clo = lo + R*(cos(th)*t1 + sin(th)*t2);
    Chi = hi + R*(cos(th)*t1 + sin(th)*t2);
    RED = [0.60 0.00 0.00];
    surf(ax, [Clo(:,1) Chi(:,1)], [Clo(:,2) Chi(:,2)], [Clo(:,3) Chi(:,3)], ...
         'FaceColor',[1 0.75 0.75], 'FaceAlpha',0.14, 'EdgeColor','none');   % 側面
    fill3(ax, Clo(:,1),Clo(:,2),Clo(:,3), [1 0.80 0.80], 'FaceAlpha',0.16, 'EdgeColor','none');
    fill3(ax, Chi(:,1),Chi(:,2),Chi(:,3), [1 0.80 0.80], 'FaceAlpha',0.16, 'EdgeColor','none');
    plot3(ax, Clo(:,1),Clo(:,2),Clo(:,3), '-', 'Color',RED, 'LineWidth',2.2);   % 下緣
    plot3(ax, Chi(:,1),Chi(:,2),Chi(:,3), '-', 'Color',RED, 'LineWidth',2.2);   % 上緣
end

% ============================================================================
function draw_box_edges3(ax, XL, YL, ZL, lw)
% 立方框 12 邊，省「離相機最近角」相連的 3 邊 → 9 邊(標準開口箱)。
%   ⚠ 與 plot_p1_pole_full.m / plot_sphere_lattice_3d.m 同一套(el>=0 分支)；改動要同步。
    [Xc,Yc,Zc] = ndgrid(XL,YL,ZL);  C = [Xc(:) Yc(:) Zc(:)];
    E = [];  for i = 1:8, for j = i+1:8
        if nnz(abs(C(i,:)-C(j,:)) > 1e-9) == 1, E = [E; i j]; end
    end, end
    cp = campos(ax);  [~,near] = min(sum((C - cp).^2, 2));
    for k = 1:size(E,1)
        a = E(k,1);  b = E(k,2);
        if a == near || b == near, continue; end
        plot3(ax, C([a b],1), C([a b],2), C([a b],3), 'k-', 'LineWidth', lw);
    end
end

% ============================================================================
function style_cbar(cb, FS)
% colorbar：標題標準數學字體 ‖b‖ (mT)、刻度數字 Helvetica 粗體（同 plot_circuit_side）。
    cb.FontSize = FS;  cb.FontWeight = 'bold';
    cb.Label.Interpreter = 'latex';  cb.Label.String = '$\mathbf{\|b\|\;(mT)}$';  cb.Label.FontSize = FS;
end
