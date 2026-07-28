function plot_sensor_scan_x(SIDE)
% plot_sensor_scan_x — sensor 截面沿 x 掃描的 ⟨B·n+⟩ 曲線：襯套導磁 vs 不導磁（coil1 = P1 自激發）。
%   把 sensor 圓柱（R0.15 × H0.10、y 中心 -0.99）沿 x 滑過感測區，每個位置算它的讀值 ⟨B·n+⟩
%   → 橫軸 x、縱軸 ⟨B·n+⟩ (mT)。取樣區 = y ∈ [-1.14,-0.84]、厚 0.10mm ＝ sensor 自己的截面。
%
%   [ADDED] SIDE：sensor 放在磁極板的哪一面（板 z ∈ [8.80, 9.05]）
%     'above'：z = 9.46（板頂上 0.41mm）、**n+ = +z**（離板朝上）→ sensor_scan_x.png
%     'below'：z = 8.39（板底下 0.41mm）、**n+ = -z**（鏡像：同樣背離板面）→ sensor_scan_x_below.png
%     'both'（預設）：兩張一起出，**四條曲線共用 ylim** → 上/下也能直接比（per figure-style 同類共用軸尺度）
%   ⚠ 上下**不對稱**：WP 在 z=9.30（板的上方那側），所以下方不是上方的鏡像結果。
%
%   取值法逐行比照 code/main_function/build_V_matrix.m：中心 0.30mm 內真實節點建 Delaunay
%   → 圓柱內面積均勻撒 NS 點 → pointLocation + cartesianToBarycentric → B(p)=Σλ_i·B_i → 投影 n+ 取平均。
%   （不直接呼叫 build_V_matrix：它一次跑 6 sensor × 6 coil、且每次呼叫重載資料，掃 177 位置會重載 177 次。）
%
%   襯套不導磁 = 現有 full_assembly（無襯套幾何）—— μr=1 的襯套 ≡ 空氣 ≡ 無襯套，物理等價。
%   出 2 檔 → figures/shared/sensor_scan_x_{sleeve,nosleeve}.png（**共用 ylim**，兩張才能直接比）

    % ======================== 參數 ========================
    % [MODIFIED] 感測區範圍 = x 0.35 ~ 17.96（0.353 = 極尖端；18.0 = 導柱 r2.5@x20.5 的內緣）。
    %   XS 是**圓柱中心**位置 → 取 0.50 ~ 17.81，使 footprint(±R0.15) 剛好落在 [0.35, 17.96] 內。
    %   若中心放到 17.96，圓柱邊緣會到 18.11 → 吃進導柱鋼裡（Bz 暴衝 −21 mT，非漏磁）。
    XS = 0.50:0.10:17.81;        % 掃描位置（圓柱中心）[mm]
    YC = -0.99;                  % sensor 中心 y [mm]
    SR = 0.15;   SH = 0.10;      % sensor 圓柱 R / H [mm]（= cnst.SENSOR_R / SENSOR_H）
    NS = 1000;                   % 每位置內插取樣點數（= build_V_matrix 的 n_uniform）
    LOC = 0.30;                  % 建 Delaunay 的局部節點球半徑 [mm]（同 build_V_matrix）
    % =====================================================

    if nargin < 1 || isempty(SIDE), SIDE = 'both'; end            % [MODIFIED] 預設一次出上下兩張
    switch lower(SIDE)
        case 'above', sides = {'above'};
        case 'below', sides = {'below'};
        case 'both',  sides = {'above','below'};
        otherwise, error('SIDE must be ''above'' | ''below'' | ''both''.');
    end

    here   = fileparts(mfilename('fullpath'));
    vbase  = fileparts(fileparts(here));
    figdir = fullfile(vbase,'figures','shared');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % import_ansys_data
    droot = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';

    variants = {'full_assembly_sleeve','full_assembly'};
    labels   = {'sleeve: ferromagnetic (\mu_r = 280)','sleeve: non-magnetic (\mu_r = 1)'};
    cols     = [0.20 0.65 0.35; 0.12 0.18 0.55];
    XSENS    = 13.15;                       % sensor 標稱位置（紅點標它）

    % ---- Pass 1：兩面 × 兩 variant 全部算出來（四條曲線）----
    Y = cell(numel(sides), 2);  seg = cell(1, numel(sides));
    for s = 1:numel(sides)
        switch sides{s}
            case 'above', ZC = 9.46;  NI = [0;0;+1];
            case 'below', ZC = 8.39;  NI = [0;0;-1];      % 板底 8.80 - 0.41
        end
        fprintf('--- SIDE = %s : sensor z = %.2f mm, n+ = [%+d %+d %+d] ---\n', sides{s}, ZC, NI);
        for v = 1:2
            Y{s,v} = scan(fullfile(droot,variants{v},'coil1'), 'coil1', XS, YC, ZC, NI, SR, SH, NS, LOC);
            [~,i0] = min(abs(XS-XSENS));
            fprintf('  %-22s: x=%.2f -> <B.n+> = %+.4f mT\n', variants{v}, XSENS, Y{s,v}(i0));
        end
        % 畫「主峰 -> B=0（或掃描右端）」：x<6.55 那段 sensor 在板緣外側空氣裡（w(x)=0.99 交點 6.547）
        [~, ipk] = max(Y{s,2});
        iz = numel(XS);
        for v = 1:2
            k = find(Y{s,v}(ipk:end) < 0, 1, 'first');
            if ~isempty(k), iz = min(iz, ipk + k - 2); end
        end
        seg{s} = ipk:iz;
        if iz < numel(XS)
            fprintf('  peak x=%.2f ; B=0 at x~%.2f -> plot %.2f ~ %.2f mm\n', XS(ipk), XS(iz+1), XS(ipk), XS(iz));
        else
            fprintf('  peak x=%.2f ; never crosses 0 -> plot %.2f ~ %.2f mm (scan end)\n', XS(ipk), XS(ipk), XS(iz));
        end
    end

    % ---- 共用 ylim（四條曲線的繪製區段合併）→ 上/下也能直接比 ----
    allY = [];
    for s = 1:numel(sides), for v = 1:2, allY = [allY; Y{s,v}(seg{s})]; end, end   %#ok<AGROW>
    YL = [min(allY) max(allY)];  pad = 0.08*diff(YL);  YL = [YL(1)-pad YL(2)+pad];
    fprintf('共用 ylim = [%.3f %.3f] mT\n', YL);

    % ---- Pass 2：逐面渲染 ----
    for s = 1:numel(sides)
        sfx = '';  if strcmp(sides{s},'below'), sfx = '_below'; end
        render_overlay(XS(seg{s}), {Y{s,1}(seg{s}), Y{s,2}(seg{s})}, cols, labels, XSENS, YL, ...
                       fullfile(figdir, ['sensor_scan_x' sfx '.png']));
    end
end

% ============================================================================
function bn = scan(dir_in, prefix, XS, YC, ZC, ni, SR, SH, NS, LOC)
    rng(1);  OVER = 2*NS;
    d = import_ansys_data(dir_in, 'all', prefix);
    X = [d.x d.y d.z]*1e3;  B3 = [d.bx d.by d.bz]*1e3;          % mm / mT
    % 預濾到條狀區（大幅加速逐位置的球搜尋）
    keep = X(:,1)>min(XS)-LOC-0.1 & X(:,1)<max(XS)+LOC+0.1 & ...
           abs(X(:,2)-YC)<LOC+SR+0.1 & abs(X(:,3)-ZC)<LOC+SH+0.1;
    X=X(keep,:); B3=B3(keep,:);
    Tb = null(ni.');  t1=Tb(:,1); t2=Tb(:,2);      % [MODIFIED] ni 由呼叫端給（above:+z / below:-z）
    bn = nan(numel(XS),1);
    for k = 1:numel(XS)
        c = [XS(k); YC; ZC];
        loc = vecnorm(X - c.', 2, 2) < LOC;
        if nnz(loc) < 8; continue; end
        P = X(loc,:);  BP = B3(loc,:);
        TR = delaunayTriangulation(P);
        a = SH*rand(OVER,1);  r = SR*sqrt(rand(OVER,1));  th = 2*pi*rand(OVER,1);   % 面積均勻
        % [MODIFIED] 對齊 build_V_matrix：a=rand*[0,H] → 圓柱從 c 往 +n 長（原本我寫 a-SH/2 = 置中，不一致）
        pts = c.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';
        ti = pointLocation(TR, pts);  good = find(~isnan(ti));
        if isempty(good); continue; end
        good = good(1:min(NS, numel(good)));
        bc = cartesianToBarycentric(TR, ti(good), pts(good,:));
        conn = TR.ConnectivityList(ti(good),:);
        Bi = zeros(numel(good),3);
        for cc = 1:4, Bi = Bi + bc(:,cc).*BP(conn(:,cc),:); end
        bn(k) = mean(Bi*ni);
    end
end

% ============================================================================
function render_overlay(XS, Y, cols, labels, XSENS, YL, outpng)
% 兩條曲線疊圖 + sensor 標稱位置紅圓點（各曲線上一顆）。風格：粗體框選項①。
    f = figure('Color','w','Position',[100 100 1240 880]);  ax = axes(f);  hold(ax,'on');
    % [MODIFIED] 移除灰色 y=0 參考線（使用者要求）
    h = gobjects(1,numel(Y));
    for v = 1:numel(Y)
        h(v) = plot(ax, XS, Y{v}, '-', 'Color',cols(v,:), 'LineWidth',2.5);
    end
    % sensor 位置紅圓點
    [~,i0] = min(abs(XS-XSENS));
    hs = gobjects(1,1);
    for v = 1:numel(Y)
        hs = plot(ax, XS(i0), Y{v}(i0), 'o', 'MarkerSize',11, 'MarkerFaceColor',[0.85 0.10 0.10], ...
                  'MarkerEdgeColor','k', 'LineWidth',1.2);
    end
    xlim(ax,[min(XS) max(XS)]);  ylim(ax, YL);      % [MODIFIED] ylim 由呼叫端給（四條共用）
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    xt=get(ax,'XTick'); if numel(xt)>3; set(ax,'XTick',xt(1:2:end)); end
    yt=get(ax,'YTick'); if numel(yt)>3; set(ax,'YTick',yt(1:2:end)); end
    xlabel(ax,'x (mm)','FontWeight','bold');
    ylabel(ax,'$\langle B\cdot n^{+}\rangle$ (mT)','Interpreter','latex','FontWeight','bold','FontSize',18);
    lg = legend([h hs], [labels, {'sensor position'}], 'Location','northeast');
    set(lg,'FontSize',14,'FontWeight','bold','Box','on');
    ax.Toolbar.Visible='off';
    exportgraphics(f, outpng, 'Resolution', 150);  close(f);
    fprintf('saved %s\n', outpng);
end
