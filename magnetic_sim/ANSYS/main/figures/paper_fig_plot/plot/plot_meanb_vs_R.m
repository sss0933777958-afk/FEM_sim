function plot_meanb_vs_R(SRC, force)
% plot_meanb_vs_R -- 工作空間平均場強 mean(||b||) 隨取樣半徑 R 的變化（真實格點）
% =========================================================================
%   定義（使用者拍板 2026-08-18）：**把球內所有格點、六顆激發的 ||b|| 全部加起來，
%   再除以總數量 N_p(R)*6**
%
%       mean||b||(R) = [ sum_j sum_{i: |p_i| <= R} ||b_ij|| ] / ( 6 * N_p(R) )
%
%   —— 即「累積球 |p| <= R」內的池化平均（六顆激發合併成單一曲線），
%   與各極先算平均再對六極平均等價（每個 R 六欄共用同一批點）。
%
%   ⚠ **真實格點、不內插**：點就是 Maxwell 匯出的 .fld 格點（WP 細格步距 0.02 mm），
%     前段沿用校正管線 extract_maxwell_data -> build_actuator_data（濾鐵 + T->mT +
%     all-source 號誌），與 select_ball 的判準一致 —— 所以每個 R 用到的點，
%     就是校正時實際會拿到的那些點。||b|| 是向量範數，不受 actuator 旋轉與號誌影響。
%
%   激發電流 = FEM 激發的 1 A（規則 fit-current-matches-sim），故縱軸即 mT @ 1 A。
%
%   快取 data/meanb_vs_R_<src>.mat 只存每個節點的 (r [µm], s = sum_j ||b_ij|| [mT])，
%   重跑數秒；重讀 .fld 要數十秒。
%
%   風格①粗體框圖：字級 36、框線 4.0、box on、無 grid、單位括號、
%   刻度奇數個等距、水平軸端點只標數字不畫 tick、首末資料點貼齊框邊。
%
%   輸出 → figures/paper_fig/Section2_E/meanb_vs_R_<src>.png（覆蓋迭代）
% =========================================================================
    clc;
    if nargin < 1 || isempty(SRC),   SRC   = 'maxwell'; end   % 'maxwell'（WP 細格 .fld）| 'apdl'（graded .dat）
    if nargin < 2 || isempty(force), force = false;     end

    DR   = 2;                                   % 曲線取點間距 [µm]（只影響解析度）
    RMIN = 0;                                   % 繪圖下限 [µm]（使用者拍板 2026-08-18：從 0 起、R=0 對齊框左緣）
                                                %   ⚠ R < 18.7 µm 球內只有「最靠近原點的那 1 個格點」（r = 1.32 µm），
                                                %     所以最左邊那段是常數 = 原點場強，不是統計平均；R = 30 才有 19 點。
    RMAX = 500;                                 % 繪圖上限 [µm]

    here   = fileparts(fileparts(mfilename('fullpath')));      % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    [rr, ss, N_I] = load_nodes(SRC, here, force);   % rr = |p| [µm]；ss = sum_j ||b_ij|| [mT]

    R  = (RMIN:DR:RMAX).';
    N  = zeros(size(R));   mb = zeros(size(R));
    for k = 1:numel(R)
        q     = rr <= R(k);
        N(k)  = nnz(q);
        % [MODIFIED 2026-08-18 使用者拍板] 曲線**從第一個有格點的 R 起畫**（r_min = 1.32 µm），
        %   不畫 R < r_min 那段（值 = 0 會在左框緣拉出一條陡升的垂直線）。軸仍自 0 起。
        if N(k) > 0, mb(k) = sum(ss(q)) / (N_I * N(k)); else, mb(k) = NaN; end   % 除以 N_p(R) * 6
    end
    ok = N > 0;   R = R(ok);   N = N(ok);   mb = mb(ok);   % 丟掉球內無格點的 R

    fprintf('%s：節點 %d、激發 %d｜mean||b|| = %.4f (R=%d) / %.4f (R=150) / %.4f (R=300) / %.4f mT (R=%d)\n', ...
            upper(SRC), numel(rr), N_I, mb(1), RMIN, mb(R==150), mb(R==300), mb(end), RMAX);
    fprintf('  N_p：%d (R=%d) -> %d (R=150) -> %d (R=%d)｜全域 %.4f ~ %.4f mT\n', ...
            N(1), RMIN, N(R==150), N(end), RMAX, min(mb), max(mb));

    render(R, mb, SRC, DR, figdir);
end

% ============================================================================
function [rr, ss, N_I] = load_nodes(SRC, here, force)
% 每個節點的 |p| [µm] 與「六顆激發的 ||b|| 總和」[mT]（已濾鐵、真實格點）
    cf = fullfile(here, 'data', sprintf('meanb_vs_R_%s.mat', lower(SRC)));
    if exist(cf,'file') && ~force
        S = load(cf);   rr = S.rr;   ss = S.ss;   N_I = S.N_I;
        fprintf('由快取載入 %s（%d 節點 × %d 激發）\n', cf, numel(rr), N_I);
        return
    end
    if strcmpi(SRC,'maxwell')
        ROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
        addpath(fullfile(ROOT,'function'));
        cfg = model_config('long2016_hexapole_halfcut','tip40um');
        raw = extract_maxwell_data(cfg, 'all', cfg.default_variant);
    else
        ROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
        addpath(fullfile(ROOT,'function'), fullfile(ROOT,'common_path'), fullfile(ROOT,'utils'));
        cfg = model_config('long2016_hexapole_halfcut','tip40um');
        raw = extract_ansys_data(cfg, 'all', 'graded');
    end
    ad  = build_actuator_data(raw, cfg);         % actuator frame + 濾鐵 + mT + all-source
    rr  = vecnorm(ad.Pa, 2, 2)*1e6;              % m → µm
    N_I = size(ad.Ba, 3);
    ss  = zeros(numel(rr), 1);
    for k = 1:N_I
        ss = ss + vecnorm(ad.Ba(:,:,k), 2, 2);   % ||b|| 逐激發累加 [mT]
    end
    save(cf, 'rr', 'ss', 'N_I');
    fprintf('已存 %s（%d 節點 × %d 激發）\n', cf, numel(rr), N_I);
end

% ============================================================================
function s = nice_step(x)
% 取「>= x 的最小 nice 步長」（figure-style：等距、疏密適中）
    n = [1 2 2.5 3 4 5 10];
    e = floor(log10(x));   c = x/10^e;
    i = find(n >= c - 1e-12, 1);   s = n(i)*10^e;
end

% ============================================================================
function [yt, s, YL] = ytick_pick(lo, hi)
% 縱軸刻度：**整數**、奇數個(3 或 5)、等距、兩端留白 = 間距（YL = [yt(1)-s, yt(end)+s]），
%   且 YL 必須完整包住資料、資料盡量置中。
%   [MODIFIED 2026-08-18 使用者拍板] 一度改成「自 0 起 + 上緣 8% 裕度」，但本圖資料
%   只在 10.1~11.3 mT，0 基準會讓下方 85% 全空 → 使用者拍板改回「貼資料的整數刻度」
%   （2026-08-18 再拍板：**固定 3 個刻度** -> 本例 YL = [9.5,11.5]、刻度 10/10.5/11、
%     填充率 58%。3 個「整數」刻度只能是 10/11/12 配 YL=[9,13]、填充率剩 29%，已被否決 ->
%     此圖採半整數步長，照 figure-style 對「值域窄的場量」保留 mT + 小數刻度的例外。）
%   ⚠ 整數是硬條件（figure-style 優先序：整數 > 奇數個 > 等距 > 兩端留白），填充率讓步。
    cand = [1 2 2.5 3 4 5 10];
    best = [];
    for nt = 3                                       % [MODIFIED 2026-08-18 使用者拍板] 固定 3 個刻度
        for e = -1:2                                 % 允許 0.5 這類半整數步長（見上：3 個整數刻度會空掉 7 成）
            for c = cand
                s  = c*10^e;
                % yt(1) 必為 s 的整數倍，且 YL=[yt1-s, yt1+nt*s] 要包住 [lo,hi]
                k  = ceil((hi - (nt-1)*s - s)/s) : floor((lo + s)/s);
                if isempty(k), continue; end
                % 在可行的位移中挑「框中心最貼近資料中心」者（取 k(end) 會把資料壓在框底）
                [~, kb] = min(abs((k*s + ((nt-1)/2)*s) - (lo+hi)/2));
                yy  = k(kb)*s + (0:nt-1)*s;
                YLc = [yy(1)-s, yy(end)+s];
                f   = (hi-lo)/diff(YLc);             % 填充率
                if isempty(best) || f > best.f
                    best = struct('f',f, 'yt',yy, 's',s, 'YL',YLc, 'nt',nt);
                end
            end
        end
    end
    yt = best.yt;   s = best.s;   YL = best.YL;
    fprintf('  YTick %d 個、間距 %g、填充率 %.0f%%\n', best.nt, s, 100*best.f);
end

% ============================================================================
function render(R, mb, SRC, DR, figdir)
    FS = 36;  LWBOX = 4.0;
    fig = figure('Color','w','Position',[80 60 1250 950]);  ax = axes(fig);  hold(ax,'on');

    plot(ax, R, mb, '-', 'Color',[0.10 0.35 1.00], 'LineWidth',4.0, 'Clipping','off');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...
           'TickLength',[.015 .015],'TickDir','out');

    % ---- 水平軸：首末資料點貼齊左右框邊；內部 3 個等距整數刻度 ----
    %   figure-style 優先序「整數 > 奇數個 > 等距 > 兩端留白」：端點 40/500 要標數字，
    %   故最外側留白本來就可讓步；取 [125 250 375]（步長 125、整數、奇數個）。
    % [MODIFIED 2026-08-18 使用者拍板] 左端不是 0，而是**第一個球內有格點的 R**
    %   （r_min = 1.32 µm，落在 DR = 2 µm 的取樣格上就是 R = 2）——軸從 0 起會宣稱
    %   有 0 ~ 2 µm 的資料，但那段球內根本沒有格點。
    XL = [R(1) R(end)];   xlim(ax, XL);
    set(ax,'XTick',[125 250 375]);

    % ---- 縱軸：整數刻度、奇數個等距、兩端留白 = 間距（使用者拍板：不用 0 基準）----
    [yt, sy, YL] = ytick_pick(min(mb), max(mb));
    ylim(ax, YL);   set(ax,'YTick',yt);

    % 水平軸端點只標數字、不畫 tick（figure-style 2026-08-06）
    yoff = YL(1) - 0.025*diff(YL);
    for j = 1:2
        text(ax, XL(j), yoff, sprintf('%g',XL(j)), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end

    xlabel(ax, '$\mathbf{R\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, '$\mathbf{\overline{\|b\|}\;(mT)}$', 'Interpreter','latex', 'FontSize',FS);
    fprintf('  XL = [%g %g]；YL = [%g %g]，YTick 間距 %g（%d 個）\n', XL, YL, sy, numel(yt));

    % [MODIFIED 2026-08-18 使用者拍板] 圖例拿掉（單一資料系列，資料源寫在檔頭與圖說即可）。
    ax.Toolbar.Visible = 'off';   hold(ax,'off');

    out = fullfile(figdir, sprintf('meanb_vs_R_%s.png', lower(SRC)));
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s  (DR = %g µm)\n', out, DR);
end
