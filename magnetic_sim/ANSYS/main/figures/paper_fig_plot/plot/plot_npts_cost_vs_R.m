function plot_npts_cost_vs_R(USE_BIAS, SRC, SAMPLING)
% plot_npts_cost_vs_R -- paper 圖：取樣點數 N 與擬合 cost J 隨取樣半徑 R 的趨勢
% =========================================================================
%   上 panel = 取樣點數 N(R)（無單位 → 不標單位）
%   下 panel = cost J(R) [mT²]
%       J = Σ_{j=1..6} ‖ S(ℓ̂,Pc)·g_j − b_j ‖²，g_j = (SᵀS)⁻¹Sᵀb_j（變數投影消去電荷強度）
%       即 fitting.m 最小化並回傳的目標函數；求和涵蓋取樣球內所有點×3 分量×6 激發。
%   ⚠ J 是外延量：會隨 N 一起長大（上 panel 正是給讀者看這件事）。
%
%   **本圖只讀 plot_ell_gain_vs_R 產生的 sweep 快取、不重算**（快取需含 J_R；
%   若舊快取沒有 J_R，先跑一次 plot_ell_gain_vs_R 重建）。
%   風格 = 選項①粗體框 @ paper scale，與 ell_gain_vs_R 完全一致（同字體/框/刻度規則）。
%   輸出 → figures/paper_fig/Section2_E/npts_cost_vs_R[_bias][_maxwell].png（覆蓋迭代）。
% =========================================================================
    clc;
    if nargin < 1, USE_BIAS = false; end   % false = fix(無 e)；true = 18-param bias
    if nargin < 2, SRC = 'maxwell';  end   % 'apdl' | 'maxwell'
    if nargin < 3 || isempty(SAMPLING)     % 與 ell_gain_vs_R 同慣例：maxwell 預設 interp
        if strcmpi(SRC,'maxwell'), SAMPLING = 'interp'; else, SAMPLING = 'raw'; end
    end
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    bstr = ''; if USE_BIAS, bstr = '_bias'; end
    sstr = ''; if strcmpi(SRC,'maxwell'),      sstr = '_maxwell'; end
    istr = ''; if strcmpi(SAMPLING,'interp'),  istr = '_interp';  end

    % ---- 讀 sweep 快取（與 plot_ell_gain_vs_R 同一顆）----
    cachef = fullfile(here, 'data', sprintf('ell_gain_sweep%s%s%s.mat', sstr, istr, bstr));
    assert(exist(cachef,'file')==2, ...
        'sweep 快取不存在：%s\n先跑 plot_ell_gain_vs_R(%d,''%s'')。', cachef, USE_BIAS, SRC);
    S = load(cachef);
    assert(isfield(S,'J_R'), ...
        '快取缺 J_R（舊版）：%s\n刪掉它並重跑 plot_ell_gain_vs_R(%d,''%s'') 重建。', cachef, USE_BIAS, SRC);
    Rum = S.Rum;  npts_R = S.npts_R;  J_R = S.J_R;
    fprintf('loaded %s\n', cachef);

    % ---- 繪圖範圍（與 ell_gain_vs_R 同：maxwell 從 80µm 起，小 R 格點太少擬合病態）----
    if strcmpi(SRC,'maxwell'), Rmin = 80; else, Rmin = 40; end
    m = Rum >= Rmin;   Rum = Rum(m);  npts_R = npts_R(m);  J_R = J_R(m);
    % [MODIFIED 2026-08-03] x 軸同樣「兩端留白 = tick 間距」（使用者拍板；原本 xlim 從 100 起，
    %   既讓首尾 tick 貼框(留白 0)，又把 R=80 那點切在框外）。
    %   資料 80–500 → 5 個奇數等距 tick 100..500、s=100 → xlim=[0,600]，R=80 回到圖內。
    [~, XT] = axlim_auto(Rum(1), Rum(end), 5);
    % [MODIFIED 2026-08-03] 水平軸 = 資料範圍本身：起點 80、**終點 500**（使用者拍板）。
    %   ⚠ 因此 x 不套「兩端留白 = 間距」（該規則只適用縱軸）：末 tick 500 落在右框線上、
    %   左側留 20。tick 仍為 5 個等距（100..500），與縱軸 tick 數一致。
    XL = [Rum(1), Rum(end)];

    % [MODIFIED 2026-08-03] single 與 eighteen 兩張圖**共用同一組 y 尺度**（使用者拍板；
    %   figure-style「同類比較圖共用軸尺度」）——否則各自 auto-scale，讀者無法直接比 cost 大小。
    %   做法：一併讀另一個 model 的快取，y 範圍/刻度都用兩者的聯集去算（缺檔就退回只用自己）。
    [npts_all, J_all] = union_with_other(here, sstr, istr, bstr, Rmin, npts_R, J_R);

    % [ADDED 2026-08-03] 上下兩 panel 的 y tick「數量必須一致」（使用者拍板）：
    %   先挑一個兩邊都可行的奇數 tick 數（取「最差 panel 的填充率最高」者），再各自算範圍。
    NY = pick_common_n(npts_all, J_all, [3 5]);

    fprintf('  R=%d: N=%d, J=%.4g mT^2   →   R=%d: N=%d, J=%.4g mT^2\n', ...
            Rum(1), npts_R(1), J_R(1), Rum(end), npts_R(end), J_R(end));

    % ---- 畫圖(2×1 panel,選項①粗體框；與 ell_gain_vs_R 逐項一致)----
    cN = [0.05 0.10 0.95];   cJ = [0.85 0.10 0.10];     % N 亮藍 / J 紅
    FS = 30;
    fig = figure('Color','w','Position',[100 80 980 940]);
    t = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

    % 上：N(R)
    ax1 = nexttile(t);  hold(ax1,'on');
    plot(ax1, Rum, npts_R, '-o', 'Color',cN, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cN);
    style_panel(ax1, FS);
    xlim(ax1, XL);   set(ax1,'XTick',XT);
    [yl,yt] = axlim_auto(min(npts_all), max(npts_all), NY);   ylim(ax1,yl);  set(ax1,'YTick',yt);   % [MODIFIED] 共用尺度
    ax1.XTickLabel = {};                                % 上 panel x 數字隱藏(共用軸)
    ylabel(ax1, '$\mathbf{Number\;of\;points}$', 'Interpreter','latex', 'FontSize',36);   % 無單位 → 不標
    hold(ax1,'off');

    % 下：J(R)
    ax2 = nexttile(t);  hold(ax2,'on');
    plot(ax2, Rum, J_R, '-s', 'Color',cJ, 'LineWidth',3.0, 'MarkerSize',8, 'MarkerFaceColor',cJ);
    style_panel(ax2, FS);
    xlim(ax2, XL);   set(ax2,'XTick',XT);
    [yl,yt] = axlim_auto(min(J_all), max(J_all), NY);         ylim(ax2,yl);  set(ax2,'YTick',yt);   % [MODIFIED] 共用尺度
    ylabel(ax2, '$\mathbf{Cost\;(mT^{2})}$', 'Interpreter','latex', 'FontSize',36);
    xlabel(ax2, '$\mathbf{Sampling\;range\;(micro\;meter)}$', 'Interpreter','latex', 'FontSize',36);   % [MODIFIED 2026-08-05] 單位改拼字 micro meter
    hold(ax2,'off');

    % [MODIFIED 2026-08-03] 檔名一律明確標「求解器_模型」（使用者拍板）
    mstr = 'single'; if USE_BIAS, mstr = 'eighteen'; end
    out = fullfile(figdir, sprintf('npts_cost_vs_R_%s_%s.png', lower(SRC), mstr));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('\nwrote %s\n', out);
end

% ============================================================================
function [npts_all, J_all] = union_with_other(here, sstr, istr, bstr, Rmin, npts_R, J_R)
% [ADDED 2026-08-03] 把「另一個 model（single↔eighteen）」的 sweep 併進來，讓兩張圖共用 y 尺度。
%   另一顆快取不存在時，靜默退回只用自己的資料（圖仍可產、只是不保證跨圖同尺度）。
    npts_all = npts_R(:);   J_all = J_R(:);
    if isempty(bstr), other = '_bias'; else, other = ''; end
    f = fullfile(here, 'data', sprintf('ell_gain_sweep%s%s%s.mat', sstr, istr, other));
    if exist(f,'file') ~= 2
        fprintf('  [warn] 另一 model 快取不存在(%s) → y 尺度只依本圖資料\n', f);   return;
    end
    O = load(f);
    if ~isfield(O,'J_R') || ~isfield(O,'npts_R') || ~isfield(O,'Rum')
        fprintf('  [warn] 另一 model 快取缺欄位 → y 尺度只依本圖資料\n');          return;
    end
    mo = O.Rum >= Rmin;                                  % 同一段 R 才可比
    npts_all = [npts_all; O.npts_R(mo).'];
    J_all    = [J_all;    O.J_R(mo).'];
    fprintf('  y 尺度共用：併入 %s（J max 本圖 %.4g / 另一圖 %.4g mT^2）\n', ...
            other_name(other), max(J_R), max(O.J_R(mo)));
end

function s = other_name(other)
    if isempty(other), s = 'single'; else, s = 'eighteen'; end
end

% ============================================================================
function n = pick_common_n(v1, v2, nlist)
% [ADDED 2026-08-03] 兩 panel 共用的 tick 數：取「兩邊都可行」且「較差那個 panel 的填充率最高」者。
%   填充率 = 資料跨度 / 軸總跨距（越高代表留白越少、曲線越舒展）。
    best = nlist(1);  bestFill = -inf;
    for n_ = nlist
        [L1, ~] = axlim_auto(min(v1), max(v1), n_);
        [L2, ~] = axlim_auto(min(v2), max(v2), n_);
        f = min( (max(v1)-min(v1))/diff(L1), (max(v2)-min(v2))/diff(L2) );
        if f > bestFill, bestFill = f;  best = n_; end
    end
    n = best;
end

% ============================================================================
function [lim, tk] = axlim_auto(lo, hi, nlist)
% [ADDED 2026-08-03] 通用軸範圍/刻度：**奇數個等距 tick，且兩端留白 = tick 間距**（使用者拍板）。
%   ⇒ 軸總跨距 = (n+1)·s。反過來先定 s：在 nice 候選(1/2/2.5/3/4/5/10 ×10^k) × 給定的奇數 tick 數中，
%   取「總跨距最小」且容得下 [lo,hi] 的組合（總跨距最小 ⇒ 自動避開負值區與多餘留白）。
%   x、y 共用此邏輯：y 給 nlist=3、x 給 5（長軸可多放）。
    cand = [1 2 2.5 3 4 5 10];
    mid  = (lo+hi)/2;   rng = max(hi-lo, realmin);
    best = {};   bestSpan = inf;
    for n = nlist
        k0 = floor(log10(rng/(n+1)));
        for k = k0:(k0+4)
            hit = false;
            for c = cand
                s   = c*10^k;
                ctr = round(mid/s)*s;                    % tick 對齊 s 的整數倍（數值好讀）
                T = { ctr + (-(n-1)/2 : (n-1)/2)*s };
                % [ADDED 2026-08-03] 資料從 ~0 起時，優先讓「最低 tick = 0」——避免計數/成本軸
                %   出現負值刻度，也就是使用者說的「資料最低點停在 0 tick 上」那個樣子。
                if lo >= 0 && lo < s, T = [{(0:n-1)*s}, T]; end %#ok<AGROW>
                for it = 1:numel(T)
                t   = T{it};
                L   = [t(1)-s, t(end)+s];                % 兩端留白 = s
                % [MODIFIED 2026-08-03] 不只要「在範圍內」，還要**與上下框線有淨空**——
                %   否則像 cost 這種小 R 時 ≈0 的資料會整條貼在底框上（使用者打槍過）。
                %   淨空門檻 0.15·s：資料最低/最高點離框至少 0.15 格。
                clr = 0.15*s;
                if lo >= L(1)+clr && hi <= L(2)-clr
                    if (n+1)*s < bestSpan, bestSpan = (n+1)*s;  best = {L, t}; end
                    hit = true;  break;
                end
                end   % T 迴圈
                if hit, break; end
            end
            if hit, break; end
        end
    end
    if isempty(best)                                     % 保險（理論上不會走到）
        n = nlist(1);  s = rng/(n+1);
        t = mid + (-(n-1)/2 : (n-1)/2)*s;   best = {[t(1)-s, t(end)+s], t};
    end
    lim = best{1};   tk = best{2};
end

% （原 local `nice_step` 與 `ylim_auto` 已於 2026-08-03 移除：軸範圍/刻度全部走 axlim_auto，
%   兩者已無任何呼叫點，留著只會讓人誤以為改它會生效。）

% ============================================================================
function style_panel(ax, FS)
% 選項①粗體框 @ paper scale(tick 由呼叫端明設,不減半)。
    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.02 .02]);   % 刻度數字 Helvetica 粗體(維持原字體)
    ax.Toolbar.Visible = 'off';
end
