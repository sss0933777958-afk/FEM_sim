function plot_err_hist_overlay(SRC, Rsel_um, CONV, MODEL, GEOM, VARIANT, ZSUF, NCSET, BINWUT, XCAP)
% plot_err_hist_overlay -- long2016 半切六極 R=150µm：fix vs bias 絕對殘差「疊圖」
% =========================================================================
%   同一取樣球(R=150µm)、同一絕對殘差定義 |b_FEM − S·ĝ_I·K̄·I| (mT)：
%     藍 = 無 bias(fix 單一 ℓ)、紅 = 有 bias(18-param)。
%   共用 fix 的 x 刻度(xlim[0 1]、tick 0.2~0.8 + 起終點)；共用 edges(nb=180)。
%   percentage 縱軸、圖例標顏色對應、**不標平均線**、無軸標題。
%   輸出 → figures/paper_fig/Section2_E/err_hist_overlay.png(覆蓋迭代)。
% =========================================================================
    clc;
    if nargin < 1, SRC = 'apdl'; end       % [ADDED] 'apdl'(ANSYS .dat) | 'maxwell'(.fld) → 檔名加 _maxwell
    if nargin < 2 || isempty(Rsel_um), Rsel_um = 150; end   % [ADDED] 取樣球半徑 [µm]；≠150 時檔名加 _R<NNN>
    % [ADDED 2026-08-13] CONV=true：**改用該 R 的 N_c 減量設計校正**（判準＝l_hat 穩定
    %   ∧ g_I 穩定 ∧ K_I 符合物理，各持續 10 步），再拿該 R 內全部真實格點當評估集。
    %   兩個模型各自有自己的 N_c。檔名加 _conv。
    if nargin < 3 || isempty(CONV), CONV = false; end
    % [ADDED 2026-08-27] model/geom/variant + fullref-cache suffix, so the same figure
    %   can be produced for zhi_peng. Defaults reproduce the long2016 behaviour exactly.
    %   VARIANT must be given explicitly for zhi_peng: its cfg.default_variant is still
    %   the superseded 'maxwell' (gI ~ 7.0) instead of 'maxwell_split'.
    if nargin < 4 || isempty(MODEL),   MODEL   = 'long2016_hexapole_halfcut'; end
    if nargin < 5,                     GEOM    = 'tip40um';                   end
    if nargin < 6,                     VARIANT = '';                          end
    if nargin < 7 || isempty(ZSUF),    ZSUF    = '';                           end
    % [ADDED 2026-08-27] NCSET = [N_single N_eighteen] overrides the point count taken from
    %   the *_fullref sweep. Needed when the wanted R is off the swept grid (that sweep runs
    %   R = 20:20:500, so R=150 is not on it and the neighbours can disagree).
    if nargin < 8 || isempty(NCSET),   NCSET   = [NaN NaN];                     end
    % [ADDED 2026-08-27] Shared bin width in uT. 2.8 suits long2016 (residuals reach
    %   ~0.94 mT). zhi_peng reaches ~9.5 mT, so 2.8 would give 3382 bins and leave the
    %   single-parameter peak bin with only ~17 counts (24% relative noise) -> use 17 uT
    %   there. NOTE: two figures drawn with different bin widths are NOT height-comparable.
    if nargin < 9 || isempty(BINWUT),  BINWUT  = 2.8;                            end
    % [ADDED 2026-08-27] XCAP: force the residual-axis upper limit [mT] instead of letting
    %   xlim_pick derive it. Anything past XCAP is clipped from view (the count is printed).
    if nargin < 10, XCAP = []; end
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section2_E');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    CAL = solver_path(SRC);                                          % [ADDED] 依 SRC 切分支(內含 rmpath 另一分支防遮蔽)
    sstr = ''; if strcmpi(SRC,'maxwell'), sstr = '_maxwell'; end     % [ADDED] 來源後綴

    % ---- 前段 pipeline(只做一次)----
    cfg = model_config(MODEL, GEOM);
    if isempty(VARIANT), VARIANT = cfg.default_variant; end
    raw = load_raw(SRC, cfg, VARIANT);                                        % [MODIFIED] apdl='graded' .dat / maxwell=.fld
    ad  = build_actuator_data(raw, cfg);   Pc_base = ad.Pc_base;
    F   = zeros(6, cfg.N_I);  for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
    Rum = Rsel_um;                                          % [MODIFIED] 取樣半徑可調（原寫死 150）
    [P, Bstack, npts] = cfg.select_ball(ad, Rum*1e-6);
    fprintf('取樣範圍 R <= %d um：npts = %d\n', Rum, npts);

    % ---- 兩模型各擬合一次 → 逐節點×激發 絕對殘差 (mT) ----
    if CONV
        addpath(fullfile(CAL,'utils'), fullfile(CAL,'utils','long2016_hexapole_halfcut'));
        n1_ = NCSET(1);   if ~isfinite(n1_), n1_ = []; end
        n2_ = NCSET(2);   if ~isfinite(n2_), n2_ = []; end
        ca = conv_set(Rum, cfg, false, CAL, n1_, MODEL, GEOM, VARIANT, ZSUF);
        cb = conv_set(Rum, cfg, true,  CAL, n2_, MODEL, GEOM, VARIANT, ZSUF);
        fprintf('N_c：single %d 點｜eighteen %d 點（評估集 = %d 個真實格點）\n', ...
                ca.npts, cb.npts, npts);
        err0 = fit_abs_resid(P, Bstack, Pc_base, F, npts, false, ca);
        err1 = fit_abs_resid(P, Bstack, Pc_base, F, npts, true,  cb);
    else
        err0 = fit_abs_resid(P, Bstack, Pc_base, F, npts, false);   % single_param (無 bias)
        err1 = fit_abs_resid(P, Bstack, Pc_base, F, npts, true);    % eighteen_param (有 bias)
    end
    cv0 = std(err0)/mean(err0)*100;   cv1 = std(err1)/mean(err1)*100;   % CV = σ/μ [%]（不受 bin 寬影響）
    fprintf('single_param  : mean=%.4f mT  max=%.4f mT  CV=%.1f%%\n', mean(err0), max(err0), cv0);
    fprintf('eighteen_param: mean=%.4f mT  max=%.4f mT  CV=%.1f%%\n', mean(err1), max(err1), cv1);

    % ---- 兩組**共用同一個 bin 寬**（使用者拍板 2026-08-19）--------------------
    %   原本兩組各自用自己的值域鋪滿 nb=180 個 bin，bin 寬差 8.57 倍
    %   （single 5.258 µT / eighteen 0.613 µT）→ 兩個直方圖的**視覺面積不相等**
    %   （面積 = 100% × bin 寬），高度就不能互比，這正是審查意見指出的問題。
    %   共用 bin 寬後兩者面積都是 100%×BINW，高度差才真正代表集中程度。
    %   BINW 取 2.8 µT：eighteen 分到約 40 個 bin（形狀看得出單峰 + 右尾），
    %   single 峰值仍有約 100 筆樣本（雜訊 ~10%）。兩組樣本數相同且固定
    %   （1771 格點 × 6 激發 = 10626 筆），所以「紅色 bin 變多」與「藍色變平滑」
    %   本質上互相排斥，2.8 µT 是折衷點。
    BINW = BINWUT*1e-3;                                   % 共用 bin 寬 [mT]
    hi   = max([err0(:); err1(:)]);
    edg  = 0 : BINW : (ceil(hi/BINW)*BINW);          % 共用 edges
    ctr0 = (edg(1:end-1)+edg(2:end))/2;   ctr1 = ctr0;
    pct0 = histcounts(err0, edg) / numel(err0) * 100;
    pct1 = histcounts(err1, edg) / numel(err1) * 100;
    fprintf('共用 bin 寬 %.2f uT：全域 %d 個 bin｜eighteen 佔 %d 個｜single 峰值 %.2f%% (~%.0f 筆)、eighteen 峰值 %.2f%% (~%.0f 筆)\n', BINW*1e3, numel(edg)-1, ceil(max(err1)/BINW), max(pct0), max(pct0)/100*numel(err0), max(pct1), max(pct1)/100*numel(err1));

    [cB, cR] = pick_bar_colors(Rum);   % [MODIFIED] 長條配色依取樣半徑分組（同 plot_err_hist 的 pick_colors）
    % [MODIFIED 2026-09-02 套繪圖規則 1/2/3/6] 刻度 60、圖例 45、軸標題 36、框線 5.0；
    %   畫布改成**正方形英吋**（14.5 in，與 ell_vs_npts / gain_vs_npts 那批同尺寸，
    %   60 pt 的刻度數字相對大小才一致）。像素尺寸會被 MATLAB 靜默縮放，故用英吋。
    FS = 60;  FSLAB = 36;  FSLEG = 45;  LWBOX = 5.0;  CANV = 14.5;
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig);  hold(ax,'on');
    % [MODIFIED 2026-08-20 使用者拍板] **黑色線段全部拿掉**：兩組都不畫邊，
    %   原本紅色的 stairs 頂部輪廓 + 每 4 根一條的垂直黑邊也一併移除。
    h0 = bar(ax, ctr0, pct0, 1, 'FaceColor',cB, 'FaceAlpha',0.60, 'EdgeColor','none');   % 亮藍(single)
    h1 = bar(ax, ctr1, pct1, 1, 'FaceColor',cR, 'FaceAlpha',0.60, 'EdgeColor','none');   % 紅(eighteen)
    % [MODIFIED 2026-08-17 使用者拍板] **不再畫 mean 虛線**（連同其圖例條目一併移除）。
    %   mean 仍算出來、印在 console 供追溯；要恢復就把下面兩行 xline 取消註解並加回圖例。
    mu0 = mean(err0);  mu1 = mean(err1);
    fprintf('  mean：single %.4f mT / eighteen %.4f mT（圖上不再標虛線）\n', mu0, mu1);
    % ml0 = xline(ax, mu0, '--', 'Color',[0.00 0.00 0.00], 'LineWidth',2.8);   % single mean（黑）
    % ml1 = xline(ax, mu1, '--', 'Color',[0.00 0.60 0.00], 'LineWidth',2.8);   % eighteen mean（深綠）

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX,'TickLength',[.015 .015],'TickDir','out');
    % [MODIFIED 2026-08-20 使用者拍板] 橫軸端點必須是**乾淨的數**（整數或 0.N）。
    %   ⬆ **作廢**：2026-08-19 曾改成「右端畫到資料結束處 `edg(end)`」-> 端點變成
    %   **0.946**，是直接從資料搬來的尾數，違反端點規則。改回用 xlim_pick 的
    %   nice 上界（maxE=0.946 -> [0, 1]），內部刻度仍是 0.2/0.4/0.6/0.8。
    [xr_, xt_] = xlim_pick(SRC, max(err0), XCAP);
    if ~isempty(XCAP)
        nc0_ = sum(err0 > XCAP);   nc1_ = sum(err1 > XCAP);
        fprintf(['  x-axis capped at %g mT: %d/%d single (%.2f%%) and %d/%d eighteen '...
                 '(%.2f%%) residuals fall outside' newline], XCAP, ...
                nc0_, numel(err0), 100*nc0_/numel(err0), ...
                nc1_, numel(err1), 100*nc1_/numel(err1));
    end
    xlim(ax, xr_);   set(ax,'XTick', xt_);
    % [REVERTED 2026-08-03] 縱軸維持原作法（使用者拍板；同 plot_err_hist）：
    %   ylim 貼齊資料最大值(進位到 0.1)、tick 取 linspace 的 3 個內縮值。
    % [MODIFIED 2026-08-27] linspace(0,ytop,5) then inner 3 gives 1.8/3.6/5.4 whenever ytop
    %   is not a nice multiple of 4, which breaks figure-style's integer-tick rule. Use the
    %   documented 'axis from zero' recipe: nice step, 8% headroom, odd tick count.
    % [MODIFIED 2026-09-02 使用者指定] 縱軸**四根刻度，且與框底 0 / 框頂 T 等距**。
    %   舊寫法是「上緣留 8% 裕度、刻度取 (1:ny)*sy」-> 末刻度到框頂的距離跟刻度間距
    %   完全不同（志鵬實測 1.2 -> 1.318，只有 0.118，而間距是 0.4），違反規則 5。
    %   正解：[0,T] 平分 n+1 段，刻度 = (1:n)*s、T = (n+1)*s，n = 4。
    %   s 取 0.1 的整數倍（刻度數字才會是整數或一位小數），在 [smin, 1.35*smin]
    %   窗口裡挑最漂亮者：整數 > 0.5 的倍數 > 其餘。
    ymax = max([pct0 pct1]);   ny = 4;
    smin = 1.02*ymax/(ny+1);
    k0 = max(1, ceil(smin/0.1 - 1e-9));   k1 = max(k0, ceil(1.35*smin/0.1));
    %   ⚠ 不可「整數優先」—— 那會在 ymax=7.06 時挑到 s=2（框頂 10、填充率僅 71%，
    %     上方空掉近三成）。改成「窗口內**最小**的 0.5 倍數」（整數自然含在內），
    %     沒有就退回窗口下界：ymax=7.06 -> s=1.5（框頂 7.5、填充 94%）、
    %     ymax=1.22 -> s=0.3（框頂 1.5、填充 81%）。
    bk = k0;
    for kk = k0:k1
        if mod(kk,5) == 0, bk = kk;   break; end
    end
    sy   = bk*0.1;
    ytop = (ny+1)*sy;
    ylim(ax,[0 ytop]);   set(ax,'YTick', round((1:ny)*sy*10)/10);
    % 起點 0 + 終點：只標數字、不畫 tick(左右角落補字)
    xr = xlim(ax);
    text(ax, xr(1), -0.022*ytop, sprintf('%g',xr(1)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');
    text(ax, xr(2), -0.022*ytop, sprintf('%g',xr(2)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');
    hLegR = bar(ax, NaN, NaN, 'FaceColor',cR, 'FaceAlpha',0.60, 'EdgeColor','none');   % 只給圖例的紅色色塊（也不加黑邊）
    % [MODIFIED 2026-08-19 使用者拍板] 圖例改放**座標框內右上角**（原本在框外上方）；
    %   每列一個系列、mean 併入該列（不畫 mean 虛線）。右上角本來就是空的，省版面。
    lg = legend([h0 hLegR], ...
                ... % [MODIFIED 2026-08-20 使用者拍板] 圖例拿掉 Mean（同 shell）；數值仍印 console
                {'Single parameter', 'Eighteen parameters'}, ...
                'Interpreter','tex', 'Location','northeast', 'NumColumns',1);
    lg.FontSize = FSLEG;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;   % [規則2/7]
    lg.Color = 'w';
    xlabel(ax, '$\mathbf{Residual\;(mT)}$', ...
           'Interpreter','latex', 'FontSize',FSLAB);   % 絕對殘差軸標題
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$', ...
           'Interpreter','latex', 'FontSize',FSLAB);   % 百分比縱軸標題
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    % [MODIFIED 2026-08-03] 檔名一律明確標「求解器_半徑」（overlay 本身含 single+eighteen，不標模型）
    cstr = ''; if CONV, cstr = '_conv'; end
    out = fullfile(figdir, sprintf('err_hist_overlay%s_%s_R%d%s.png', ...
                                   cstr, lower(SRC), Rum, ZSUF));
    % [MODIFIED 2026-09-02 規則 6] exportgraphics 會裁掉畫布白邊 -> 正方形匯出後不等邊。
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, out, '-dpng', '-r200');            % 14.5 in x 200 dpi = 2900 px 見方
    fprintf('wrote %s\n', out);
end

% ============================================================================
function [c0, c1] = pick_bar_colors(Rum)
% [ADDED 2026-08-03] 長條配色每個取樣半徑一組（使用者拍板，與 plot_err_hist 的 pick_colors 同組）：
%   R150：藍(single) / 紅(eighteen)
%   其他(R300…)：紫 #7B52AB (single) / 橘 #E69F00 (eighteen)
%   ⚠ mean 虛線不隨之改——一律黑(single)/深綠(eighteen)，中性色在兩組長條上都看得清。
    if Rum == 150
        c0 = [0.05 0.10 0.95];      c1 = [0.85 0.10 0.10];   % [MODIFIED 2026-08-20] 藍統一成深藍
    else
        c0 = [0.482 0.322 0.671];   c1 = [0.902 0.624 0.000];
    end
end

% ============================================================================
function CAL = solver_path(SRC)
% [ADDED] 依 SRC 掛上對應分支、並「移除另一分支」——兩分支有同名函式，只 addpath 會被遮蔽。
    % [MODIFIED 2026-08-27] Tree moved under matlab\Flux\; the old hardcoded paths are gone.
    MAIN = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));  % ...\ANSYS\main
    APDL = fullfile(MAIN,'matlab','Flux','APDL','Calibration_using_FEM_modeling');
    MW   = fullfile(MAIN,'matlab','Flux','Maxwell');
    switch lower(SRC)
        case 'apdl',    CAL = APDL;  OTHER = MW;
        case 'maxwell', CAL = MW;    OTHER = APDL;
        otherwise, error('SRC 必為 ''apdl'' | ''maxwell''');
    end
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(OTHER,'function'));  rmpath(fullfile(OTHER,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));   addpath(fullfile(CAL,'common_path'));
    fprintf('[SRC=%s] model_config -> %s\n', SRC, which('model_config'));
end

% ============================================================================
function raw = load_raw(SRC, cfg, VARIANT)
% [ADDED] apdl = graded .dat / maxwell = .fld
    if strcmpi(SRC,'maxwell')
        raw = extract_maxwell_data(cfg, 'all', VARIANT);
    else
        raw = extract_ansys_data(cfg, 'all', 'graded');
    end
end

% ============================================================================
function [xr, xt] = xlim_pick(SRC, maxE, XCAP)
    % [ADDED 2026-08-27] explicit upper limit: 4 equidistant internal ticks, centred.
    %   8/5 = 1.6 is not an integer step, so per figure-style's priority order
    %   (integer > count > equidistant > end padding) we keep integers and give up
    %   'padding == spacing': XCAP=8 -> ticks 1/3/5/7 (step 2, padding 1).
    if nargin >= 3 && ~isempty(XCAP)
        % [MODIFIED 2026-09-02 規則 5 + 9] 原本固定 4 根、步長由 XCAP/5 取 nice ->
        %   XCAP=8 得到 1/3/5/7（步長 2、但與兩端只隔 1，違反「與端點間距相等」）。
        %   改成 [0,XCAP] 平分 n+1 段：s = XCAP/(n+1)，**n 優先四根**（使用者指定，
        %   與其他圖的水平軸例外一致）；XCAP=8 -> s=1.6 -> 1.6/3.2/4.8/6.4（一位小數，合規則 9）；
        %   s 必須是 **0.1 的整數倍**（規則 9：水平軸小數只到一位），整數 s 更優先。
        xr = [0 XCAP];   xt = [];   best = -1;
        for n = [4 3 5 7 6]
            s = XCAP/(n+1);
            if abs(s*10 - round(s*10)) > 1e-9, continue; end     % 非 0.1 的倍數 -> 跳過
            sc = 2*(abs(s-round(s)) < 1e-9);                     % 整數 s 加分（根數已由順序決定）
            if sc > best, best = sc;   xt = round((1:n)*s*10)/10;   end
            if n == 4 && ~isempty(xt), break; end                % 四根優先，找到就用
        end
        if isempty(xt), xt = (1:4)*(XCAP/5); end                  % 保險
        return;
    end
% [ADDED] 橫軸(殘差 mT)：apdl 沿用定案 fix 刻度；maxwell 依 single 殘差最大值自動等距 4 內縮 tick。
    if strcmpi(SRC,'apdl'), xr = [0 1.0];  xt = [0.2 0.4 0.6 0.8];  return; end
    % [FIXED 2026-08-03] 原本用 `while xr(2)<maxE, s=nice_step(s*1.3); end` 找上界 —— **無窮迴圈**：
    %   nice_step 會把 s*1.3 吸附回同一格（1→1.3→1、2.5→3.25→2.5、5→6.5→5），s 永遠不變。
    %   R=300 踩到此路徑 → MATLAB 空轉數小時、無輸出。改成在候選清單上直接往上找，保證終止。
    cand = [1 2 2.5 5 10];
    k0 = floor(log10(max(maxE, realmin)/5));   s = [];
    for k = k0:(k0+3)
        for c = cand
            if 5*c*10^k >= maxE, s = c*10^k; break; end
        end
        if ~isempty(s), break; end
    end
    if isempty(s), s = maxE/5; end             % 保險（理論上不會走到）
    xr = [0, 5*s];
    xt = (1:4)*s;
    % [ADDED 2026-08-03] 長尾資料時 5×步長 會讓右端留白過大（違反 figure-style「留白不可大於間距」）：
    %   上界 >1.25×maxE 才收緊 → 改「n 個等距 tick + 半格留白」。R150 等既有圖不受影響（比值僅 1.05）。
    if xr(2) > 1.25*maxE
        s2 = nice_step(maxE/3.5);
        n  = max(3, floor(maxE/s2));
        up = (n + 0.5)*s2;
        while up < maxE, n = n + 1;  up = (n + 0.5)*s2;  end   % n 遞增，必終止
        xr = [0, up];   xt = (1:n)*s2;
    end
end

function s = nice_step(x)
% [ADDED] nice 等距步長（1/2/2.5/5/10 × 10^k）
    k = floor(log10(x));   m = x/10^k;
    cand = [1 2 2.5 5 10];
    [~,i] = min(abs(cand - m));
    s = cand(i)*10^k;
end


% ---- 擬合 + 逐節點×激發 絕對殘差 |S·G − Bstack| (mT) ----
function cal = conv_set(Rum, cfg, USE_BIAS, CAL, NFORCE, MODEL, GEOM, VARIANT, ZSUF)
% [MODIFIED 2026-08-28] Sampler switched to sample_axes_shells (the six actuator-axis
%   directions on Nr equally spaced shells + the centre -> N = 6*Nr+1). The mixed ehgap
%   ladder (ax / rg / gf / eh) it used before is retired; design_by_npts below is kept
%   only so an older cache can still be replayed, and is no longer called from here.
%   N_c comes from NFORCE (the settled values: long2016 single 31 / eighteen 25;
%   zhi_peng 13 / 7). Without NFORCE we fall back to the axsh R-sweep cache.
    here_ = fileparts(fileparts(mfilename('fullpath')));
    N_ = [];
    if nargin >= 5 && ~isempty(NFORCE), N_ = NFORCE; end
    if isempty(N_)
        df_ = fullfile(here_, 'data', ['full_vs_conv_vs_R_maxwell_axsh' ZSUF '.mat']);
        assert(exist(df_,'file')==2, ...
               ['no NCSET given and %s is missing -- pass the point count explicitly ' ...
                '(long2016 31/25, zhi_peng 13/7)'], df_);
        D_  = load(df_);
        NCv = D_.n_c1;   FBv = D_.fb_c1;
        if USE_BIAS, NCv = D_.n_c2;  FBv = D_.fb_c2; end
        Rv_ = D_.R_um(:).';   ok_ = isfinite(NCv(:).') & ~logical(FBv(:).');
        if ~any(ok_), ok_ = isfinite(NCv(:).'); end
        cand_ = find(ok_);   [~, jj_] = min(abs(Rv_(cand_) - Rum));
        N_ = NCv(cand_(jj_));
    end
    assert(isfinite(N_) && mod(N_-1,6)==0, ...
           'axes-shells point count must be 6*Nr+1; got N=%g', N_);
    Nr_ = (N_ - 1)/6;
    fprintf(['  conv design: R=%g um, axes-shells Nr=%d -> N_c=%d (bias=%s)' newline], ...
            Rum, Nr_, N_, mat2str(USE_BIAS));

    addpath(fullfile(fileparts(fileparts(here_)), 'temp_code', 'scripts'));
    o_ = struct('model',MODEL, 'geom',GEOM, ...
                'variant',VARIANT, 'frame','actuator', 'quiet',true);
    Pq_ = conv_design_ws(Nr_, Rum*1e-6, struct('points_only',true,'R_act',cfg.R_act,'quiet',true));
    evalc('[Pd_, Bd_] = conv_design_ws([], Rum*1e-6, setfield(o_,''query'',Pq_));');
    assert(size(Pd_,1) == N_, 'axes-shells returned %d points, expected %d', size(Pd_,1), N_);

    F_ = zeros(6, cfg.N_I);
    for j_ = 1:cfg.N_I, F_(cfg.apdl_to_paper_idx(j_), j_) = 1; end
    [e_, l_] = fitting(Pd_, Bd_, cfg.Pc_base, 0.5e-3, USE_BIAS);
    [~, ~, G_] = solve_current(l_, e_, cfg.Pc_base, Pd_, Bd_, F_);
    cal = struct('e',e_, 'l_hat',l_, 'G',G_, 'npts',size(Pd_,1));
end

% ============================================================================
function [P, Bs] = design_by_npts(KIND, SPEC, nc, R, o, cfg)   %#ok<DEFNU>
% Recover a design from the mixed ladder (ax / rg / gf / eh) by its POINT COUNT.
% [RETIRED 2026-08-28] conv_set now builds the design with sample_axes_shells; this is
%   kept only for replaying an old *_fullref cache.
    P = [];  Bs = [];
    for q = 1:numel(KIND)
        switch KIND{q}
            case 'eh',  Pq = sample_equal_h(R, SPEC{q}, struct('quiet',true));
            case 'rg',  Pq = sample_rings(  R, SPEC{q}, struct('quiet',true));
            case 'ax',  Pq = sample_axes4(R, cfg.R_act, ...
                                 struct('quiet',true,'naxes',SPEC{q},'plus',true));
            % [REMOVED 2026-09-02] case 'gf'（等測度球格）已隨等測度取樣法一併廢除。
            otherwise,  continue;
        end
        if isempty(Pq) || size(Pq,1) ~= nc, continue; end
        try,  evalc('[Pk,Bk] = conv_design_ws(Pq, R, o);');  catch, continue;  end
        if size(Pk,1) == nc,  P = Pk;  Bs = Bk;  return;  end
    end
end

% ============================================================================
function err = fit_abs_resid(P, Bstack, Pc_base, F, npts, USE_BIAS, cal)
% [MODIFIED 2026-08-23] 可選第 7 引數 cal：**直接用 main.m 產的校正結果**
%   （e / l_hat / G，由 N_c 減量設計擬合而來），再拿 (P, Bstack)（該 R 內全部
%   真實格點）當評估集。不給就是原行為（同一組點校正 + 評估）。
    if nargin >= 7 && ~isempty(cal)
        e = cal.e;   l_hat = cal.l_hat;   G = cal.G;
    else
        [e, l_hat, ~] = fitting(P, Bstack, Pc_base, 0.5e-3, USE_BIAS);
        [~, ~, G]     = solve_current(l_hat, e, Pc_base, P, Bstack, F);
    end
    Pc = make_Pc(e, Pc_base);
    pbar = P / l_hat;  S = zeros(3*npts, 6);
    for k = 1:6, d = pbar - Pc(:,k).'; r3 = sum(d.^2,2).^1.5; S(:,k) = reshape((d./r3).', 3*npts, 1); end
    resid = S*G - Bstack;
    err = zeros(npts, 6);
    for j = 1:6, rj = reshape(resid(:,j), 3, npts); err(:,j) = sqrt(sum(rj.^2,1)).'; end
    err = err(:);
end

% ---- 電荷格 Pc = Pc_base + E(e)(與 solve_current 內一致；fix e=0 → Pc_base)----
function Pc = make_Pc(e17, Pc_base)
    if isempty(e17) || all(e17(:)==0), Pc = Pc_base; return; end
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
