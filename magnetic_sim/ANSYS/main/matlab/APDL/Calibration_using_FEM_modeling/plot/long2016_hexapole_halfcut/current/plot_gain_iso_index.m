function plot_gain_iso_index()
% PLOT_GAIN_ISO_INDEX  逐節點性能量 C / kappa 的**直方圖**（fix 與 bias × gain/iso × 兩半徑）。
%   球內每個真實 FEM 節點 p 各自算致動增益矩陣 ^BG_I(p)=S(p)·Ĥ_I 的 SVD（σ₁≥σ₂≥σ₃），
%   依 singular_value.pdf 定義：
%     C(p)     = ∏σ_k = σ₁σ₂σ₃   （flux-generating ellipsoid 體積，(mT/A)³）
%     kappa(p) = σ_min/σ_max = σ₃/σ₁   （isotropy，≤1、→1 等向）
%   顏色：R≤150（含掃描 R_k）深藍 navy [0.05 0.12 0.40]、R≤500 棕 [0.55 0.35 0.17]；R150 fit 參數所有半徑共用。
%   效率：只 select_ball(500e-6) 算一次逐點 C/κ + 每點半徑 r；任何 R_k（含 R150）= 用 r<=R_k 過濾（巢狀子集，等價 select_ball(R_k)）。
%   標準直方圖（原始資料、零事後平滑）：長條=真實計數；固定 bin 數 135；左上黑框 N/mean/CV(%)、紅虛線在 mean、白邊；x 軸各自尺度（nice_axis 乾淨整數刻度）。
%   疊圖（overlay）：淺藍(R_k [0.45 0.68 0.90])+淺紅(R500 固定 [0.95 0.55 0.55]) 疊一起、x 軸固定用 R500 尺度、共用 edges、半透明+白邊；legend 動態標半徑；「僅此」= 無統計框、無均值線。
%   掃描疊圖：R500 固定，R_k=150:50:500（8 檔 _R{k}um）→ 逐張看 R_k 分布隨範圍長大到 = R500（x 軸恆定可比）。
%   輸出（每 dir）：{gain,iso}_index.png（R150 navy）、_R500um.png（R500 棕）、_overlay.png（headline=R_k150）、_overlay_R{150..500}um.png（掃描 8 張）。
%   ★ 兩模型都在 actuator 框同一組節點算（R_act·dhat=Pc_base）：fix 電荷=Pc_base（在軸）、bias=make_Pc(ê)（離軸）。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
    here  = fileparts(mfilename('fullpath'));
    nofix = fileparts(fileparts(here));
    calroot = fileparts(nofix);
    addpath(fullfile(nofix,'code','main_function'));

    cnst = mt_constants();
    apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

    %% ---- 載入資料 + 固定校正參數（R150 fit，兩半徑共用）----
    D = load_coils_actuator('long2016_hexapole_halfcut', cnst, apdl_to_paper_idx, 'all', 'gap_200um');
    Pc_base = D.Pc_base;

    Sf = load(fullfile(calroot,'current_base','data','fit_fixl_R150um_gap200um_mueq.mat'), 'ell','gB','Khat');
    Sb = load(fullfile(nofix,        'data','fit_bias_R150um_gap200um_mueq.mat'), 'ell','gB','Khat','e_hat');
    ellf = Sf.ell*1e-6;  Hf = Sf.gB*Sf.Khat;  Pcf = Pc_base;                     % fix：在軸電荷
    ellb = Sb.ell*1e-6;  Hb = Sb.gB*Sb.Khat;  Pcb = make_Pc(Sb.e_hat, Pc_base);  % bias：離軸電荷（R150 fit 參數，兩半徑沿用）

    fixfig = fullfile(calroot,'current_base','figures','shared');
    nofig  = fullfile(nofix,'figures','shared');
    Clab  = '$\mathcal{C}\;[(\mathrm{mT/A})^{3}]$';
    NBINS = 135;                                                 % 統一 bin 數（粗細一致）
    NAVY  = [0.05 0.12 0.40];   BROWN = [0.55 0.35 0.17];        % 標準圖：R≤150 深藍、R≤500 棕
    OV_R150 = [0.45 0.68 0.90]; OV_R500 = [0.95 0.55 0.55];      % 疊圖專用：R_k 淺藍、R500 淺紅
    RLAB  = @(um) sprintf('R \\leq %d \\mum', um);               % legend 標籤（tex）

    %% ---- 一次算最大球（R500）逐點 C/κ；小半徑用 r 過濾（巢狀子集，等價 select_ball(R_k)）----
    [P, ~, npts] = select_ball(D, 500e-6);
    r  = vecnorm(P, 2, 2);                                        % 每點半徑 [m]
    Cf = zeros(npts,1); Kf = zeros(npts,1);  Cb = zeros(npts,1); Kb = zeros(npts,1);
    for i = 1:npts
        p = P(i,:).';
        svf = svd(((p/ellf - Pcf) ./ (vecnorm(p/ellf - Pcf).^3)) * Hf);
        svb = svd(((p/ellb - Pcb) ./ (vecnorm(p/ellb - Pcb).^3)) * Hb);
        Cf(i)=prod(svf); Kf(i)=svf(3)/svf(1);   % C=∏σ_k [(mT/A)³]、kappa=σ₃/σ₁（singular_value.pdf）
        Cb(i)=prod(svb); Kb(i)=svb(3)/svb(1);
    end
    i150 = r <= 150e-6;                                           % R150 子集（= select_ball(150e-6)）
    fprintf('\n=== R≤500µm 球 %d 節點；R≤150µm 子集 %d 節點 ===\n', npts, nnz(i150));
    fprintf('mean C : R150 fix=%.4g bias=%.4g / R500 fix=%.4g bias=%.4g (mT/A)^3\n', mean(Cf(i150)),mean(Cb(i150)),mean(Cf),mean(Cb));
    fprintf('mean k : R150 fix=%.4f bias=%.4f / R500 fix=%.4f bias=%.4f\n',           mean(Kf(i150)),mean(Kb(i150)),mean(Kf),mean(Kb));

    %% ---- 標準直方圖 8 張（R150 navy 子集 + R500 brown 全部；各自尺度）----
    render_hist(Cf(i150), Clab,       '(mT/A)^3', mean(Cf(i150)), fullfile(fixfig,'gain_index.png'),        NBINS, NAVY);
    render_hist(Kf(i150), '$\kappa$', '',         mean(Kf(i150)), fullfile(fixfig,'iso_index.png'),         NBINS, NAVY);
    render_hist(Cb(i150), Clab,       '(mT/A)^3', mean(Cb(i150)), fullfile(nofig, 'gain_index.png'),        NBINS, NAVY);
    render_hist(Kb(i150), '$\kappa$', '',         mean(Kb(i150)), fullfile(nofig, 'iso_index.png'),         NBINS, NAVY);
    render_hist(Cf,       Clab,       '(mT/A)^3', mean(Cf),       fullfile(fixfig,'gain_index_R500um.png'), NBINS, BROWN);
    render_hist(Kf,       '$\kappa$', '',         mean(Kf),       fullfile(fixfig,'iso_index_R500um.png'),  NBINS, BROWN);
    render_hist(Cb,       Clab,       '(mT/A)^3', mean(Cb),       fullfile(nofig, 'gain_index_R500um.png'), NBINS, BROWN);
    render_hist(Kb,       '$\kappa$', '',         mean(Kb),       fullfile(nofig, 'iso_index_R500um.png'),  NBINS, BROWN);

    %% ---- 掃描疊圖：R500 固定(棕) + R_k(navy) 掃 150:50:500；x 軸固定用 R500 尺度、legend 動態 ----
    %  每 dir gain/iso 各 8 張（_R{k}um）；x 軸恆為 nice_axis(R500) → 逐張看 R_k 分布長大。
    lab500 = RLAB(500);
    for Rk = 150:50:500
        ik = r <= Rk*1e-6;   labk = RLAB(Rk);   suf = sprintf('_R%dum', Rk);   % R_k 子集
        render_overlay(Cf(ik), Cf, Clab,       fullfile(fixfig,['gain_index_overlay' suf '.png']), OV_R150, OV_R500, NBINS, labk, lab500);
        render_overlay(Kf(ik), Kf, '$\kappa$', fullfile(fixfig,['iso_index_overlay'  suf '.png']), OV_R150, OV_R500, NBINS, labk, lab500);
        render_overlay(Cb(ik), Cb, Clab,       fullfile(nofig, ['gain_index_overlay' suf '.png']), OV_R150, OV_R500, NBINS, labk, lab500);
        render_overlay(Kb(ik), Kb, '$\kappa$', fullfile(nofig, ['iso_index_overlay'  suf '.png']), OV_R150, OV_R500, NBINS, labk, lab500);
    end

    %% ---- 無後綴 headline 疊圖（R150 vs R500）= R_k=150 frame，navy 重生 ----
    render_overlay(Cf(i150), Cf, Clab,       fullfile(fixfig,'gain_index_overlay.png'), OV_R150, OV_R500, NBINS, RLAB(150), lab500);
    render_overlay(Kf(i150), Kf, '$\kappa$', fullfile(fixfig,'iso_index_overlay.png'),  OV_R150, OV_R500, NBINS, RLAB(150), lab500);
    render_overlay(Cb(i150), Cb, Clab,       fullfile(nofig, 'gain_index_overlay.png'), OV_R150, OV_R500, NBINS, RLAB(150), lab500);
    render_overlay(Kb(i150), Kb, '$\kappa$', fullfile(nofig, 'iso_index_overlay.png'),  OV_R150, OV_R500, NBINS, RLAB(150), lab500);
end

function [lo, hi, tk] = nice_axis(vals)
% 從資料 union 算「乾淨」x 軸範圍 + 整數刻度：floor/ceil 到 1/2/5 nice step（不裁資料），刻度過多則隔一取一。
    dmin = min(vals);  dmax = max(vals);
    raw  = (dmax - dmin)/4;                    % 目標約 4 段
    p    = 10^floor(log10(raw));  r = raw/p;
    if     r >= 5, step = 5*p;
    elseif r >= 2, step = 2*p;
    else           step = p;
    end
    lo = floor(dmin/step)*step;   hi = ceil(dmax/step)*step;
    tk = lo:step:hi;
    if numel(tk) > 6, tk = tk(1:2:end); end    % 刻度太多 → 隔一取一（乾淨、避免邊緣擠刻度）
end

function render_hist(y, xlab, sunit, vline, out, nbins, facecolor, xrange, xticks)
% 真實直方圖（原始資料、零事後平滑）：長條高度 = histcounts 真實計數，不再 smoothdata。
%   nbins = 固定 bin 數（統一 → 長條根數/粗細一致）；nbins=[] → 預設 220。
%   facecolor = 長條填色（R≤150 藍 [0.20 0.40 0.70]、R≤500 棕 [0.55 0.35 0.17]）；預設藍。
%   xrange = [lo hi] 明確 x 軸範圍（R150/R500 共用、統一跨度）；xticks = 明確整數刻度（避免邊緣怪刻度）。皆選填。
%   長條寬度仍 = 各自資料 range/nbins（xrange 只改視窗，不改 binning）；左上黑框 N/mean/CV%、vline 紅虛線在 mean。
    n = numel(y);  m = mean(y);  s = std(y);
    if nargin < 6 || isempty(nbins), nbins = 220; end
    if nargin < 7 || isempty(facecolor), facecolor = [0.20 0.40 0.70]; end
    [cnt, edges] = histcounts(y, nbins);                          % 固定 bin 數（同根數 → 粗細一致）
    ctr = edges(1:end-1) + diff(edges)/2;
    fig = figure('Color','w','Position',[100 100 760 560]);
    ax  = axes(fig);
    bar(ax, ctr, cnt, 1, 'FaceColor',facecolor, 'EdgeColor','w', 'LineWidth',0.3); hold(ax,'on');   % 真實計數、白邊細長條
    if nargin >= 8 && ~isempty(xrange), xlo = xrange(1); xhi = xrange(2); else, [xlo, xhi, ~] = nice_axis(y); end   % 預設各自尺度
    if nargin >= 9 && ~isempty(xticks), xtk = xticks;                    else, [~, ~, xtk] = nice_axis(y); end       % 乾淨整數刻度（解邊緣怪刻度）
    xlim(ax, [xlo xhi]);
    xline(ax, vline, '--', 'Color',[0.85 0.20 0.20], 'LineWidth', 2);   % 標線（標在 mean）
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    set(ax,'XTick',xtk);
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xint = 'tex'; if startsWith(strtrim(xlab),'$'), xint = 'latex'; end
    xlabel(ax, xlab, 'FontWeight','bold','Interpreter',xint);
    ylabel(ax,'Count','FontWeight','bold');
    if isempty(sunit)
        txt = sprintf('N = %d\nmean = %.4g\nCV = %.2f%%', n, m, 100*s/m);        % CV = 100·std/mean（無因次，百分比）
    else
        txt = sprintf('N = %d\nmean = %.4g %s\nCV = %.2f%%', n, m, sunit, 100*s/m);
    end
    text(ax, 0.03, 0.95, txt, 'Units','normalized', 'HorizontalAlignment','left', ...
         'VerticalAlignment','top', 'FontSize',14, 'FontWeight','bold', ...
         'BackgroundColor','w', 'EdgeColor','k', 'LineWidth',1.5, 'Margin',5);
    ax.Toolbar.Visible = 'off';
    exportgraphics(fig, out, 'Resolution',600);
    fprintf('saved %s\n', out);
    close(fig);
end

function render_overlay(y1, y2, xlab, out, col1, col2, nbins, lab1, lab2)
% 疊圖：y1(col1 深藍, R_k) + y2(col2 棕, R500 固定) 疊在一起。x 軸統一用 y2(=R500) 尺度、兩者共用 bin edges（對齊）。
%   真實計數（半透明看重疊）；legend 用 lab1/lab2 標「哪個顏色 = 哪個半徑」；「僅此」= 無統計框、無均值線。
    if isempty(nbins), nbins = 135; end
    if nargin < 8 || isempty(lab1), lab1 = 'R \leq 150 \mum'; end
    if nargin < 9 || isempty(lab2), lab2 = 'R \leq 500 \mum'; end
    [rlo, rhi, rtk] = nice_axis(y2);                                   % x 軸尺度 = R500
    edges = linspace(rlo, rhi, nbins+1);                              % 共用 edges（R150/R500 對齊）
    fig = figure('Color','w','Position',[100 100 760 560]);
    ax  = axes(fig); hold(ax,'on');
    histogram(ax, y1, edges, 'FaceColor',col1, 'FaceAlpha',0.6, 'EdgeColor','w', 'LineWidth',0.3);   % R≤150（藍、中央尖峰）；白邊同標準圖
    histogram(ax, y2, edges, 'FaceColor',col2, 'FaceAlpha',0.6, 'EdgeColor','w', 'LineWidth',0.3);   % R≤500（棕、寬峰）
    xlim(ax, [rlo rhi]);  set(ax,'XTick',rtk);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xint = 'tex'; if startsWith(strtrim(xlab),'$'), xint = 'latex'; end
    xlabel(ax, xlab, 'FontWeight','bold','Interpreter',xint);
    ylabel(ax,'Count','FontWeight','bold');
    lg = legend(ax, {lab1, lab2}, 'Location','northeast', 'Interpreter','tex');
    set(lg,'FontSize',14,'FontWeight','bold','Box','on');
    ax.Toolbar.Visible = 'off';
    exportgraphics(fig, out, 'Resolution',600);
    fprintf('saved %s\n', out);
    close(fig);
end

% ---- local make_Pc (self-contained; charge lattice Pc = Pc_base + E(e17)) ----
function Pc = make_Pc(e17, Pc_base)
    E = zeros(3,6);
    E(:,1)=e17(1:3); E(:,2)=e17(4:6); E(:,3)=e17(7:9); E(:,4)=e17(10:12); E(:,5)=e17(13:15);
    E(1,6)=e17(16);  E(2,6)=e17(17);  E(3,6)=e17(1)-e17(4)+e17(8)-e17(11)+e17(15);
    Pc = Pc_base + E;
end
