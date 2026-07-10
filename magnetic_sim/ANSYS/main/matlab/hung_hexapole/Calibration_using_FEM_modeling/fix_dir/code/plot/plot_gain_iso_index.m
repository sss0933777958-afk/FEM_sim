function plot_gain_iso_index()
% PLOT_GAIN_ISO_INDEX  hung 逐節點性能量 C / kappa 直方圖（fix-l charge model，R500 球）。
%   球內每個真實 FEM 節點 p 各自算致動增益矩陣 ^BG_I(p)=S(p)·Ĥ_I 的 SVD（σ₁≥σ₂≥σ₃）：
%     C(p)     = ∏σ_k = σ₁σ₂σ₃   （(mT/A)³）
%     kappa(p) = σ₃/σ₁            （isotropy，≤1、→1 等向）
%   charge model：S(p)=(p/ℓ̂ − Pc_base)./|p/ℓ̂ − Pc_base|³（Coulomb kernel），Ĥ_I=gB·K̂。
%   在 actuator 框同一組真實節點上算（Pc_base=R_act·d̂，在軸單位晶格）。
%   資料：hung 6-coil FEM（gap200um_mueq，1 A、all-source flip-sink、actuator frame）。
%   輸出：fix_dir/figures/{gain,iso}_index_R500um.png（棕、135 bins、N/mean/CV 框、紅虛線=mean）。

    here   = fileparts(mfilename('fullpath'));
    fixdir = fileparts(fileparts(here));                 % .../fix_dir
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis'); % filter_iron_nodes
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hung\analysis\core');          % hung mt_constants + import_ansys_data
    addpath(fullfile(fixdir,'code','function'));
    addpath(fullfile(fixdir,'code','main_function'));

    cnst = mt_constants();                                % hung
    idx  = cnst.apdl_to_paper_idx;                        % [1,2,3,4,5,6]
    results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\hung_hexapole\data';

    %% ---- 載 6 coil（air 節點、actuator frame）+ 固定校正參數（R150 fit）----
    C = load_coils(results_root, cnst, idx, 'gap_200um');
    P = C(1).P;  r = vecnorm(P,2,2);  P = P(r<=500e-6, :);   % R500 球內真實節點
    npts = size(P,1);

    tip   = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    dhat  = tip ./ vecnorm(tip);
    R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';
    Pc_base = R_act * dhat;                                % actuator 框在軸電荷晶格

    S = load(fullfile(fixdir,'data','fit_fixl_R150um_gap200um_mueq.mat'), 'ell','gB','Khat');
    ellf = S.ell*1e-6;  Hf = S.gB*S.Khat;

    %% ---- 逐節點 C / kappa ----
    Cf = zeros(npts,1);  Kf = zeros(npts,1);
    for i = 1:npts
        p  = P(i,:).';
        sv = svd(((p/ellf - Pc_base) ./ (vecnorm(p/ellf - Pc_base).^3)) * Hf);
        Cf(i) = prod(sv);  Kf(i) = sv(3)/sv(1);
    end
    fprintf('R<=500um 球 %d 節點 | mean C=%.4g (mT/A)^3 | mean kappa=%.4f\n', npts, mean(Cf), mean(Kf));

    %% ---- 直方圖 ----
    figdir = fullfile(fixdir,'figures');  if ~exist(figdir,'dir'); mkdir(figdir); end
    Clab  = '$\mathcal{C}\;[(\mathrm{mT/A})^{3}]$';
    NBINS = 135;  BROWN = [0.55 0.35 0.17];
    render_hist(Cf, Clab,       '(mT/A)^3', mean(Cf), fullfile(figdir,'gain_index_R500um.png'), NBINS, BROWN);
    render_hist(Kf, '$\kappa$', '',         mean(Kf), fullfile(figdir,'iso_index_R500um.png'),  NBINS, BROWN);
end

function [lo, hi, tk] = nice_axis(vals)
% 從資料算「乾淨」x 軸範圍 + 整數刻度：floor/ceil 到 1/2/5 nice step；刻度過多隔一取一。
    dmin = min(vals);  dmax = max(vals);
    raw  = (dmax - dmin)/4;
    p    = 10^floor(log10(raw));  rr = raw/p;
    if     rr >= 5, step = 5*p;
    elseif rr >= 2, step = 2*p;
    else            step = p;
    end
    lo = floor(dmin/step)*step;   hi = ceil(dmax/step)*step;
    tk = lo:step:hi;
    if numel(tk) > 6, tk = tk(1:2:end); end
end

function render_hist(y, xlab, sunit, vline, out, nbins, facecolor)
% 真實直方圖（原始資料、零事後平滑）：長條高度 = histcounts 真實計數。
    n = numel(y);  m = mean(y);  s = std(y);
    if nargin < 6 || isempty(nbins), nbins = 220; end
    if nargin < 7 || isempty(facecolor), facecolor = [0.20 0.40 0.70]; end
    [cnt, edges] = histcounts(y, nbins);
    ctr = edges(1:end-1) + diff(edges)/2;
    fig = figure('Color','w','Position',[100 100 760 560]);
    ax  = axes(fig);
    bar(ax, ctr, cnt, 1, 'FaceColor',facecolor, 'EdgeColor','w', 'LineWidth',0.3); hold(ax,'on');
    [xlo, xhi, xtk] = nice_axis(y);
    xlim(ax, [xlo xhi]);
    xline(ax, vline, '--', 'Color',[0.85 0.20 0.20], 'LineWidth', 2);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    set(ax,'XTick',xtk);
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xint = 'tex'; if startsWith(strtrim(xlab),'$'), xint = 'latex'; end
    xlabel(ax, xlab, 'FontWeight','bold','Interpreter',xint);
    ylabel(ax,'Count','FontWeight','bold');
    if isempty(sunit)
        txt = sprintf('N = %d\nmean = %.4g\nCV = %.2f%%', n, m, 100*s/m);
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
