%% run_hung_current.m -- hung R700 + R300 current(bias) 校正 + gain/iso 疊圖
% =========================================================================
%  用共用管線（extract_ansys_data / build_D / fitting / solve_current / emit_results）對 hung 的
%  兩個工作半徑變體 R700 / R300 做 current 側 18-param(bias) 校正：
%    ① 校正結果 PDF → Calibration_using_FEM_modeling/results/hung_hexapole/eighteen/
%    ② R≤150µm 的逐節點 gain 𝒞 / iso κ（均勻 3D 網格連續取樣，非 FEM 節點；避免密度偏差）
%       → R300 vs R700 疊圖（style 同 plot_gain_iso_overlay：nb=180、共用 edges、mean 虛線、legend 含 mean）
%       → hung_hexapole/Calibration_using_FEM_modeling/current_base/figures/shared/
%  𝒞=∏σ、κ=σ₃/σ₁，σ=svd(S·Ĥ_I)，Ĥ_I=ĝ_I·K̄_I、S 用 biased 電荷格 Pc=Pc_base+E(e)（與 solve_current 一致）。
% =========================================================================

clear; clc;

%% ---- per-run 調參 --------------------------------------------------------
USE_BIAS = true;               % 18-param bias → results/.../eighteen
R_select = 150e-6;             % fit 取點球半徑 [m]
R_iso    = 150e-6;             % gain/iso 取樣半徑 [m]（R≤150µm）
l0       = 0.5e-3;             % l_hat 初值 [m]
I_actual = 1;                  % 驅動電流 [A] = FEM 激發
MODEL    = 'hung_hexapole';
VARIANTS = {'R700','R300'};    % 兩個工作半徑變體（共用 R_norm=500µm 幾何、WP-centered）

%% ---- paths（相對自身定位）-----------------------------------------------
here  = fileparts(mfilename('fullpath'));                    % utils/hung_hexapole
CAL   = fileparts(fileparts(here));                          % Calibration_using_FEM_modeling
APDL  = fileparts(CAL);                                      % .../main/matlab/APDL
ANSYS = fileparts(fileparts(fileparts(APDL)));               % .../ANSYS
addpath(fullfile(CAL, 'function'));
addpath(fullfile(CAL, 'common_path'));                                  % ansys_path（共用路徑 resolver）
% import_ansys_data / filter_iron_nodes 已入 function/（shared 自足、不再依賴 backup）
figdir = fullfile(CAL, 'figures', 'hung_hexapole', 'current', 'common');   % 共用夾 figures（R300/R700 疊圖屬 common）
if ~exist(figdir,'dir'), mkdir(figdir); end

cfg = model_config(MODEL);

%% ---- 校正兩個變體 + 收集 gain/iso 所需（ell / Ĥ_I / Pc）------------------
S = struct('variant',{},'ell',{},'Hhat',{},'Pc',{});
for c = 1:numel(VARIANTS)
    variant = VARIANTS{c};
    raw = extract_ansys_data(cfg, 'all', variant);
    D   = build_D(raw, cfg);
    [P, Bstack, npts] = cfg.select_ball(D, R_select);

    [e, l_hat, J] = fitting(P, Bstack, D.Pc_base, l0, USE_BIAS);
    F = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
    [KI_bar, gI_hat, G, rm] = solve_current(l_hat, e, D.Pc_base, P, Bstack, F);

    rec = struct('base','current', 'l_hat',l_hat, 'e',e, 'J',J, 'npts',npts, ...
                 'model',MODEL, 'VARIANT',variant, 'DATASET','all', 'USE_BIAS',USE_BIAS, ...
                 'R_select',R_select, 'l0',l0, 'I_actual',I_actual, ...
                 'apdl_to_paper_idx',cfg.apdl_to_paper_idx, 's_source',cfg.s_source, ...
                 'KI_bar',KI_bar, 'gI_hat',gI_hat, 'G',G, 'Fmap',F, ...
                 'C_mean',rm.C_mean, 'kappa_mean',rm.kappa_mean, 'C_min',rm.C_min, 'kappa_worst',rm.kappa_worst);
    matdir = fullfile(CAL, 'data', MODEL, '.mat');
    if ~exist(matdir,'dir'), mkdir(matdir); end
    matfile = fullfile(matdir, sprintf('calib_current_%s_R%03d_eighteen.mat', variant, round(R_select*1e6)));
    save(matfile, '-struct', 'rec');

    errpct = 100*sqrt(J / sum(Bstack(:).^2));
    fprintf('[%s] l_hat=%.1f um | gI=%.4f mT/A | K(1,1)=%.4f | diag=[%s] | eps=%.3f%% | npts=%d\n', ...
            variant, l_hat*1e6, gI_hat, KI_bar(1,1), sprintf('%.3f ',diag(KI_bar)), errpct, npts);
    emit_results(matfile);

    Pc = make_Pc_local(e, D.Pc_base);
    S(c) = struct('variant',variant, 'ell',l_hat, 'Hhat',gI_hat*KI_bar, 'Pc',Pc);
end

%% ---- gain 𝒞 / iso κ：均勻球內網格連續取樣 → R300 vs R700 疊圖 -------------
Pg = ball_grid(R_iso, 32);                                   % 均勻 3D 網格（WP 原點）
[C7,K7] = ck_nodes(Pg, S(1).ell, S(1).Hhat, S(1).Pc);       % R700
[C3,K3] = ck_nodes(Pg, S(2).ell, S(2).Hhat, S(2).Pc);       % R300
fprintf('gain/iso (N=%d, R<=%gµm)  R700: meanC=%.4g meanK=%.4f | R300: meanC=%.4g meanK=%.4f\n', ...
        size(Pg,1), R_iso*1e6, mean(C7), mean(K7), mean(C3), mean(K3));

RED = [0.60 0.08 0.08];  BLUE = [0.12 0.18 0.55];
% R300 深紅 / R700 深藍（半透明、共用 bins、mean 虛線、legend 含 mean）
% gain 用 log x（R300/R700 差 ~15×、線性軸下 R700 擠成尖峰）；iso 線性
render_overlay(C3, C7, '$\mathcal{C}\;[(\mathrm{mT/A})^{3}]$', RED, BLUE, figdir, ...
               sprintf('gain_R300_vs_R700_R%dum.png', round(R_iso*1e6)), {'R300','R700'}, true);
render_overlay(K3, K7, '$\kappa$', RED, BLUE, figdir, ...
               sprintf('iso_R300_vs_R700_R%dum.png', round(R_iso*1e6)), {'R300','R700'}, false);

%% ======================= local functions ==================================
function Pc = make_Pc_local(e17, Pc_base)
    E = zeros(3,6);
    E(:,1)=e17(1:3); E(:,2)=e17(4:6); E(:,3)=e17(7:9); E(:,4)=e17(10:12); E(:,5)=e17(13:15);
    E(1,6)=e17(16);  E(2,6)=e17(17);  E(3,6)=e17(1)-e17(4)+e17(8)-e17(11)+e17(15);
    Pc = Pc_base + E;
end

function P = ball_grid(R, nr)
    h = R/nr;  v = -R:h:R;
    [X,Y,Z] = ndgrid(v,v,v);
    in = (X.^2 + Y.^2 + Z.^2) <= R^2;
    P = [X(in), Y(in), Z(in)];
end

function [Cv,Kv] = ck_nodes(P, ell_m, Hhat, Pc)
    n = size(P,1);  Cv = zeros(n,1);  Kv = zeros(n,1);
    for i = 1:n
        Dk = P(i,:).'/ell_m - Pc;                     % 3×6
        sv = svd((Dk ./ (vecnorm(Dk).^3)) * Hhat);
        Cv(i) = prod(sv);  Kv(i) = sv(3)/sv(1);
    end
end

function render_overlay(yA, yB, xlab, colA, colB, figdir, fname, labels, logx)
    if nargin < 9 || isempty(logx), logx = false; end
    nb = 180;
    if logx
        edges = logspace(log10(min([yA;yB])), log10(max([yA;yB])), nb+1);   % log 等距 bins
    else
        edges = linspace(min([yA;yB]), max([yA;yB]), nb+1);
    end
    pA = histcounts(yA, edges) / numel(yA) * 100;      % 百分比（各組 count/N×100）
    pB = histcounts(yB, edges) / numel(yB) * 100;
    fig = figure('Color','w','Position',[100 100 820 600]);
    ax  = axes(fig); hold(ax,'on');
    hA = histogram(ax, 'BinEdges',edges, 'BinCounts',pA, 'FaceColor',colA, 'FaceAlpha',0.55, 'EdgeColor','w', 'LineWidth',0.3);
    hB = histogram(ax, 'BinEdges',edges, 'BinCounts',pB, 'FaceColor',colB, 'FaceAlpha',0.55, 'EdgeColor','w', 'LineWidth',0.3);
    xline(ax, mean(yA), '--', 'Color',colA, 'LineWidth',2);
    xline(ax, mean(yB), '--', 'Color',colB, 'LineWidth',2);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    if logx, set(ax,'XScale','log'); end
    box(ax,'on'); grid(ax,'off');
    yl = ylim(ax); ylim(ax,[0 yl(2)*1.20]);
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xint = 'tex'; if startsWith(strtrim(xlab),'$'), xint = 'latex'; end
    xlabel(ax, xlab, 'FontWeight','bold','Interpreter',xint);
    ylabel(ax,'Percentage (%)','FontWeight','bold');
    lg = legend(ax, [hA hB], {sprintf('%s  (mean = %.4g)', labels{1}, mean(yA)), ...
                     sprintf('%s  (mean = %.4g)', labels{2}, mean(yB))}, 'Location','northeast');
    set(lg,'FontSize',15,'FontWeight','bold','Box','on');
    ax.Toolbar.Visible = 'off';
    exportgraphics(fig, fullfile(figdir, fname), 'Resolution',600);
    fprintf('saved %s\n', fullfile(figdir, fname));
    close(fig);
end
