function plot_gain_iso_overlay(variant, variant_hung, labels, out_tag)
% PLOT_GAIN_ISO_OVERLAY  hung vs long 的 gain(C)/iso(κ) 疊圖（cross-model 對照）。
%   plot_gain_iso_overlay()               -> gap 版（'gap_200um'，向後相容）
%   plot_gain_iso_overlay('no_gap') -> no-gap 版（檔名加 '_nogap'）
%   plot_gain_iso_overlay('gap_200um','gap_200um') -> [ADDED] long/hung 吃不同 variant：
%       long=第1參數(fit+lattice)、hung=第2參數(僅 fit)。檔名 suffix 依「第1參數」(long) 決定，
%       故 ('gap_200um','gap_200um') 仍覆蓋 canonical gain_hung_vs_long_R*.png（hung=新 100mm² gap）。
    if nargin < 1, variant = 'gap_200um'; end                 % [ADDED] variant 切換（long）
    if nargin < 2 || isempty(variant_hung), variant_hung = variant; end  % [ADDED] hung variant（預設 = long）
    if nargin < 3 || isempty(labels),  labels  = {'long','hung'}; end   % [ADDED] legend 名稱（long/hung → 可換 design1/design2）
    if nargin < 4 || isempty(out_tag), out_tag = ''; end                % [ADDED] 檔名後綴（另存新檔用，如 '_design'）
    if     strcmp(variant,'gap_200um'),  vsfx = '';           % gap 保留原檔名
    elseif strcmp(variant,'no_gap'),   vsfx = '_nogap';
    else,                                    vsfx = ['_' variant]; end
%   C、κ 是位置相關的「平滑擬合函數」（charge model：C=∏σ_k、κ=σ₃/σ₁；非 raw FEM 場）→
%   **在球內「均勻 3D 網格」上連續取樣**（不用 FEM 節點；節點密度不均，尤其 WP 核心球會
%   過度加權中心 → 直方圖有偏）。均勻體積取樣 = 無偏的 C/κ 空間分布。
%   long=深紅 / hung=深藍（半透明、共用 bins、mean 虛線、legend headroom 不壓圖）。
%   同 gap200um_mueq、70 匝，R150 fit 固定；取樣半徑 R∈{150,500}µm。
%   出 4 檔：{gain,iso}_hung_vs_long_R{150,500}um.png，**同存 long 與 hung 兩個 figures/**。

    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ 改指 live config；
    %   本檔靠「後 addpath 覆蓋」同時取得 long 與 hung 兩套 mt_constants，該語意保持不變。
    % ⚠ 本檔仍有**其他**死相依（LONG_CAL/HUNG_CAL 指向已刪除的舊 per-model 樹、
    %   load_coils_actuator 已不存在）→ 目前跑不起來，需另外重寫才能復活。
    CALROOT       = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    LONG_ANALYSIS = fullfile(CALROOT,'config','long2016_hexapole_halfcut','tip40um');
    HUNG_CORE     = fullfile(CALROOT,'config','hung_hexapole');
    LONG_CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\hung_hexapole\Calibration_using_FEM_modeling';
    HUNG_CAL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\hung_hexapole\Calibration_using_FEM_modeling';

    %% ---- long：只取電荷 lattice Pc_base + fit（不需 FEM 節點）----
    addpath(LONG_ANALYSIS);
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');  % ansys_path resolver
    addpath(fullfile(LONG_CAL,'current_base','code','main_function'));      % load_coils_actuator
    cnstL = mt_constants();                                               % long
    DL = load_coils_actuator('hung_hexapole', cnstL, [1,3,6,5,2,4], 'all', variant);
    PcL = DL.Pc_base;
    SL = load(fullfile(LONG_CAL,'current_base','data',sprintf('fit_fixl_R150um_%s.mat',variant)), 'ell','gB','Khat');
    ellL = SL.ell*1e-6;  HhatL = SL.gB*SL.Khat;

    %% ---- hung：electric lattice PcH + fit（後 addpath hung core 蓋掉 long）----
    addpath(HUNG_CORE);
    addpath(fullfile(HUNG_CAL,'current_base','code','main_function'));
    cnstH = mt_constants();                                               % hung（後 addpath 優先）
    tip = [cnstH.pole_tip_x; cnstH.pole_tip_y; cnstH.pole_tip_z_wp];  dhat = tip./vecnorm(tip);
    R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';  PcH = R_act*dhat;
    SH = load(fullfile(HUNG_CAL,'current_base','data',sprintf('fit_fixl_R150um_%s.mat',variant_hung)), 'ell','gB','Khat');  % [MODIFIED] hung 用 variant_hung
    ellH = SH.ell*1e-6;  HhatH = SH.gB*SH.Khat;

    %% ---- figures dirs, colors ----
    figL = fullfile(LONG_CAL,'current_base','figures','shared');
    figH = fullfile(HUNG_CAL,'current_base','figures','shared');
    RED  = [0.60 0.08 0.08];   BLUE = [0.12 0.18 0.55];

    %% ---- 迴圈 R=150,500µm：球內均勻網格連續取樣 → 疊圖 ----
    for R_um = [150 500]
        R = R_um*1e-6;
        P = ball_grid(R, 32);                                            % 均勻 3D 網格（WP 在原點）
        [CL,KL] = ck_nodes(P, ellL, HhatL, PcL);                         % long
        [CH,KH] = ck_nodes(P, ellH, HhatH, PcH);                         % hung
        fprintf('R%d (N=%d)  long: meanC=%.4g meanK=%.4f | hung: meanC=%.4g meanK=%.4f | hung/long C=%.4f\n', ...
                R_um, size(P,1), mean(CL), mean(KL), mean(CH), mean(KH), mean(CH)/mean(CL));
        render_overlay(CL, CH, '$\mathcal{C}\;[(\mathrm{mT/A})^{3}]$', RED, BLUE, {figL,figH}, sprintf('gain_hung_vs_long%s_R%dum%s.png', vsfx, R_um, out_tag), labels);
        render_overlay(KL, KH, '$\kappa$',                            RED, BLUE, {figL,figH}, sprintf('iso_hung_vs_long%s_R%dum%s.png', vsfx, R_um, out_tag), labels);
    end
end

function P = ball_grid(R, nr)
% 均勻 Cartesian 網格取「半徑 R 球內」點（nr = 半徑方向格數；點數 ~ (π/6)(2nr+1)³）。
    h = R/nr;  v = -R:h:R;
    [X,Y,Z] = ndgrid(v,v,v);
    in = (X.^2 + Y.^2 + Z.^2) <= R^2;
    P = [X(in), Y(in), Z(in)];
end

function [Cv,Kv] = ck_nodes(P, ell_m, Hhat, Pc)
    n = size(P,1);  Cv = zeros(n,1);  Kv = zeros(n,1);
    for i = 1:n
        p  = P(i,:).';
        Dk = p/ell_m - Pc;                                    % 3×6
        sv = svd((Dk ./ (vecnorm(Dk).^3)) * Hhat);
        Cv(i) = prod(sv);  Kv(i) = sv(3)/sv(1);
    end
end

function render_overlay(yL, yH, xlab, colL, colH, figdirs, fname, labels)
    if nargin < 8 || isempty(labels), labels = {'long','hung'}; end   % [ADDED] legend 名稱
    nb = 180;                                              % 離散直方圖、bin 夠密
    edges = linspace(min([yL;yH]), max([yL;yH]), nb+1);
    fig = figure('Color','w','Position',[100 100 820 600]);
    ax  = axes(fig); hold(ax,'on');
    hL = histogram(ax, yL, edges, 'FaceColor',colL, 'FaceAlpha',0.55, 'EdgeColor','w', 'LineWidth',0.3);   % long 深紅
    hH = histogram(ax, yH, edges, 'FaceColor',colH, 'FaceAlpha',0.55, 'EdgeColor','w', 'LineWidth',0.3);   % hung 深藍
    xline(ax, mean(yL), '--', 'Color',colL, 'LineWidth',2);                             % 兩條均值虛線（同色）
    xline(ax, mean(yH), '--', 'Color',colH, 'LineWidth',2);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    yl = ylim(ax); ylim(ax,[0 yl(2)*1.20]);   % headroom：留白讓 legend 不壓到 bar
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xint = 'tex'; if startsWith(strtrim(xlab),'$'), xint = 'latex'; end
    xlabel(ax, xlab, 'FontWeight','bold','Interpreter',xint);
    ylabel(ax,'Count','FontWeight','bold');
    lg = legend(ax, [hL hH], {sprintf('%s  (mean = %.4g)', labels{1}, mean(yL)), ...
                     sprintf('%s  (mean = %.4g)', labels{2}, mean(yH))}, 'Location','northeast');   % legend 含 mean
    set(lg,'FontSize',15,'FontWeight','bold','Box','on');
    ax.Toolbar.Visible = 'off';
    for k = 1:numel(figdirs)
        if ~exist(figdirs{k},'dir'); mkdir(figdirs{k}); end
        exportgraphics(fig, fullfile(figdirs{k}, fname), 'Resolution',600);
        fprintf('saved %s\n', fullfile(figdirs{k}, fname));
    end
    close(fig);
end
