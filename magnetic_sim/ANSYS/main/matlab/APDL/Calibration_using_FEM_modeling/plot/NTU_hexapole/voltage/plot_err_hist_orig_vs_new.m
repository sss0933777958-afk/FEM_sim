function plot_err_hist_orig_vs_new(variant)
%PLOT_ERR_HIST_ORIG_VS_NEW  用「原 sensor 位置」與「新 sensor 位置(找到的 50% 點)」各自校正，
%   出場誤差 |B_model − B_FEM| 直方圖疊圖 → figures/eighteen_param/error_hist_orig_vs_new.png。
%   兩位置：①原標稱 sensor、②vmat_search 找到的 50% 新位置。各自 build V̄ → Ĥ_V=G·V̄⁻¹ →
%   從該位置 V̄ 反推電荷 Ĝ=Ĥ_V·V̄ → 重建場 A·Ĝ → 與 FEM Bstack 逐點逐激發比。18-param bias。
    if nargin<1 || isempty(variant), variant = 'full_assembly'; end

    CAL = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
           'NTU_hexapole\Calibration_using_FEM_modeling'];
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live 樹。
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
    addpath(fullfile(CAL,'current_base','code','main_function'));
    addpath(fullfile(CAL,'voltage_base','code','main_function'));
    results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';
    here = fileparts(mfilename('fullpath'));  base = fileparts(fileparts(here));

    S = load(fullfile(base,'data',sprintf('calib_D_%s.mat',variant)));   % R_select/S_hall/n_uniform/apdl/paper
    apdl_to_paper = [1,3,6,5,2,4];
    cnst = model_config('NTU_hexapole');

    %% ---- 電荷擬合（WP 球，bias）→ A, Dv(G), Bstack ----
    D = load_coils_actuator('NTU_hexapole', cnst, apdl_to_paper, 'all', variant);
    [P, Bstack, npts] = select_ball(D, S.R_select);
    s_sink = ones(1,6);  for j=1:6, if ismember(apdl_to_paper(j),[1 3 6]), s_sink(j)=-1; end; end
    Bstack = (-Bstack).*s_sink;                                          % all-source [mT]
    [ell, ~, Pc, ~] = fitting(P, Bstack, D.Pc_base, 0.5e-3, true);
    A  = build_S_matrix(ell, Pc, P);   Dv = (A.'*A)\(A.'*Bstack);        % 3Np×6, 6×6

    %% ---- 兩位置 sensor 幾何 → 各自 V̄ ----
    [sp, sn] = build_sensor_geometry(cnst);                             % 原標稱位置
    found_mm = [ 11.94  1.59  9.46;    % P1  找到的 50% 新位置（world mm，取自 vmat_search PDF）
                -14.05  2.55 10.21;    % P2
                 -6.72  8.94  9.46;    % P3
                  4.39 -11.15 10.21;   % P4
                  4.76  12.44 10.21;   % P5
                 -4.68 -11.85  9.46];  % P6
    found_pos = found_mm.'*1e-3;                                         % 3×6 [m]
    Vorig = build_V_matrix(results_root, cnst, apdl_to_paper, sp,        sn, S.S_hall, [], S.n_uniform, [], 0.10e-3, variant);
    Vnew  = build_V_matrix(results_root, cnst, apdl_to_paper, found_pos, sn, S.S_hall, [], S.n_uniform, [], 0.10e-3, variant);

    %% ---- 校正 Ĥ_V=G·V̄⁻¹（用原位置 V̄，固定）→ 兩位置 V̄ 各自反推電荷 → 重建場 → 誤差 ----
    Hv   = Dv / Vorig;                                                   % Ĥ_V 校正於原 sensor 位置（一次）
    erro = vecnorm(reshape(A*(Hv*Vorig) - Bstack, 3, npts*6), 2, 1).';   % 讀原位置 V̄_pdf → Ĝ=G（正確）
    errn = vecnorm(reshape(A*(Hv*Vnew ) - Bstack, 3, npts*6), 2, 1).';   % 讀新位置 V̄_found → Ĝ=Ĥ_V·V̄_found≠G（sensor 移位，重建偏掉）
    fprintf('orig: median=%.4f mT | new: median=%.4f mT | max|Vorig-Vnew|=%.3f mV\n', ...
            median(erro), median(errn), max(abs(Vorig(:)-Vnew(:))));

    %% ---- 疊圖 → figures/eighteen_param/ ----
    figdir = fullfile(base,'figures','eighteen_param');  if ~exist(figdir,'dir'); mkdir(figdir); end
    nb = 180;  XMAX = prctile([erro;errn], 99.5);  edges = linspace(0, XMAX, nb+1);
    colO = [0.85 0.45 0.30];  colN = [0.35 0.50 0.80];
    fig = figure('Color','w','Position',[100 100 900 620]);  ax = axes(fig);  hold(ax,'on');
    hO = histogram(ax, erro, edges, 'FaceColor',colO, 'FaceAlpha',0.55, 'EdgeColor','w','LineWidth',0.3);
    hN = histogram(ax, errn, edges, 'FaceColor',colN, 'FaceAlpha',0.55, 'EdgeColor','w','LineWidth',0.3);
    xline(ax, median(erro), '--', 'Color',colO, 'LineWidth',2);
    xline(ax, median(errn), '--', 'Color',colN, 'LineWidth',2);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);  box(ax,'on'); grid(ax,'off');
    xlim(ax,[0 XMAX]);  yl=ylim(ax); ylim(ax,[0 yl(2)*1.20]);
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));  yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel(ax,'|B_{model} - B_{FEM}| (mT)','FontWeight','bold');  ylabel(ax,'Count','FontWeight','bold');
    lg = legend(ax,[hO hN], {sprintf('original position  (median %.3f mT)', median(erro)), ...
                             sprintf('new position  (median %.3f mT)', median(errn))}, ...
                'Location','northeast','Interpreter','tex');  set(lg,'FontSize',14,'FontWeight','bold','Box','on');
    outpng = fullfile(figdir,'error_hist_orig_vs_new.png');
    exportgraphics(fig, outpng, 'Resolution',300);  close(fig);
    fprintf('saved %s\n', outpng);
end
