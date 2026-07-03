function plot_gain_iso_index()
% PLOT_GAIN_ISO_INDEX  逐節點性能量 C / kappa 的**直方圖**（fix 與 bias 各兩張，共 4 張）。
%   R≤150µm 球內每個真實 FEM 節點 p 各自算致動增益矩陣 ^BG_I(p)=S(p)·Ĥ_I 的 SVD（σ₁≥σ₂≥σ₃），
%   依 singular_value.pdf 定義：
%     C(p)     = ∏σ_k = σ₁σ₂σ₃   （flux-generating ellipsoid 體積，(mT/A)³）
%     kappa(p) = σ_min/σ_max = σ₃/σ₁   （isotropy，≤1、→1 等向）
%   直方圖（風格同 err-hist 圖）：橫軸=C / kappa 值、縱軸=Count、左上黑框標 N/mean/std、紅虛線（C→mean、kappa→1）；NB=160 細長條、鋼藍白邊。
%   輸出：fix_dir/figures/{gain_index,iso_index}.png、no_fix_dir/figures/{gain_index,iso_index}.png。
%   ★ 兩模型都在 actuator 框同一組節點算（R_act·dhat=Pc_base）：fix 電荷=Pc_base（在軸）、bias=make_Pc(ê)（離軸）。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
    here  = fileparts(mfilename('fullpath'));
    nofix = fileparts(fileparts(here));
    calroot = fileparts(nofix);
    addpath(fullfile(nofix,'code','function'));
    addpath(fullfile(nofix,'code','main_function'));

    cnst = mt_constants();
    apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

    %% ---- R≤150µm 球內真實節點（actuator 框）----
    D = load_coils_actuator('long2016_hexapole_halfcut', cnst, apdl_to_paper_idx, 'all', 'gap200um_mueq');
    [P, ~, npts] = select_ball(D, 150e-6);
    Pc_base = D.Pc_base;

    Sf = load(fullfile(calroot,'fix_dir','data','fit_fixl_R150um_gap200um_mueq.mat'), 'ell','gB','Khat');
    Sb = load(fullfile(nofix,        'data','fit_bias_R150um_gap200um_mueq.mat'), 'ell','gB','Khat','e_hat');
    ellf = Sf.ell*1e-6;  Hf = Sf.gB*Sf.Khat;  Pcf = Pc_base;                     % fix：在軸電荷
    ellb = Sb.ell*1e-6;  Hb = Sb.gB*Sb.Khat;  Pcb = make_Pc(Sb.e_hat, Pc_base);  % bias：離軸電荷

    Cvf = zeros(npts,1); kapf = zeros(npts,1);  Cvb = zeros(npts,1); kapb = zeros(npts,1);
    for i = 1:npts
        p = P(i,:).';
        svf = svd(((p/ellf - Pcf) ./ (vecnorm(p/ellf - Pcf).^3)) * Hf);
        svb = svd(((p/ellb - Pcb) ./ (vecnorm(p/ellb - Pcb).^3)) * Hb);
        Cvf(i)=prod(svf); kapf(i)=svf(3)/svf(1);   % C=∏σ_k [(mT/A)³]、kappa=σ₃/σ₁（singular_value.pdf）
        Cvb(i)=prod(svb); kapb(i)=svb(3)/svb(1);
    end

    fixfig = fullfile(calroot,'fix_dir','figures');
    nofig  = fullfile(nofix,'figures');
    Clab = '$\mathcal{C}\;[(\mathrm{mT/A})^{3}]$';
    render_hist(Cvf,  Clab,       '(mT/A)^3', mean(Cvf), fullfile(fixfig,'gain_index.png'));   % C→mean 線
    render_hist(kapf, '$\kappa$', '',         1,         fullfile(fixfig,'iso_index.png'));    % kappa→1 線
    render_hist(Cvb,  Clab,       '(mT/A)^3', mean(Cvb), fullfile(nofig, 'gain_index.png'));
    render_hist(kapb, '$\kappa$', '',         1,         fullfile(nofig, 'iso_index.png'));

    fprintf('\n=== 逐節點 C / kappa（R≤150µm 球 %d 節點）===\n', npts);
    fprintf('mean C     : fix=%.4g  bias=%.4g (mT/A)^3\n', mean(Cvf), mean(Cvb));
    fprintf('mean kappa : fix=%.4f  bias=%.4f\n',           mean(kapf), mean(kapb));
end

function render_hist(y, xlab, sunit, vline, out)
% 平滑密度填充圖（連續、非離散長條）：細箱 histcounts + gaussian 平滑 → area 實心填充。
%   左上黑框註記、Count 縱軸、vline 紅虛線（C→mean、kappa→1）。
    n = numel(y);  m = mean(y);  s = std(y);
    NB = 500;                                                       % 高解析取樣（曲線用 500 點畫 → 平滑無折角）
    [cnt, edges] = histcounts(y, NB);
    ctr   = edges(1:end-1) + diff(edges)/2;
    cnt_s = smoothdata(cnt, 'gaussian', round(NB/6));               % 強高斯平滑：徹底去 bin 噪聲、看不出離散
    fig = figure('Color','w','Position',[100 100 760 560]);
    ax  = axes(fig);
    area(ax, ctr, cnt_s, 'FaceColor',[0.20 0.40 0.70], 'EdgeColor','none'); hold(ax,'on');   % 實心、連續
    xl = xlim(ax);                                                  % 若 vline 貼/超右界 → 延伸右界讓虛線落在框內
    if vline >= xl(2) - 0.005*(xl(2)-xl(1))
        xlim(ax, [xl(1), vline + 0.03*(xl(2)-xl(1))]);
    end
    xline(ax, vline, '--', 'Color',[0.85 0.20 0.20], 'LineWidth', 2);   % 標線（C→mean、kappa→1）
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on'); grid(ax,'off');
    xt = get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xint = 'tex'; if startsWith(strtrim(xlab),'$'), xint = 'latex'; end
    xlabel(ax, xlab, 'FontWeight','bold','Interpreter',xint);
    ylabel(ax,'Count','FontWeight','bold');
    if isempty(sunit)
        txt = sprintf('N = %d\nmean = %.4g\nstd = %.4g', n, m, s);
    else
        txt = sprintf('N = %d\nmean = %.4g %s\nstd = %.4g %s', n, m, sunit, s, sunit);
    end
    text(ax, 0.03, 0.95, txt, 'Units','normalized', 'HorizontalAlignment','left', ...
         'VerticalAlignment','top', 'FontSize',14, 'FontWeight','bold', ...
         'BackgroundColor','w', 'EdgeColor','k', 'LineWidth',1.5, 'Margin',5);
    ax.Toolbar.Visible = 'off';
    exportgraphics(fig, out, 'Resolution',600);
    fprintf('saved %s\n', out);
    close(fig);
end
