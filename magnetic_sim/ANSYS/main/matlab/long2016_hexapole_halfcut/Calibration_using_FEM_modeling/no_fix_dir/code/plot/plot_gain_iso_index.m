function plot_gain_iso_index()
% PLOT_GAIN_ISO_INDEX  逐節點 gain / iso vs node index（fix 與 bias 各兩張，共 4 張）。
%   R≤150µm 球內每個真實 FEM 節點 p 各自算 gain(p)=‖T‖_F、iso(p)=σ_max/σ_min，T=S(p)·Ĥ_I。
%   橫軸=node index、縱軸=gain / iso；圖上標註該量的**變異數 Var**（+ mean 參考線/值）。
%   輸出：fix_dir/figures/{gain_index,iso_index}.png、no_fix_dir/figures/{gain_index,iso_index}.png。
%   ★ 兩模型都在 actuator 框同一組節點算（R_act·dhat=Pc_base）：fix 電荷=Pc_base（在軸）、bias=make_Pc(ê)（離軸）。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
    here  = fileparts(mfilename('fullpath'));
    nofix = fileparts(fileparts(here));
    calroot = fileparts(nofix);
    addpath(fullfile(nofix,'code','function'));

    cnst = mt_constants();
    apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];

    %% ---- R≤150µm 球內真實節點（actuator 框）----
    D = load_coils_actuator('long2016_hexapole_halfcut', cnst, apdl_to_paper_idx, 'all', 'gap200um_mueq');
    [P, ~, npts] = select_ball(D, 150e-6);
    Pc_base = D.Pc_base;

    Sf = load(fullfile(calroot,'fix_dir','data','fit_fixl_R150um_gap200um_mueq.mat'), 'ell','gB','Khat');
    Sb = load(fullfile(nofix,        'data','fit_bias_R150um_gap200um_mueq.mat'), 'ell','gB','Khat','e_hat');
    ellf = Sf.ell*1e-6;  Hf = Sf.gB*Sf.Khat;  Cf = Pc_base;
    ellb = Sb.ell*1e-6;  Hb = Sb.gB*Sb.Khat;  Cb = make_Pc(Sb.e_hat, Pc_base);

    gf = zeros(npts,1); isf = zeros(npts,1);  gb = zeros(npts,1); isb = zeros(npts,1);
    for i = 1:npts
        p = P(i,:).';
        svf = svd(((p/ellf - Cf) ./ (vecnorm(p/ellf - Cf).^3)) * Hf);
        svb = svd(((p/ellb - Cb) ./ (vecnorm(p/ellb - Cb).^3)) * Hb);
        gf(i)=norm(svf); isf(i)=svf(1)/svf(3);
        gb(i)=norm(svb); isb(i)=svb(1)/svb(3);
    end

    fixfig = fullfile(calroot,'fix_dir','figures');
    nofig  = fullfile(nofix,'figures');
    render_index(gf,  'Gain (mT/A)', '(mT/A)^2', fullfile(fixfig,'gain_index.png'));
    render_index(isf, 'iso',         '',         fullfile(fixfig,'iso_index.png'));
    render_index(gb,  'Gain (mT/A)', '(mT/A)^2', fullfile(nofig, 'gain_index.png'));
    render_index(isb, 'iso',         '',         fullfile(nofig, 'iso_index.png'));

    fprintf('\n=== 逐節點 gain/iso（R≤150µm 球 %d 節點）===\n', npts);
    fprintf('Var  gain: fix=%.5g  bias=%.5g (mT/A)^2 | iso: fix=%.5g  bias=%.5g\n', ...
            var(gf), var(gb), var(isf), var(isb));
    fprintf('Mean gain: fix=%.4f  bias=%.4f mT/A     | iso: fix=%.4f  bias=%.4f\n', ...
            mean(gf), mean(gb), mean(isf), mean(isb));
end

function render_index(y, ylab, vunit, out)
    v = var(y);  m = mean(y);  n = numel(y);
    fig = figure('Color','w','Position',[80 80 1000 620]); hold on;
    plot(1:n, y, '.', 'MarkerSize', 4, 'Color', [0.10 0.35 0.75]);
    yline(m, '--', 'Color', [0.85 0.20 0.20], 'LineWidth', 2);
    grid off; box on;
    xlim([1 n]);
    set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.012 .012]);
    xlabel('node index','FontWeight','bold');
    ylabel(ylab,'FontWeight','bold');
    if isempty(vunit), vtxt = sprintf('Var = %.4g', v);
    else,              vtxt = sprintf('Var = %.4g %s', v, vunit); end
    txt = sprintf('%s\nmean = %.4g', vtxt, m);
    text(0.030, 0.94, txt, 'Units','normalized', 'FontSize',15,'FontWeight','bold', ...
         'VerticalAlignment','top','BackgroundColor',[1 1 1],'EdgeColor',[0 0 0],'Margin',5);
    ax = gca; ax.Toolbar.Visible = 'off';
    hold off;
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('saved %s\n', out);
    close(fig);
end
