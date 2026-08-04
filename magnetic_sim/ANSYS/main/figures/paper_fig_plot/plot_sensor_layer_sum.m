function plot_sensor_layer_sum(NPL)
% plot_sensor_layer_sum -- 感測圓柱逐層 Σ(b·n̂⁺)：P1 上平面 vs P2 上錐面
% =========================================================================
%   圓柱（R=0.15mm、H=0.10mm）沿 n̂⁺ 從**底面**往上切 10 層：層 k = [(k-1)*0.01, k*0.01] mm。
%   **每層在該層體積內均勻撒 NPL 點**（預設 50）：a~U(層下緣, 層上緣)、r=R√U、θ~U(0,2π)。
%   層值 = 該層 NPL 點的 Σ(b·n̂⁺) [mT]；兩個區域共用同一組局部座標 → 逐層配對可比。
%   ⚠ 橫軸把層 k 標在**層上緣** k*0.01（故軸為 0.01…0.10）；該層點的平均高度是層中心，
%     比標示值低 0.005mm —— 由圖讀交叉點時要扣掉這 0.005mm（腳本末會直接印正確值）。
%
%   ⚠ **本圖為內插**：Maxwell sensor 區匯出格距 0.1mm，柱內真格點僅 ~7-14 個，
%     故以 scatteredInterpolant（鄰域 1.5mm，同 build_V_matrix 的 'scattered'）內插。
%
%   風格：figure-style 選項①粗體框（字級 36、box、無 grid、單位 ()、tick 奇數等距、
%         兩端留白 = 間距）。輸出 → figures/paper_fig/Section3_A/sensor_layer_sum_maxwell.png
% =========================================================================
    if nargin < 1 || isempty(NPL), NPL = 50; end         % 每層撒點數
    clc;
    here   = fileparts(mfilename('fullpath'));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    MW = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(MW,'function'));

    cfg  = model_config('long2016_hexapole_halfcut','tip40um');
    SOFF = 4.572e-3;  RS = 0.15e-3;  HS = 0.10e-3;  DL = 0.01e-3;   % 層厚 [m]
    RLOC = 1.5e-3;    NL = round(HS/DL);                            % 層數 = 10
    FS   = 36;

    cfn = fullfile(here, sprintf('sensor_layer_sum_npl%d_unif.mat', NPL));
    if exist(cfn,'file')
        L = load(cfn,'S');  S = L.S;  fprintf('loaded cache %s\n', cfn);
    else
        % ---- 每層在該層「體積」內均勻撒 NPL 點（層 k = [(k-1)*DL, k*DL]，從底面往上）----
        rng(0);
        aL = zeros(NPL*NL,1);  rrL = aL;  thL = aL;  lay = aL;
        for k = 1:NL
            ix = (k-1)*NPL + (1:NPL);
            aL(ix)  = (k-1)*DL + DL*rand(NPL,1);             % 層內高度均勻
            rrL(ix) = RS*sqrt(rand(NPL,1));                  % √U → 面積均勻
            thL(ix) = 2*pi*rand(NPL,1);
            lay(ix) = k;
        end
        fprintf('每層 %d 點 × %d 層 = %d 點（層內體積均勻撒）\n', NPL, NL, NPL*NL);

        S = zeros(NL,2);
        for c = 1:2                                          % 1 = P1 平切上表面, 2 = P2 上錐面
            if c == 1, pidx = 1; FACE = 'flat'; else, pidx = 2; FACE = 'cone'; end
            [sc, nh] = sensor_geom(cfg, pidx, SOFF, FACE);
            ci = sc + [0;0;cfg.SPH_OFST];                    % WP → Maxwell 框
            d  = import_maxwell_fld(fullfile(cfg.fld_dir, cfg.fld_files_voltage{pidx}));
            m  = abs(d.x-ci(1))<=RLOC & abs(d.y-ci(2))<=RLOC & abs(d.z-ci(3))<=RLOC;
            P  = [d.x(m), d.y(m), d.z(m)];
            Fx = scatteredInterpolant(P, d.bx(m), 'linear','none');
            Fy = scatteredInterpolant(P, d.by(m), 'linear','none');
            Fz = scatteredInterpolant(P, d.bz(m), 'linear','none');
            clear d P;
            t1 = [-nh(2); nh(1); 0];  if norm(t1) < 1e-9, t1 = [1;0;0]; end
            t1 = t1/norm(t1);  t2 = cross(nh, t1);
            Q  = ci.' + aL.*nh.' + (rrL.*cos(thL)).*t1.' + (rrL.*sin(thL)).*t2.';
            B  = cfg.s_source(pidx) * [Fx(Q), Fy(Q), Fz(Q)] * 1e3;        % mT
            bn = B*nh;   bn(~isfinite(bn)) = 0;
            S(:,c) = accumarray(lay, bn, [NL 1]);
            fprintf('  區域%d: 層1 sum=%.1f  層%d sum=%.1f mT\n', c, S(1,c), NL, S(NL,c));
        end
        save(cfn,'S');  fprintf('saved cache %s\n', cfn);
    end

    %% ---- 繪圖（選項①粗體框）----
    kk = (1:NL).' * DL * 1e3;                  % 各層標在其上緣 [mm]：0.01 … 0.10
    cA = [0.10 0.35 1.00];   cB = [0.85 0.10 0.10];
    fig = figure('Color','w','Position',[80 60 1180 900]);  ax = axes(fig);  hold(ax,'on');
    h1 = plot(ax, kk, S(:,1), '-o', 'Color',cA, 'LineWidth',3.0, 'MarkerSize',10, 'MarkerFaceColor',cA);
    h2 = plot(ax, kk, S(:,2), '-s', 'Color',cB, 'LineWidth',3.0, 'MarkerSize',10, 'MarkerFaceColor',cB);

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.018 .018]);
    % 水平軸 = 沿 n̂⁺ 高度，起點 0.01、終點 0.10 mm；3 個奇數等距內縮 tick（s=0.02）
    xlim(ax,[0.01 0.1]);  set(ax,'XTick',[0.04 0.06 0.08]);
    [yl, yt] = axlim_auto(min(S(:)), max(S(:)));
    ylim(ax, yl);  set(ax,'YTick', yt);
    xlabel(ax, '$\mathbf{Height\;along\;\hat{n}^{+}\;(mm)}$', 'Interpreter','latex','FontSize',FS);
    ylabel(ax, '$\mathbf{\sum b \cdot \hat{n}^{+}\;(mT)}$', 'Interpreter','latex','FontSize',FS);
    % 起點 0.01 與終點 0.1 只標數字、不畫 tick（figure-style 慣例 #5）
    for xe = [0.01 0.1]
        text(ax, xe, yl(1)-0.022*diff(yl), sprintf('%g',xe), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
    end
    lg = legend([h1 h2], {'P1 top flat surface','P2 cone surface'}, 'Location','northeast');
    lg.FontSize = 28;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.0;
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    out = fullfile(figdir, 'sensor_layer_sum_maxwell.png');
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
    fprintf('\nP1 上平面: 層1 %.1f → 層10 %.1f mT (%.2f%%)\n', S(1,1), S(NL,1), 100*(S(NL,1)/S(1,1)-1));
    fprintf('P2 上錐面: 層1 %.1f → 層10 %.1f mT (%.2f%%)\n', S(1,2), S(NL,2), 100*(S(NL,2)/S(1,2)-1));
    % 交叉點：用各層「平均高度」= 層中心（比圖上標的層上緣低 DL/2）
    ac = ((1:NL).'-0.5)*DL*1e3;   q1 = polyfit(ac,S(:,1),1);   q2 = polyfit(ac,S(:,2),1);
    xc = (q2(2)-q1(2))/(q1(1)-q2(1));   Hmm = NL*DL*1e3;
    fprintf('\n交叉點（距柱底）= %.4f mm   距磁極表面 = %.4f mm\n', xc, 0.41+xc);
    fprintf('P2 > P1 佔 %.1f%% ；P2 < P1 佔 %.1f%% 的柱高\n', 100*xc/Hmm, 100*(Hmm-xc)/Hmm);
end

% ============================================================================
function [pos, n] = sensor_geom(cfg, i, SOFF, FACE)
% 與 build_V_matrix / plot_sensor_cylinder_B_3d 同式。
    be = atan2(cfg.POLE_R, cfg.POLE_CONE_LEN);
    ps = atan2(cfg.R_norm_z, cfg.R_norm_xy);
    iu = cfg.upper_incline;   AIR = 0.41e-3;
    dr = @(el,az) [cos(el)*cos(az); cos(el)*sin(az); sin(el)];
    th = cfg.pole_angles(i)*pi/180;
    if strcmpi(FACE,'flat')
        e1 = dr(-ps, th);  e2 = dr(0, th);          n = [0;0;1];
    elseif cfg.pole_is_lower(i)
        e1 = dr(-ps, th);  e2 = dr(-be, th);        n = dr(-be-pi/2, th);
    else
        e1 = dr(+ps, th);  e2 = dr(iu+be, th);      n = dr(iu+be+pi/2, th);
    end
    pos = cfg.R_norm*e1 + SOFF*e2 + AIR*n;
end

% ============================================================================
function [lim, tk] = axlim_auto(lo, hi)
% 奇數個等距 tick + **兩端留白 = tick 間距**（figure-style）；取容得下資料的最小跨距。
    cand = [1 2 2.5 3 4 5 10];   best = {};   bestSpan = inf;
    for n = [3 5]
        k0 = floor(log10(max(hi-lo,realmin)/(n+1)));
        for k = (k0-1):(k0+3)
            for c = cand
                s = c*10^k;   ctr = round((lo+hi)/2/s)*s;
                t = ctr + (-(n-1)/2:(n-1)/2)*s;   L = [t(1)-s, t(end)+s];
                if lo >= L(1)+0.12*s && hi <= L(2)-0.12*s && (n+1)*s < bestSpan
                    bestSpan = (n+1)*s;   best = {L, t};
                end
            end
        end
    end
    if isempty(best), s = (hi-lo)/4; t = (lo+hi)/2 + (-1:1)*s; best = {[t(1)-s t(end)+s], t}; end
    lim = best{1};   tk = best{2};
end
