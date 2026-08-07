function plot_flux_vs_s(WHICH)
% plot_flux_vs_s -- 磁極截面通量沿極軸的分布 Phi(s)（Maxwell）
% =========================================================================
%   在**垂直極軸的平面**上逐站位積分 Phi(s) = ∫ B·â dA（â 朝極尖為正）：
%     P1 下極（半切）-> 半圓截面 A = pi·R²/2
%     P2 上極（完整）-> 整圓截面 A = pi·R²
%   R(s) = s·tan(beta)，beta = atan(3/15) = 11.31°（**無 R_TIP 偏移**，已用原始格點驗證）。
%   站位錨在 s = 4.4832 mm（= 沿錐面斜距 4.572 mm，sensor 站位），步長 0.1 mm。
%
%   每站位：0.02 x 0.02 mm 直角格、格心取值、邊界格用 10x10 子點求「鋼內面積比 + 形心」
%   加權（面積誤差中位 -0.001%）。B 由**鋼件側**格點內插（r <= R(s)、|B| > 40 mT；
%   P1 另加 upc <= 0 只留鋼側半邊）。
%
%   ⚠ 這是**內插**（Maxwell 匯出 0.1 mm 規則格）。
%   ⚠ s < 2.5 不可用：截面內格點數 ≈ 12.6·s²，s=1 時只有 13 點，撐不起 0.02 mm 積分格。
%   ⚠ 單位 Phi[µWb] = B[T] x A[mm²]。
%
%   風格①粗體框圖；橫軸依 figure-style 2026-08-06：首末點貼框、端點只標數字不畫 tick。
% =========================================================================
    clc;
    here   = fileparts(fileparts(mfilename('fullpath')));
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    if nargin < 1 || isempty(WHICH), WHICH = 'inner'; end
    C = load(fullfile(here,'data','flux_vs_s_maxwell.mat'));
    s = C.SS2(:);
    if strcmpi(WHICH,'surface')     % 每站位一個 BW 寬環帶的表面通量；逐點相加 = ΔΦ_inner
        P1 = C.PHS(:,1);  P2 = C.PHS(:,2);
        ylab = '$\mathbf{\Phi_{surface}\;(\mu Wb)}$';  fn = 'flux_surface_vs_s_maxwell';
        fprintf('Φ_surface：每站位一個 %.2f mm 寬環帶（δ = %.2f mm）的側表面通量\n', C.BW, C.DLR);
        for pk = 1:2
            fprintf('  P%d: Σ Φ_surface = %.4f | ΔΦ_inner = %.4f uWb | 差 %+.2f%%\n', pk, ...
                sum(P1*(pk==1)+P2*(pk==2)) - (P1(1)*(pk==1)+P2(1)*(pk==2)), ...
                C.PHI3(end,pk)-C.PHI3(1,pk), ...
                ((sum(C.PHS(2:end,pk)))/(C.PHI3(end,pk)-C.PHI3(1,pk))-1)*100);
        end
    else
        P1 = C.PHI3(:,1);  P2 = C.PHI3(:,2);
        ylab = '$\mathbf{\Phi_{inner}\;(\mu Wb)}$';  fn = 'flux_inner_vs_s_maxwell';
    end

    i0 = find(abs(s-4.4832) < 1e-6, 1);
    fprintf('s = %.2f ~ %.2f mm，%d 站位，步長 %.2f mm\n', s(1), s(end), numel(s), s(2)-s(1));
    fprintf('  錨點 s = %.4f : P1 %.5f | P2 %.5f   (P1/P2 = %.4f)\n', s(i0), P1(i0), P2(i0), P1(i0)/P2(i0));
    fprintf('  端點 s = %.2f : P1 %.5f | P2 %.5f   |  s = %.2f : P1 %.5f | P2 %.5f\n', ...
            s(1), P1(1), P2(1), s(end), P1(end), P2(end));

    FS  = 36;
    BLU = [0.10 0.35 1.00];   RED = [0.85 0.10 0.10];

    if strcmpi(WHICH,'both')                 % 雙軸：左 Φ_inner（圓）+ 右 <b·â>（方）
        A1 = pi*(C.GEO(1).R0 + C.GEO(1).sl*s).^2/2;    % 下極半圓截面 [mm²]
        A2 = pi*(C.GEO(2).R0 + C.GEO(2).sl*s).^2;      % 上極整圓
        D1 = C.PHI3(:,1)./A1*1000;  D2 = C.PHI3(:,2)./A2*1000;   % <b·â> [mT]
        fprintf('<b·â>: P1 %.1f -> %.1f mT | P2 %.1f -> %.1f mT\n', D1(1),D1(end),D2(1),D2(end));
        fprintf('  比值 P1/P2: %.3f (s=%.1f) ~ %.3f (s=%.1f)\n', D1(1)/D2(1),s(1),D1(end)/D2(end),s(end));
        mi = 1:10:numel(s);
        fig = figure('Color','w','Position',[100 100 1250 900]);  ax = axes(fig);
        yyaxis(ax,'left');   hold(ax,'on');
        h1 = plot(ax,s,C.PHI3(:,1),'-o','Color',BLU,'LineWidth',3.0,'MarkerSize',11, ...
                  'MarkerFaceColor',BLU,'MarkerEdgeColor','k','MarkerIndices',mi,'Clipping','off');
        h2 = plot(ax,s,C.PHI3(:,2),'-o','Color',RED,'LineWidth',3.0,'MarkerSize',11, ...
                  'MarkerFaceColor',RED,'MarkerEdgeColor','k','MarkerIndices',mi,'Clipping','off');
        [yrL,ytL]=ylim_pick([C.PHI3(:,1);C.PHI3(:,2)]);  ylim(ax,yrL); set(ax,'YTick',ytL);
        ylabel(ax,'$\mathbf{\Phi_{inner}\;(\mu Wb)}$','Interpreter','latex','FontSize',FS);
        yyaxis(ax,'right');  hold(ax,'on');
        h3 = plot(ax,s,D1,'-s','Color',BLU,'LineWidth',3.0,'MarkerSize',12, ...
                  'MarkerFaceColor',BLU,'MarkerEdgeColor','k','MarkerIndices',mi,'Clipping','off');
        h4 = plot(ax,s,D2,'-s','Color',RED,'LineWidth',3.0,'MarkerSize',12, ...
                  'MarkerFaceColor',RED,'MarkerEdgeColor','k','MarkerIndices',mi,'Clipping','off');
        [yrR,ytR]=ylim_pick([D1;D2]);  ylim(ax,yrR); set(ax,'YTick',ytR);
        ylabel(ax,'$\mathbf{\langle b\cdot\hat{a}\rangle\;(mT)}$','Interpreter','latex','FontSize',FS);
        ax.YAxis(1).Color='k'; ax.YAxis(2).Color='k';
        box(ax,'on'); grid(ax,'off');
        set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');
        xe=[s(1) s(end)]; xlim(ax,xe); set(ax,'XTick',[4 8 12]);
        yoff = yrR(1) - 0.022*diff(yrR);
        for kk=1:2, text(ax,xe(kk),yoff,sprintf('%.2f',xe(kk)),'HorizontalAlignment','center', ...
              'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off'); end
        lg = legend([h1 h2 h3 h4], {'P1 flux','P2 flux','P1 density','P2 density'}, ...
                    'Interpreter','tex','Location','northwest','NumColumns',2);
        lg.FontSize=24; lg.FontWeight='bold'; lg.Box='on'; lg.EdgeColor='k'; lg.LineWidth=2.5;
        xlabel(ax,'$\mathbf{Position\;along\;pole\;axis\;(mm)}$','Interpreter','latex','FontSize',FS);
        ax.Toolbar.Visible='off';
        out = fullfile(figdir,'flux_inner_density_vs_s_maxwell.png');
        exportgraphics(fig,out,'Resolution',150);  fprintf('wrote %s\n',out);  return;
    end

    fig = figure('Color','w','Position',[100 100 1250 900]);  ax = axes(fig);  hold(ax,'on');
    h1 = plot(ax, s, P1, '-', 'Color',BLU, 'LineWidth',3.5, 'Clipping','off');
    h2 = plot(ax, s, P2, '-', 'Color',RED, 'LineWidth',3.5, 'Clipping','off');

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');

    xe = [s(1) s(end)];                       % 首末點貼框
    xlim(ax, xe);  set(ax,'XTick',[4 8 12]);  % 3 個(奇數)等距、端點不進 tick
    [yr, yt] = ylim_pick([P1;P2]);
    ylim(ax, yr);  set(ax,'YTick', yt);

    yoff = yr(1) - 0.022*diff(yr);            % 端點只標數字、不畫 tick
    for k = 1:2
        text(ax, xe(k), yoff, sprintf('%.2f',xe(k)), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end

    lg = legend([h1 h2], {'P1 (half-cut)','P2 (full cone)'}, 'Interpreter','tex', 'Location','northwest');
    lg.FontSize = 26;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;

    xlabel(ax, '$\mathbf{Position\;along\;pole\;axis\;(mm)}$', 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, ylab, 'Interpreter','latex', 'FontSize',FS);
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    out = fullfile(figdir, [fn '.png']);
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function [yr, yt] = ylim_pick(v)
% 奇數個等距 nice tick + 兩端留白 = 間距；資料全為正時不讓下界落到負值
    lo = min(v);  hi = max(v);   nonneg = lo >= 0;
    for n = [3 5]
        s = nice_step((hi-lo)/(n+1));
        % 需涵蓋的是「框」不是「tick 帶」：框寬 = (n+1)*s（tick 帶 + 兩端各留白 s）
        while (n+1)*s < (hi-lo), s = next_step(s); end
        for shift = 0:2                                  % 微調中心，找出下界不為負的擺法
            ctr = (round((lo+hi)/2/s) + shift)*s;
            yt  = ctr + (-(n-1)/2:(n-1)/2)*s;
            yr  = [yt(1)-s, yt(end)+s];
            ok  = yr(1) <= lo && yr(2) >= hi && (hi-lo)/diff(yr) > 0.35;
            if ok && (~nonneg || yr(1) >= -1e-12), return; end
        end
    end
end

function s = nice_step(x)
    k = floor(log10(max(x, realmin)));  m = x/10^k;
    cand = [1 2 2.5 5 10];  [~,i] = min(abs(cand - m));  s = cand(i)*10^k;
end

function s2 = next_step(s)
    k = floor(log10(s));  m = round(s/10^k, 4);
    cand = [1 2 2.5 5 10];  i = find(cand > m + 1e-9, 1);
    if isempty(i), s2 = 10^(k+1); else, s2 = cand(i)*10^k; end
end
