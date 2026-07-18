function plot_sensor_Bn_hist()
% PLOT_SENSOR_BN_HIST  NTU hexapole：sensor 有效範圍內「內插取樣點」的 B·n+ 分布疊圖直方圖。
%   兩張圖：pole_assembly(=singlepole 資料) 與 upper_assembly，各把 S1、S2 疊在同一張。
%   x = B·n+ [mT]；y = (count/總數)×100%（'percentage'）；離散 nb=180 bars、共用 edges、
%   mean 虛線(同色)、legend 顯示 mean；①粗體框圖。★此圖為「內插取樣點」分布（非真實節點）。
%
%   內插機制逐字沿用 field_viz/code/scripts/sensor_average.m：
%     import 'all' → sensor 中心 0.30mm 內真實 FEM 節點建 Delaunay → 圓柱(R=0.15,H=0.10mm 沿 n+)
%     面積均勻撒點 → pointLocation + barycentric 內插 B(p) → Bproj = (B·n+)·1e3 [mT] 每點。
%   S1 center (13.15,-0.99,-0.41) n+=-z；S2 (13.15,-0.99,+0.66) n+=+z。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % import_ansys_data
    RR  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';
    OUT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\NTU_hexapole\field_viz\figures';
    if ~exist(OUT,'dir'); mkdir(OUT); end

    N = 1000;  R = 0.15e-3;  H = 0.10e-3;  OVER = 2*N;     % 內插取樣：每 sensor 1000 點（同 sensor_average）
    cen = [13.15 -0.99 -0.41; 13.15 -0.99 0.66]*1e-3;      % S1,S2 中心 [m]
    nrm = [0 0 -1; 0 0 +1];                                % n+：S1=-z, S2=+z
    colS = [0.72 0.11 0.11; 0.12 0.35 0.70];               % S1 深紅 / S2 深藍

    % model 資料夾 -> 使用者命名(圖檔) 對應
    variants = { 'singlepole',     'pole_assembly',  'sensor_Bn_hist_pole_assembly.png';
                 'upper_assembly', 'upper_assembly', 'sensor_Bn_hist_upper_assembly.png' };

    for v = 1:size(variants,1)
        model = variants{v,1};
        fprintf('\n=== %s (%s) ===\n', model, variants{v,2});
        d = import_ansys_data(fullfile(RR,model,'coil1'),'all',model);       % SI [m],[T]
        X = [d.x d.y d.z];  B3 = [d.bx d.by d.bz];
        Bn = cell(1,2);
        for k = 1:2
            Bn{k} = interp_Bn(X, B3, cen(k,:).', nrm(k,:).', N, R, H, OVER);
            fprintf('  S%d (n+=%s): %d pts  <B.n+>=%.4f (sd %.4f) mT\n', ...
                    k, tern(nrm(k,3)<0,'-z','+z'), numel(Bn{k}), mean(Bn{k}), std(Bn{k}));
        end
        render_overlay(Bn{1}, Bn{2}, colS, N, variants{v,2}, fullfile(OUT, variants{v,3}));

        % --- upper_assembly：另出「原位 + 外移 -y1mm」四條疊圖 ---
        if strcmp(model,'upper_assembly')
            cen_mv = cen;  cen_mv(:,2) = cen(:,2) - 1e-3;      % 兩顆 sensor 皆 -y 1mm
            b1m = interp_Bn(X, B3, cen_mv(1,:).', nrm(1,:).', N, R, H, OVER);
            b2m = interp_Bn(X, B3, cen_mv(2,:).', nrm(2,:).', N, R, H, OVER);
            fprintf('  moved -y1mm: S1 %d pts <B.n+>=%.4f mT | S2 %d pts <B.n+>=%.4f mT\n', ...
                    numel(b1m), mean(b1m), numel(b2m), mean(b2m));
            render_overlay4(Bn{1}, Bn{2}, b1m, b2m, ...
                            fullfile(OUT,'sensor_Bn_hist_upper_assembly_ym1mm.png'));
        end
    end
end

% ---- 四條疊圖：原位 S1/S2(深) + 外移 -y1mm S1/S2(淺) ----
function render_overlay4(s1, s2, s1m, s2m, outpath)
    nb = 180;
    all4 = [s1;s2;s1m;s2m];
    edges = linspace(min(all4), max(all4), nb+1);
    c1=[0.72 0.11 0.11]; c1m=[0.96 0.55 0.30];   % S1 深紅 / 外移 橘紅
    c2=[0.12 0.35 0.70]; c2m=[0.40 0.72 0.98];   % S2 深藍 / 外移 淺藍
    fig = figure('Color','w','Position',[100 100 1000 700]);  ax = axes(fig); hold(ax,'on');
    h1 =histogram(ax,s1, edges,'Normalization','percentage','FaceColor',c1, 'FaceAlpha',0.55,'EdgeColor','w','LineWidth',0.3);
    h1m=histogram(ax,s1m,edges,'Normalization','percentage','FaceColor',c1m,'FaceAlpha',0.55,'EdgeColor','w','LineWidth',0.3);
    h2 =histogram(ax,s2, edges,'Normalization','percentage','FaceColor',c2, 'FaceAlpha',0.55,'EdgeColor','w','LineWidth',0.3);
    h2m=histogram(ax,s2m,edges,'Normalization','percentage','FaceColor',c2m,'FaceAlpha',0.55,'EdgeColor','w','LineWidth',0.3);
    xline(ax,mean(s1), '--','Color',c1, 'LineWidth',2);  xline(ax,mean(s1m),'--','Color',c1m,'LineWidth',2);
    xline(ax,mean(s2), '--','Color',c2, 'LineWidth',2);  xline(ax,mean(s2m),'--','Color',c2m,'LineWidth',2);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    box(ax,'on');  grid(ax,'off');
    yl=ylim(ax);  ylim(ax,[0 yl(2)*1.20]);
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));  yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel(ax,'B\cdot n_+ (mT)','FontWeight','bold','Interpreter','tex');
    ylabel(ax,'Percentage (%)','FontWeight','bold');
    lg=legend(ax,[h1 h1m h2 h2m], { ...
        sprintf('S1  (mean = %.3f mT)',   mean(s1)),  sprintf('S1''  (mean = %.3f mT)', mean(s1m)), ...
        sprintf('S2  (mean = %.3f mT)',   mean(s2)),  sprintf('S2''  (mean = %.3f mT)', mean(s2m))}, ...
        'Location','northoutside','NumColumns',2);
    set(lg,'FontSize',13,'FontWeight','bold','Box','on');
    ax.Toolbar.Visible='off';
    exportgraphics(fig, outpath, 'Resolution',150);
    fprintf('  saved %s\n', outpath);
    close(fig);
end

% ---- 內插取樣點 B·n+ [mT]（沿用 sensor_average.m 機制）----
function Bproj = interp_Bn(X, B3, c, ni, N, R, H, OVER)
    rng(1);
    loc = vecnorm(X - c.', 2, 2) < 0.30e-3;                % 局部真實節點（球內+一點外）
    P = X(loc,:);  BP = B3(loc,:);
    TR = delaunayTriangulation(P);
    Tb = null(ni.');  t1 = Tb(:,1);  t2 = Tb(:,2);         % 圓柱橫切基底 ⊥ n+
    a = H*rand(OVER,1);  r = R*sqrt(rand(OVER,1));  th = 2*pi*rand(OVER,1);   % 面積均勻
    pts = c.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';
    ti = pointLocation(TR, pts);  good = find(~isnan(ti));
    good = good(1:min(N,numel(good)));
    bc = cartesianToBarycentric(TR, ti(good), pts(good,:));
    conn = TR.ConnectivityList(ti(good),:);
    Bint = zeros(numel(good),3);
    for cc = 1:4, Bint = Bint + bc(:,cc).*BP(conn(:,cc),:); end
    Bproj = (Bint*ni)*1e3;                                 % B·n+ [mT]，每點
end

% ---- 疊圖直方圖（render_overlay 樣式：nb=180、percentage、mean 虛線、legend mean、①粗體框）----
function render_overlay(y1, y2, colS, N, tag, outpath)
    nb = 180;
    edges = linspace(min([y1;y2]), max([y1;y2]), nb+1);    % 兩 sensor 共用 edges
    fig = figure('Color','w','Position',[100 100 900 640]);
    ax = axes(fig); hold(ax,'on');
    h1 = histogram(ax, y1, edges, 'Normalization','percentage', ...
                   'FaceColor',colS(1,:), 'FaceAlpha',0.55, 'EdgeColor','w', 'LineWidth',0.3);
    h2 = histogram(ax, y2, edges, 'Normalization','percentage', ...
                   'FaceColor',colS(2,:), 'FaceAlpha',0.55, 'EdgeColor','w', 'LineWidth',0.3);
    m1 = mean(y1);  m2 = mean(y2);
    xline(ax, m1, '--', 'Color',colS(1,:), 'LineWidth',2);
    xline(ax, m2, '--', 'Color',colS(2,:), 'LineWidth',2);
    % ①粗體框圖
    set(ax, 'FontSize',16, 'FontWeight','bold', 'LineWidth',2, 'TickLength',[.018 .018]);
    box(ax,'on');  grid(ax,'off');
    yl = ylim(ax);  ylim(ax,[0 yl(2)*1.20]);               % headroom
    xt = get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
    yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel(ax, 'B\cdot n_+ (mT)', 'FontWeight','bold', 'Interpreter','tex');
    ylabel(ax, 'Percentage (%)', 'FontWeight','bold');     % = (count/total)×100
    lg = legend(ax, [h1 h2], {sprintf('S1  (mean = %.3f mT)', m1), ...
                              sprintf('S2  (mean = %.3f mT)', m2)}, 'Location','northeast');
    set(lg, 'FontSize',15, 'FontWeight','bold', 'Box','on');
    ax.Toolbar.Visible = 'off';
    exportgraphics(fig, outpath, 'Resolution',150);
    fprintf('  saved %s\n', outpath);
    close(fig);
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
