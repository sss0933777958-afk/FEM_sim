function plot_tol_overlay(matname, xlab, unitfmt)
%% plot_tol_overlay -- 通用 region_error 疊圖（P1 vs all-6）+ <2% 縱向虛線容忍框
% =========================================================================
% 讀 data/<matname>.mat（欄位 x / eps_P1 / eps_all / errpct）→ 疊 P1(藍) + all-6(紅) 曲線，
% 以 all-6（最壞情況）region_error<2% 的窗畫縱向虛線框（避開頂部圖例）。風格 ①粗體框圖。
% 輸出 → figures/eighteen_param/<matname>_overlay.png（覆蓋迭代）。
%   matname : .mat 檔名（不含 .mat）；xlab: x 軸標籤；unitfmt: 框標籤 sprintf 格式（含兩個 %，如 '%.2f-%.2f mm'）。
% =========================================================================
    here = fileparts(mfilename('fullpath'));  TREE = fileparts(fileparts(here));
    S = load(fullfile(TREE,'data',[matname '.mat']));
    x = S.x(:);  eP1 = S.eps_P1(:);  eAll = S.eps_all(:);
    fig_dir = fullfile(TREE,'figures','eighteen_param'); if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

    col1=[0 .20 .60]; col6=[.85 .10 .10];
    fig=figure('Color','w','Units','inches','Position',[1 1 9 5.4]); ax=axes(fig); hold(ax,'on');
    h6=plot(ax,x,eAll,'-o','Color',col6,'LineWidth',2.5,'MarkerFaceColor',col6,'MarkerSize',5);
    h1=plot(ax,x,eP1, '-o','Color',col1,'LineWidth',2.5,'MarkerFaceColor',col1,'MarkerSize',5);
    emax=max([eP1;eAll]); cands=[1 2 5 10 20 25 50];
    ystep=cands(find(ceil(emax./cands)<=12,1)); ytop=ceil(emax/ystep)*ystep;
    [~,imn]=min(eAll); xmn=x(imn); xc=local_cross(x,eAll,2);              % <2% 容忍窗（all-6）
    xLc=xc(xc<xmn); if isempty(xLc), xL=min(x); else, xL=max(xLc); end
    xRc=xc(xc>xmn); if isempty(xRc), xR=max(x); else, xR=min(xRc); end
    rectangle(ax,'Position',[xL 0 xR-xL 0.65*ytop],'EdgeColor',[.12 .12 .12],'LineStyle','--','LineWidth',2.2);
    text(ax,xmn,0.72*ytop,sprintf(['< 2%%\n' unitfmt],xL,xR),'HorizontalAlignment','center', ...
         'FontSize',12.5,'FontWeight','bold','Color',[.12 .12 .12],'Interpreter','none');
    hold(ax,'off');
    xlim(ax,[min(x) max(x)]); ylim(ax,[0 ytop]);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]); box(ax,'on'); grid(ax,'off');
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end)); set(ax,'YTick',0:ystep:ytop);
    xlabel(ax,xlab,'FontWeight','bold'); ylabel(ax,'relative RMS error  (%)','FontWeight','bold');
    lg=legend([h1 h6],{'P1','all 6'},'Location','north'); lg.FontSize=14; lg.FontWeight='bold';
    out=fullfile(fig_dir,[matname '_overlay.png']); exportgraphics(fig,out,'Resolution',150);
    fprintf('已輸出 %s（<2%% 窗 = [%.3f, %.3f]）\n', out, xL, xR);
end

function xc = local_cross(x,y,L)
    xc=[];
    for i=1:numel(x)-1
        if (y(i)-L)*(y(i+1)-L)<0, xc(end+1,1)=x(i)+(x(i+1)-x(i))*(L-y(i))/(y(i+1)-y(i)); end %#ok<AGROW>
    end
end
