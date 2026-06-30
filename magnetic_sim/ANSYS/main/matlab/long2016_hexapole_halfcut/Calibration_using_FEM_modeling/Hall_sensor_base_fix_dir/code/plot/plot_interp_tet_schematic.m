function plot_interp_tet_schematic(PREVIEW)
% PLOT_INTERP_TET_SCHEMATIC  重心法內插「示意圖」：一顆四面體 + 內部一個取樣點。
% -------------------------------------------------------------------------
%   說明 standard 粗網格內插的原理（非真實資料、純示意）：
%   取樣點 p 落在一顆四面體內，用 4 個角節點的場 B1..B4 依重心權重 λ 內插：
%       B(p) = λ1·B1 + λ2·B2 + λ3·B3 + λ4·B4 ,  Σλ_i = 1
%   λ_i = (把第 i 角換成 p 的子四面體體積) / (整顆 tet 體積)。
%   PREVIEW=true → 暫存 png；false → 存 figures/。

    if nargin < 1 || isempty(PREVIEW), PREVIEW = true; end
    DPI = 200;

    % ---- 一顆漂亮的四面體 ----
    V = [ 0.00 0.00 0.00 ;     % B1
          1.00 0.00 0.00 ;     % B2
          0.50 0.866 0.00 ;    % B3
          0.50 0.289 0.816 ];  % B4
    lam = [0.34 0.24 0.20 0.22]; lam = lam/sum(lam);   % 取樣點的重心權重
    p   = lam * V;                                     % p = Σ λ_i V_i

    col = [0.20 0.45 0.85; 0.90 0.55 0.10; 0.20 0.65 0.30; 0.65 0.30 0.70]; % 4 角配色

    fig = figure('Position',[80 80 980 860],'Color','w');
    ax = axes(fig,'Position',[0.06 0.08 0.88 0.84]); hold(ax,'on');

    % 四面體面（半透明）+ 邊
    F = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
    patch('Vertices',V,'Faces',F,'FaceColor',[0.6 0.7 0.9],'FaceAlpha',0.12, ...
          'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.6);

    % p → 4 角的虛線 + λ 標註
    nm = {'_1','_2','_3','_4'};
    for i = 1:4
        plot3([p(1) V(i,1)],[p(2) V(i,2)],[p(3) V(i,3)],'--', ...
              'Color',[col(i,:) ],'LineWidth',1.6);
        mid = 0.55*p + 0.45*V(i,:);
        text(mid(1),mid(2),mid(3), sprintf('\\lambda%s=%.2f',nm{i},lam(i)), ...
             'Color',col(i,:),'FontSize',13,'FontWeight','bold');
    end

    % 4 角節點
    for i = 1:4
        plot3(V(i,1),V(i,2),V(i,3),'o','MarkerSize',13, ...
              'MarkerFaceColor',col(i,:),'MarkerEdgeColor','k','LineWidth',1.2);
        off = (V(i,:)-p); off = 0.13*off/norm(off);
        text(V(i,1)+off(1),V(i,2)+off(2),V(i,3)+off(3), sprintf('B%s',nm{i}), ...
             'Color',col(i,:),'FontSize',16,'FontWeight','bold','Interpreter','tex');
    end

    % 取樣點 p
    plot3(p(1),p(2),p(3),'p','MarkerSize',20,'MarkerFaceColor',[0.85 0.1 0.1], ...
          'MarkerEdgeColor','k','LineWidth',1.2);
    text(p(1)+0.03,p(2)-0.06,p(3)+0.05,' p','Color',[0.85 0.1 0.1], ...
         'FontSize',16,'FontWeight','bold');

    % 公式
    text(0.02, 0.06, ...
        'B(p) = \lambda_1B_1 + \lambda_2B_2 + \lambda_3B_3 + \lambda_4B_4 ,   \Sigma\lambda_i = 1', ...
        'Units','normalized','FontSize',15,'FontWeight','bold','Interpreter','tex');

    axis(ax,'equal'); grid(ax,'on'); box(ax,'on'); view(-37,16);
    set(ax,'XTick',[],'YTick',[],'ZTick',[]);
    title('重心法內插示意：一顆四面體內一個取樣點 p（B(p) 由 4 角節點依重心權重內插）', ...
          'Interpreter','none','FontSize',13);
    ax.Toolbar.Visible = 'off';

    if PREVIEW
        out = fullfile(tempdir,'interp_tet_schematic_preview.png');
    else
        out = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\figures\interp_tet_schematic.png'];
    end
    exportgraphics(fig, out, 'Resolution', DPI);
    fprintf('saved: %s\n', out);
end
