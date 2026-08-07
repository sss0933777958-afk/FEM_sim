function plot_pole_circuit_side(force)
%PLOT_POLE_CIRCUIT_SIDE  六極 P1 下極磁路側視圖（y=0, xz 平面）：填滿 vs 半切，上下兩 panel。
%   plot_pole_circuit_side()      % 有快取就用
%   plot_pole_circuit_side(true)  % 強制重算
%
%   上 panel = lower_filled（完整錐）、下 panel = graded（半切）。
%   **共用單一色階與 colorbar**（figure-style「同類比較共用色階」）。
%
%   箭頭取樣（**混合**，使用者 2026-08-05 指定）：
%     ① 每格優先取「最近 y=0 的真實 FEM 節點」（不內插）
%     ② **空格改用內插補**（磁路在空氣區網格粗、很多格沒節點 → 圖面破碎）
%   ⚠ **本圖含內插**：空格由 scatteredInterpolant('linear','none') 於格心補值。
%     鐵/空氣**分別建內插器**（B 的切向分量在界面不連續，混用會把界面抹平），
%     依解析幾何判該格屬鐵或空氣後查對應內插器。內插點在報告與圖說須標明。
%
%   箭頭：方向 = 單位 B_xz；長度 ∝ (|B|/CAP)^0.25；顏色 = turbo 依 |B| (mT) 28 bin。
%   號誌：下極 raw = sink → 套 s_source = -1 轉 all-source（磁通朝尖端流、由尖端射出）。
%   已數值驗證：錐體內 100% 的箭頭 B_x < 0（朝尖端）；尖端外空氣徑向外分量平均 cos ≈ +0.41。
%
%   幾何（節點雲實測）：錐軸水平 y=0、z = ZAX = -13.000；apex x = 0.408；
%     錐 R(x) = 0.04 + (x-0.408)*0.2028 至 x = 15.235（R = 3.047），
%     其後為**等徑圓柱 R = 3.047**，鋼件約至 x ≈ 17.9 才結束（故輪廓不在 15.235 封口）。
%     半切：整根保留 z <= ZAX（錐段與圓柱段皆是）。
%
%   風格①粗體框圖（框線 3.5）+ |B| colorbar 標準樣式；tick 奇數等距。
%   輸出 figures/paper_fig/Section3_A/pole_circuit_side.png

    if nargin < 1 || isempty(force), force = false; end
    HERE   = fileparts(fileparts(mfilename('fullpath')));
    ROOT   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
    FIGDIR = fullfile(fileparts(HERE), 'paper_fig', 'Section3_A');
    CACHE  = fullfile(HERE, 'data', 'pole_circuit_side.mat');
    if ~exist(FIGDIR,'dir'), mkdir(FIGDIR); end

    ZAX = -13.000;  XT = 0.408;  XC = 15.235;  RT = 0.04;  SLP = 0.2028;
    RMAX = RT + (XC-XT)*SLP;                     % 3.047 mm（圓柱段半徑）
    XSTEEL = 17.9;                               % 鋼件末端（節點雲實測 ~17.8-18.0）
    SGN = -1;
    XL  = [-1.0 17.0];   ZL = [-17.0 -9.0];      % 視野 [mm]
    CELL = 0.26;  SLAB = 0.5;  SLAB_I = 1.5;     % 取樣薄板 / 內插源薄板
    VAR = {'lower_filled','graded'};  LBL = {'Filled cone','Half-cut cone'};
    geo = struct('ZAX',ZAX,'XT',XT,'XC',XC,'RT',RT,'SLP',SLP,'RMAX',RMAX,'XSTEEL',XSTEEL);

    if ~force && exist(CACHE,'file')
        S = load(CACHE);
        if isequal(S.XL,XL) && isequal(S.ZL,ZL) && S.CELL==CELL && isfield(S,'ISREAL')
            fprintf('cache hit: %s\n', CACHE);
        else, S = []; end
    else, S = []; end

    if isempty(S)
        D = struct('X',{},'Z',{},'BX',{},'BZ',{},'BM',{},'real',{});
        for k = 1:2
            [x,y,z,B] = load_node_B(fullfile(ROOT,VAR{k},'coil1'));
            B = SGN*B;                                                   % → all-source
            isH = strcmp(VAR{k},'graded');
            nx = round(diff(XL)/CELL);  nz = round(diff(ZL)/CELL);
            xe = linspace(XL(1),XL(2),nx+1);  ze = linspace(ZL(1),ZL(2),nz+1);
            xcen = (xe(1:end-1)+xe(2:end))/2;   zcen = (ze(1:end-1)+ze(2:end))/2;

            % ---- ① 真實節點（每格取最近 y=0）----
            in = x>=XL(1)&x<=XL(2)&z>=ZL(1)&z<=ZL(2)&abs(y)<SLAB;
            gi = find(in);
            ix = discretize(x(gi),xe);  iz = discretize(z(gi),ze);
            ok = ~isnan(ix)&~isnan(iz);  gi=gi(ok); ix=ix(ok); iz=iz(ok);
            cid = (ix-1)*nz + iz;   ay = abs(y(gi));
            g = findgroups(cid);
            pick = accumarray(g,(1:numel(cid))',[],@(r) r(find(ay(r)==min(ay(r)),1)));
            sel = gi(pick);   cell_real = unique(cid);
            Xr=x(sel); Zr=z(sel); Br=B(sel,:);

            % ---- ② 空格用內插補（鐵/空氣分開）----
            allc = (1:nx*nz)';   empt = setdiff(allc, cell_real);
            % cid = (ix-1)*nz + iz  →  ix = ceil(cid/nz)、iz = mod(cid-1,nz)+1
            xq = xcen(ceil(empt/nz)).';  zq = zcen(mod(empt-1,nz)+1).';
            % [MODIFIED 2026-08-05] 源範圍由 ±1 放大到 ±4 mm：'none' 不外插，視野角落若落在
            %   源點凸包外會回 NaN → 圖面出現白塊（實測 graded 有 70 格）。放大源範圍讓整個
            %   視野都在凸包內即可補滿，且**仍不外插**（不捏造凸包外的值）。
            src = abs(y) < SLAB_I & x>=XL(1)-4 & x<=XL(2)+4 & z>=ZL(1)-4 & z<=ZL(2)+4;
            ironN = src & is_steel(x,y,z,geo,isH);
            airN  = src & ~is_steel(x,y,z,geo,isH);
            Bq = nan(numel(xq),3);
            ironQ = is_steel(xq, zeros(size(xq)), zq, geo, isH);
            for w = 1:2
                m = (w==1);  qm = ironQ==m;   nd = (w==1)*ironN + (w==2)*airN;  nd = logical(nd);
                if ~any(qm) || sum(nd)<100, continue, end
                F = scatteredInterpolant(x(nd),y(nd),z(nd),B(nd,1),'linear','none');
                Fy = F; Fy.Values=B(nd,2);  Fz = F; Fz.Values=B(nd,3);
                Bq(qm,1)=F(xq(qm),zeros(sum(qm),1),zq(qm));
                Bq(qm,2)=Fy(xq(qm),zeros(sum(qm),1),zq(qm));
                Bq(qm,3)=Fz(xq(qm),zeros(sum(qm),1),zq(qm));
            end
            good = ~any(isnan(Bq),2);
            D(k) = struct('X',[Xr; xq(good)], 'Z',[Zr; zq(good)], ...
                          'BX',[Br(:,1); Bq(good,1)], 'BZ',[Br(:,3); Bq(good,3)], ...
                          'BM',[vecnorm(Br,2,2); vecnorm(Bq(good,:),2,2)]*1e3, ...
                          'real',[true(numel(Xr),1); false(sum(good),1)]);
            fprintf('%-13s 格 %d | 真實 %d + 內插 %d (%.0f%% 內插) | |B| p92 %.1f mT\n', ...
                    VAR{k}, nx*nz, numel(Xr), sum(good), 100*sum(good)/(numel(Xr)+sum(good)), ...
                    prctile(D(k).BM,92));
        end
        S = struct('D',D,'XL',XL,'ZL',ZL,'CELL',CELL,'VAR',{VAR},'LBL',{LBL},'ISREAL',true);
        save(CACHE,'-struct','S');  fprintf('saved %s\n',CACHE);
    end

    %% ---- 共用色階 ----
    CAP = ceil(max(arrayfun(@(d) prctile(d.BM,92), S.D))/50)*50;
    fprintf('shared CAP = %g mT\n', CAP);

    nb = 28;  edges = linspace(0,CAP,nb+1);  cmap = turbo(nb);
    FS = 36;  LW = 3.5;  W = 1500;  Hp = 430;  y0 = 120;  gapv = 70;  x0 = 160;
    fig = figure('Position',[40 40 W 2*Hp+gapv+y0+80],'Color','w');

    for k = 1:2
        ax = axes(fig,'Units','pixels');
        ax.Position = [x0, y0 + (2-k)*(Hp+gapv), W-x0-200, Hp];
        hold(ax,'on');
        d = S.D(k);
        bxz = hypot(d.BX,d.BZ);  bxz(bxz<eps)=eps;
        len = 0.90*CELL*(min(d.BM,CAP)/CAP).^0.25;
        u = d.BX./bxz.*len;  v = d.BZ./bxz.*len;
        ib = discretize(min(d.BM,CAP), edges);  ib(isnan(ib))=nb;
        for b = 1:nb
            m = ib==b;  if ~any(m), continue, end
            quiver(ax, d.X(m)-u(m)/2, d.Z(m)-v(m)/2, u(m), v(m), 0, ...
                   'Color',cmap(b,:), 'LineWidth',1.2, 'MaxHeadSize',0.6, 'AutoScale','off');
        end
        draw_outline(ax, k==2, geo, S.XL(2), LW);
        hold(ax,'off');
        daspect(ax,[1 1 1]);  xlim(ax,S.XL);  ylim(ax,S.ZL);
        set(ax,'XTick',0:4:16,'YTick',[-16 -13 -10]);       % 皆奇數個等距
        set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LW,'TickLength',[.012 .012]);
        box(ax,'on'); grid(ax,'off');
        ylabel(ax,'$\mathbf{z\;(mm)}$','Interpreter','latex','FontSize',FS);
        if k==2, xlabel(ax,'$\mathbf{x\;(mm)}$','Interpreter','latex','FontSize',FS); end
        text(ax, S.XL(1)+0.4, S.ZL(2)-1.0, S.LBL{k}, 'FontSize',30,'FontWeight','bold');
        ax.Toolbar.Visible='off';
        if k==1
            colormap(ax,turbo);  clim(ax,[0 CAP]);
            cb = colorbar(ax,'Units','pixels');
            cb.Position = [W-185, y0, 26, 2*Hp+gapv];
            style_cbar(cb, FS);        % 專案 canonical helper（只設字級/標題，**不動線寬**）
            ax.Position = [x0, y0 + (2-k)*(Hp+gapv), W-x0-200, Hp];
        end
    end

    out = fullfile(FIGDIR,'pole_circuit_side.png');
    exportgraphics(fig,out,'Resolution',150);  fprintf('wrote %s\n',out);  close(fig);
end

% ================= local =================

function style_cbar(cb, FS)
% colorbar 樣式（與 plot_circuit_side.m 逐字相同的 canonical helper）：
%   標題用標準數學字體 \mathbf（latex）、刻度數字 Helvetica 粗體；**不設 LineWidth**。
    cb.FontSize = FS;  cb.FontWeight = 'bold';
    cb.Label.Interpreter = 'latex';  cb.Label.String = '$\mathbf{|B|\;(mT)}$';  cb.Label.FontSize = FS;
end

function m = is_steel(x, y, z, g, isHalf)
% P1 下極鋼件（錐段 + 等徑圓柱段）；y 可為 0（切面查詢）
    zl = z - g.ZAX;   r = hypot(y, zl);
    Rp = g.RT + (x - g.XT)*g.SLP;   Rp = min(Rp, g.RMAX);
    m = x >= g.XT & x <= g.XSTEEL & r <= Rp;
    if isHalf, m = m & zl <= 0; end
end

function draw_outline(ax, isHalf, g, xright, LW)
% 輪廓畫到視野右緣（鋼件延續至 x≈17.9，故**不**在 15.235 封口）
    xr = min(xright, g.XSTEEL);
    xx = linspace(g.XT, min(g.XC,xr), 300);   R = g.RT + (xx-g.XT)*g.SLP;
    plot(ax, xx, g.ZAX-R, 'k-', 'LineWidth', LW);
    if xr > g.XC, plot(ax, [g.XC xr], g.ZAX-g.RMAX*[1 1], 'k-', 'LineWidth', LW); end
    if isHalf
        plot(ax, [g.XT xr], [g.ZAX g.ZAX], 'k-', 'LineWidth', LW);       % 削平面
    else
        plot(ax, xx, g.ZAX+R, 'k-', 'LineWidth', LW);
        if xr > g.XC, plot(ax, [g.XC xr], g.ZAX+g.RMAX*[1 1], 'k-', 'LineWidth', LW); end
    end
end

function [x,y,z,B] = load_node_B(dir_coil)
    C = readmatrix(fullfile(dir_coil,'coil1_coord_all.dat'),'FileType','text');
    C = C(~any(isnan(C(:,1:4)),2),:);
    t = fileread(fullfile(dir_coil,'coil1_bfield_all.dat'));
    t = regexprep(t,'(\d)([-+])','$1 $2');
    tf=[tempname '.dat']; f=fopen(tf,'w'); fprintf(f,'%s',t); fclose(f);
    M = readmatrix(tf,'FileType','text'); delete(tf);
    M = M(~any(isnan(M(:,1:5)),2),:);
    [~,ic,ib] = intersect(C(:,1),M(:,1));
    x=C(ic,2)*1e3; y=C(ic,3)*1e3; z=C(ic,4)*1e3;  B=M(ib,2:4);
end
