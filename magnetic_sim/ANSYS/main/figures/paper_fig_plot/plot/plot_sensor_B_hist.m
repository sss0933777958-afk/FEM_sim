function plot_sensor_B_hist(force, NPT, SOFF, NB, PAIR)
% plot_sensor_B_hist -- Hall 感測圓柱內 ‖b‖ 分布疊圖(Maxwell)；同時是本組圖的**取樣來源**
% =========================================================================
%   三個貼附面(每面各一顆自激發、同尺寸圓柱 R0.15×H0.10mm、離面氣隙 0.41mm、SOFF 4.572mm)：
%     'p1flat' = P1 半切**平切上表面**(n̂⁺=+z)          場源 B_voltage_p1.fld
%     'p1cone' = P1 半切件**底錐面**(n̂⁺ 朝下出鋼)      場源 B_voltage_p1.fld
%     'p2cone' = P2 **上錐面**(n̂⁺ 朝上出鋼)            場源 B_voltage_p2.fld
%   每柱撒 NPT 點 → scatteredInterpolant 內插 → ‖b‖[mT] 與對 n̂⁺ 的夾角。
%
%   PAIR = 'p1'   → P1 平切面 vs P1 底錐面（預設；同極換面）
%          'p1p2' → P1 平切面 vs P2 上錐面（實機組合）
%          'cone' → P1 底錐面 vs P2 上錐面（同面換極：控制貼附面、只看上/下磁極差異）
%   橫軸 = ‖b‖ 大小值 (mT)、縱軸 = 百分比 (%)。
%
%   ⚠ 這是**內插**(Maxwell 格距 0.1mm)，與校正管線 V_METHOD='scattered' 同一套。
%     已驗證 500 點的內插四面體**無一含鋼件頂點**（撒點離最近鋼格點 0.385~0.524mm）。
%   ⚠ [MODIFIED 2026-08-10] sensor 幾何 = matlab/Maxwell/utils/pole_sensor_geometry.m（**唯一來源**，
%     與校正管線 build_V_matrix 同一支）→ **不再需要兩邊同步**；原本此處複製的舊公式已刪除。
%   ★ 本檔的 per-pole 快取 `sensor_B_hist_P<k>_maxwell_soff<S>_n<N>.mat` 同時被
%     plot_sensor_cyl_B_3d.m 讀取（存撒點 P_<case>、場向量 B_<case>、sensor 鄰域格點 Xs/Bs）。
%
%   plot_sensor_B_hist()                    % P1 兩面，有快取就用
%   plot_sensor_B_hist([],[],[],[], 'p1p2') % P1 平切面 vs P2 上錐面
%   plot_sensor_B_hist(true)                % 強制重讀 2GB .fld 重採樣(慢)
%
%   風格①粗體框圖 + 疊圖直方圖慣例(共用 edges、mean 虛線、百分比縱軸、圖例在框外)。
% =========================================================================
    clc;
    if nargin < 1 || isempty(force), force = false;    end
    if nargin < 2 || isempty(NPT),   NPT   = 500;      end
    if nargin < 3 || isempty(SOFF),  SOFF  = 4.572e-3; end
    if nargin < 4 || isempty(NB),    NB    = 180;      end
    if nargin < 5 || isempty(PAIR),  PAIR  = 'p1';     end

    here   = fileparts(fileparts(mfilename('fullpath')));      % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    % SHOW_CV：圖例的 mean 條目要不要附 CV。
    %   'p1'（同一顆極換面）→ 附，均勻度差 1.55 倍是該圖的重點之一。
    %   'p1p2'（跨極比較）→ 不附（使用者拍板 2026-08-06）。
    switch lower(PAIR)
        case 'p1',   cs = {'p1flat','p1cone'};   lbl = {'Flat cut surface','Bottom cone surface'};
                     tag = 'P1_flat_vs_cone';    SHOW_CV = true;
        case 'p1p2', cs = {'p1flat','p2cone'};   lbl = {'P1 flat cut surface','P2 cone surface'};
                     tag = 'P1flat_vs_P2cone';   SHOW_CV = false;
        % [ADDED 2026-08-06] 同為錐面的跨極比較——把「貼附面」控制住,只看「上/下磁極」差異。
        case 'cone', cs = {'p1cone','p2cone'};   lbl = {'P1 cone surface','P2 cone surface'};
                     tag = 'P1cone_vs_P2cone';   SHOW_CV = false;
        otherwise,   error('PAIR 必為 ''p1'' | ''p1p2'' | ''cone''');
    end

    P1 = get_case(cs{1}, here, NPT, SOFF, force);
    P2 = get_case(cs{2}, here, NPT, SOFF, force);
    b1 = P1.b;   b2 = P2.b;                     % ‖b‖ [mT]

    % ---- 統計（含對 n̂⁺ 的夾角）----
    for k = 1:2
        S = {P1,P2}; S = S{k};
        fprintf(['[%-7s] n=%d  ‖b‖ mean=%.4f  min=%.4f  max=%.4f mT  σ/μ=%.2f%%  ' ...
                 '|  ⟨b·n̂⁺⟩=%.4f mT  對 n̂⁺ 夾角 mean=%.2f° (min %.2f° / max %.2f° / σ %.2f°)\n'], ...
                cs{k}, numel(S.b), mean(S.b), min(S.b), max(S.b), std(S.b)/mean(S.b)*100, ...
                mean(S.bn), mean(S.ang), min(S.ang), max(S.ang), std(S.ang));
    end
    fprintf('%s / %s：‖b‖ mean 比 = %+.2f%%   分布重疊 = %s\n', cs{1}, cs{2}, ...
            (mean(b1)/mean(b2)-1)*100, ternary_str(min(b1) > max(b2) || min(b2) > max(b1), '無(完全分離)', '有'));

    % ---- 共用 edges(兩組畫在同一 x 軸、同一 bin 寬)----
    lo = min([b1; b2]);  hi = max([b1; b2]);
    edg = linspace(lo, hi, NB+1);   ctr = (edg(1:end-1)+edg(2:end))/2;
    pct1 = histcounts(b1, edg) / numel(b1) * 100;
    pct2 = histcounts(b2, edg) / numel(b2) * 100;
    mu1 = mean(b1);  mu2 = mean(b2);

    % ---- 繪圖(風格①；配色沿用 R150 組:藍/紅)----
    cB = [0.10 0.35 1.00];   cR = [0.85 0.10 0.10];
    FS = 28;
    fig = figure('Color','w','Position',[100 100 1100 830]);  ax = axes(fig);  hold(ax,'on');
    h1 = bar(ax, ctr, pct1, 1, 'FaceColor',cB, 'FaceAlpha',0.60, 'EdgeColor','k', 'LineWidth',0.3);
    h2 = bar(ax, ctr, pct2, 1, 'FaceColor',cR, 'FaceAlpha',0.60, 'EdgeColor','k', 'LineWidth',0.3);
    ml1 = xline(ax, mu1, '--', 'Color',[0.00 0.00 0.00], 'LineWidth',2.8);
    ml2 = xline(ax, mu2, '--', 'Color',[0.00 0.60 0.00], 'LineWidth',2.8);

    box(ax,'on');  grid(ax,'off');
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.015 .015],'TickDir','out');
    [xr_, xt_] = xticks_odd(lo, hi);        % 奇數個等距 tick + 兩端留白 = 間距
    xlim(ax, xr_);  set(ax,'XTick', xt_);
    ymax = max([pct1 pct2]);  ytop = ceil(ymax*10)/10;  ylim(ax,[0 ytop]);
    yd = round(linspace(0, ytop, 5), 1);  set(ax,'YTick', yd(2:4));   % 3 內縮 tick(首末不標)
    % 水平軸起點 + 終點:只標數字、不畫 tick(左右角落補字)
    xr = xlim(ax);
    text(ax, xr(1), -0.022*ytop, sprintf('%g',xr(1)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');
    text(ax, xr(2), -0.022*ytop, sprintf('%g',xr(2)), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS,'FontWeight','bold','Clipping','off');

    cv1 = std(b1)/mu1*100;   cv2 = std(b2)/mu2*100;      % CV = σ/μ [%]（不受 bin 寬影響）
    % ⚠ 不可用「同一個 fmt + 永遠傳 cv」的寫法：sprintf 引數多於格式符時會**循環套用格式**，
    %   會在圖例後面接出一串亂碼（踩過）。要嘛連引數一起省掉。
    % [MODIFIED 2026-08-15] 圖例改成「**每個系列一行、標記一律方框**」（使用者拍板）：
    %   統計值（mean / CV）併進該系列自己那一行，不再另外用虛線標記佔一欄。
    %   mean 的虛線仍畫在圖上（ml1/ml2），只是不進圖例。
    if SHOW_CV
        s1 = sprintf('%s: mean = %.3f mT, CV = %.2f%%', lbl{1}, mu1, cv1);
        s2 = sprintf('%s: mean = %.3f mT, CV = %.2f%%', lbl{2}, mu2, cv2);
    else
        s1 = sprintf('%s: mean = %.3f mT', lbl{1}, mu1);
        s2 = sprintf('%s: mean = %.3f mT', lbl{2}, mu2);
    end
    lg = legend([h1 h2], {s1, s2}, ...
                'Interpreter','tex', 'Location','northoutside', 'NumColumns',1);
    lg.FontSize = 24;  lg.FontWeight = 'bold';  lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = 2.5;
    xlabel(ax, '$\mathbf{\|b\|\;(mT)}$', 'Interpreter','latex', 'FontSize',36);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$', 'Interpreter','latex', 'FontSize',36);
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    drawnow;
    axp = get(ax,'Position');  lgh = get(lg,'Position');  lgh = lgh(4);
    GAPN = 0.022;  newTop = 1 - lgh - GAPN - 0.006;
    axp(4) = newTop - axp(2);  set(ax,'Position',axp);
    set(lg, 'Position', [axp(1), newTop + GAPN, axp(3), lgh]);

    out = fullfile(figdir, sprintf('sensor_B_hist_%s_maxwell_soff%g_n%d.png', tag, SOFF*1e3, NPT));
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function s = ternary_str(c,a,b), if c, s=a; else, s=b; end, end
function w = first_word(s), t = strsplit(s,' '); w = t{1}; end

% ============================================================================
function [C, S] = get_case(name, here, NPT, SOFF, force)
% 取出單一貼附面的撒點結果 C（.P[mm,WP框] .B[mT] .b .bn .ang .c .n），S = 該極的完整快取。
    kpole = 1;  if strncmpi(name,'p2',2), kpole = 2; end
    S = load_pole(kpole, here, NPT, SOFF, force);
    C.P = S.(['P_' name]);   C.B = S.(['B_' name]);
    C.c = S.(['c_' name]);   C.n = S.(['n_' name]);
    C.b  = vecnorm(C.B,2,2);                       % ‖b‖ [mT]
    C.bn = C.B * C.n;                              % b·n̂⁺ [mT]（投影，V 就是它的平均×S_hall）
    C.ang = acosd( min(1,max(-1, C.bn ./ C.b)) );  % 對 n̂⁺ 夾角 [deg]
end

% ============================================================================
function S = load_pole(kpole, here, NPT, SOFF, force)
% per-pole 快取：撒點 + 場向量 + sensor 鄰域原始格點（改樣式/加面都不必重讀 2GB）。
    CACHE = fullfile(here,'data', sprintf('sensor_B_hist_P%d_maxwell_soff%g_n%d.mat', kpole, SOFF*1e3, NPT));
    if ~force && exist(CACHE,'file')
        S = load(CACHE);   fprintf('cache hit: %s\n', CACHE);
        if ~isfield(S,'has_vec')        % 舊版快取只有 ‖b‖ → 由快取內的場源補算向量（不讀 .fld）
            S = add_vectors(S, kpole);  save(CACHE,'-struct','S');
            fprintf('  (舊快取已補上撒點/向量欄位並回存)\n');
        end
    else
        S = sample_pole(kpole, NPT, SOFF);
        save(CACHE,'-struct','S');   fprintf('wrote cache %s\n', CACHE);
    end
end

% ============================================================================
function S = sample_pole(kpole, NPT, SOFF)
% 讀該極自激發的 sensor 區 .fld → 該極所有貼附面各撒 NPT 點 → 內插。
    cfg = mw_cfg();
    faces = pole_faces(kpole);
    % sensor 幾何（WP 框）→ ANSYS/Maxwell 框
    cen = struct();  nrm = struct();
    for i = 1:numel(faces)
        [p, n] = sensor_geom(cfg, kpole, faces{i}, SOFF);
        cen.(faces{i}) = p + [0;0;cfg.SPH_OFST];   nrm.(faces{i}) = n;
    end
    fp = fullfile(cfg.fld_dir, cfg.fld_files_voltage{kpole});   % Maxwell map = identity → p<k> = P<k>
    fprintf('reading %s (~2 GB, 需數分鐘) ...\n', fp);
    d = import_maxwell_fld(fp);
    X = [d.x, d.y, d.z];   B = [d.bx, d.by, d.bz];
    R_loc = 1.5e-3;  if isfield(cfg,'sensor_r_loc') && ~isempty(cfg.sensor_r_loc), R_loc = cfg.sensor_r_loc; end
    keep = false(size(X,1),1);
    for i = 1:numel(faces), keep = keep | vecnorm(X - cen.(faces{i}).',2,2) < R_loc; end
    S.Xs = X(keep,:);   S.Bs = 1e3*B(keep,:);      % [mT]
    fprintf('sensor 鄰域格點 %d / 全域 %d\n', nnz(keep), size(X,1));
    S.NPT = NPT;  S.SOFF = SOFF;  S.R_loc = R_loc;  S.src = fp;  S.kpole = kpole;  S.faces = faces;
    for i = 1:numel(faces)
        S.(['c_' faces{i}]) = cen.(faces{i});   S.(['n_' faces{i}]) = nrm.(faces{i});
    end
    S = add_vectors(S, kpole);
end

% ============================================================================
function S = add_vectors(S, kpole)
% 由快取內的 sensor 鄰域格點重建撒點與場向量（不必讀 .fld）。
%   撒點序列固定：rng(0) 後依 faces 順序逐面抽 → 同一顆極的結果可重現。
    cfg = mw_cfg();
    if ~isfield(S,'faces') || isempty(S.faces)      % 舊版快取（P1、欄位名為 cF/nF/cC/nC）
        S.faces = pole_faces(kpole);
        S.c_p1flat = S.cF;  S.n_p1flat = S.nF;
        S.c_p1cone = S.cC;  S.n_p1cone = S.nC;
    end
    FX = scatteredInterpolant(S.Xs(:,1),S.Xs(:,2),S.Xs(:,3),S.Bs(:,1),'linear','none');
    FY = scatteredInterpolant(S.Xs(:,1),S.Xs(:,2),S.Xs(:,3),S.Bs(:,2),'linear','none');
    FZ = scatteredInterpolant(S.Xs(:,1),S.Xs(:,2),S.Xs(:,3),S.Bs(:,3),'linear','none');
    rng(0);
    for i = 1:numel(S.faces)
        nm = S.faces{i};
        [P, Bv] = one_face(S.(['c_' nm]), S.(['n_' nm]), S.NPT, FX,FY,FZ, cfg.SPH_OFST);
        S.(['P_' nm]) = P;   S.(['B_' nm]) = Bv;
    end
    S.has_vec = true;
end

function [Pmm, Bv] = one_face(ci, ni, NPT, FX,FY,FZ, SPH_OFST)
% 圓柱內均勻撒點（同 build_V_matrix）→ 線性內插；回傳 WP 框 [mm] 與 B [mT]。
    sensor_r = 0.15e-3;  axial_tol = 0.10e-3;
    t1 = [-ni(2); ni(1); 0];  if norm(t1) < 1e-9, t1 = [1;0;0]; end
    t1 = t1/norm(t1);  t2 = cross(ni, t1);
    a  = axial_tol * rand(NPT,1);
    rr = sensor_r  * sqrt(rand(NPT,1));            % √U → 截面積均勻
    th = 2*pi      * rand(NPT,1);
    P  = ci.' + a.*ni.' + (rr.*cos(th)).*t1.' + (rr.*sin(th)).*t2.';   % ANSYS 框 [m]
    Bv = [FX(P), FY(P), FZ(P)];                                        % [mT]
    if any(isnan(Bv(:))), warning('one_face: %d 點落在凸包外(NaN)', nnz(any(isnan(Bv),2))); end
    Pmm = P;  Pmm(:,3) = Pmm(:,3) - SPH_OFST;   Pmm = Pmm * 1e3;       % → WP 框 [mm]
end

% ============================================================================
function f = pole_faces(kpole)
    if kpole == 1, f = {'p1flat','p1cone'}; else, f = {'p2cone'}; end
end

% ============================================================================
function cfg = mw_cfg()
% Maxwell 分支 config（rmpath APDL 分支防同名函式遮蔽）。
    MW   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    APDL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(APDL,'function'));  rmpath(fullfile(APDL,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(MW,'function'));   addpath(fullfile(MW,'common_path'));
    addpath(fullfile(MW,'utils'));      % [ADDED 2026-08-08] pole_sensor_geometry
    cfg = model_config('long2016_hexapole_halfcut','tip40um');
end

% ============================================================================
function [pos, nhat] = sensor_geom(cfg, kpole, face, soff)
% Hall 中心 + 法線 n̂⁺（WP 框）。
% [MODIFIED 2026-08-08] 改呼叫 utils/pole_sensor_geometry（**sensor 幾何唯一來源**）。
%   原本這裡複製了一份 build_V_matrix 的舊公式（名目 beta、無真實錐體半徑外推），已刪除。
    switch lower(face)
        case 'p1flat',            fl = 'flat';   % 下極平切上表面
        case {'p1cone','p2cone'}, fl = 'cone';   % 下極底錐面 / 上極錐面（上極不受 face_lower 影響）
        otherwise, error('sensor_geom: 未知 face ''%s''', face);
    end
    [P, N] = pole_sensor_geometry(cfg, struct('soff_upper', soff, ...
                                              'soff_lower', soff, 'face_lower', fl));
    pos = P(:,kpole);  nhat = N(:,kpole);
end

% ============================================================================
function [xr, xt] = xticks_odd(lo, hi)
% 橫軸刻度:**奇數個(5)等距** nice 步長，且**兩端留白 = tick 間距**(figure-style 2D 規則)。
    cand = [1 2 2.5 3 4 5 10];
    span = hi - lo;
    k0 = floor(log10(max(span,realmin))) - 1;
    for k = k0:(k0+3)
        for c = cand
            s = c*10^k;
            ctr0 = round((lo+hi)/2 / s) * s;
            t = ctr0 + (-2:2)*s;
            if t(1)-s <= lo && t(end)+s >= hi
                xr = [t(1)-s, t(end)+s];  xt = t;  return;
            end
        end
    end
    s = span/5;  t = (lo+hi)/2 + (-2:2)*s;  xr = [t(1)-s, t(end)+s];  xt = t;
end
