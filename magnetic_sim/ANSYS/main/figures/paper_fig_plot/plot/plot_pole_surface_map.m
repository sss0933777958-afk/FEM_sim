function plot_pole_surface_map(MODE, DELTA, force)
% plot_pole_surface_map -- 磁極側表面「攤平（展開）」熱圖：b·n̂ 分布（Maxwell）
% =========================================================================
%   錐面沿母線展開（isometric development）→ 以**虛擬錐頂**為心的環狀扇形：
%     ρ = σ/cosβ（σ = 自虛擬錐頂的軸向距離）、ψ = φ·sinβ
%     P2 完整錐  張角 2π·sinβ = 68.8°
%     P1 半錐    張角 π·sinβ = 35.9° ＋ 削平面（本身即平面）2β = 23.0° → 共 58.9°
%   展開為**等距映射** → 圖上面積 = 真實側表面積（總面積可核 105.40 / 124.22 mm²），
%   故格心值的面積加權平均可直接對回 ΔΦ/A = 11.33 / 10.84 mT。
%   ⚠ 展開必須在某條母線「剪開」：P1 剪在**錐面正下方**（φ=180°）→ 那條線在圖上
%     出現在扇形的**兩個外緣**（= P1 錐面 sensor 的方位）。P2 剪在正下方（φ=±180°）。
%
%   幾何（STEP 實測 per-pole，見 rule ansys-cad-alignment）：
%     下極 R(s) = 0.0327 + 0.20330·s（β = 11.492°）  上極 0.0330 + 0.19463·s（11.014°）
%     虛擬錐頂在極尖點外 e = R0/slope ≈ 0.16 mm；σ = s + e。
%   P1 削平面 = **通過極軸的水平面**（mt_constants 下極 pole_axis 水平、sensor_geometry
%   flat 分支法線 +z）→ 截面為正半圓、平面寬 2R(s)、鋼件在平面下方。
%
%   取樣：站位 s = 1.0~14.0（步長 DS）、周向格邊 ≈ DL；每格格心沿**外法線**位移
%   δ（DELTA 參數）後，由**空氣側**格點 3D 內插（'linear','none'）取 b·n̂。
%     **δ = 0.41 mm（預設）= sensor 實際氣隙** → 圖上的值即「sensor 若裝在該處會看到的場」。
%       ⚠ Hall 讀值還要再對 R0.15×H0.10mm 圓柱做體積平均，本圖是單點。
%     δ = 0.15 mm = 貼面版（0.1mm 規則格能支撐真實內插的最近距離；δ≤0.05 被邊界四面體壓平），
%       用來與 ΔΦ/A 做通量核對。兩版共存於 data/（快取檔名含 d015 / d041）。
%     空氣側判準 = 幾何（離鋼面 > 0.05mm）+ ‖B‖ < 400 mT 護欄。
%     ⚠ 不可用舊的 ‖B‖<40mT 判準：s≈1 附近**空氣側**本來就有 50+ mT，會被誤刪。
%   ψ=0（黑虛線 = 中心線）= 該極 sensor 方位：P1 削平面中心（正上方）、P2 錐面正上方。
%   P1 另以灰虛線標**切邊**（|ψ| = β）＝ 削平面與錐面的交界。
%
%   ⚠ 這是**內插**（Maxwell 匯出 0.1 mm 規則格）。⚠ s < 1 不做：近錐頂格點太稀。
%
%   風格①粗體框圖；兩極**共用 clim 與座標尺度**（figure-style「同類比較共用色階」）。
% =========================================================================
    clc;
    if nargin < 1 || isempty(MODE),  MODE  = 'abs'; end     % 'abs' | 'rel'
    if nargin < 2 || isempty(DELTA), DELTA = 0.41;  end     % 離鋼面距離 [mm]
    if nargin < 3 || isempty(force), force = false; end

    here   = fileparts(fileparts(mfilename('fullpath')));      % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    M = cell(1,2);
    for k = 1:2
        cf = fullfile(here,'data', sprintf('surfmap_P%d_d%03d_maxwell.mat', k, round(DELTA*100)));
        if exist(cf,'file') && ~force
            M{k} = load(cf);
            fprintf('P%d: 由快取載入 %s（%d 格）\n', k, cf, numel(M{k}.val));
        else
            S = build_map(k, here, DELTA);   save(cf,'-struct','S');   M{k} = S;
            fprintf('P%d: 已存 %s\n', k, cf);
        end
    end

    % ---- sensor 貼附點（貼在表面上，不含 0.41mm 氣隙）→ 展開座標 ----
    %   sensor_geometry: 貼附點 = tip + SOFF·e2。P1 削平面 e2 沿極軸 → 軸向 t = SOFF；
    %   P2 錐面 e2 沿母線 → t = SOFF·cos(β_nom)（β_nom = atan(3/15)，config 用的名目值）。
    %   兩者都在 ψ=0 的中心線上。
    SOFF = 4.572;   bn = atan2(3,15);
    for k = 1:2
        bet = M{k}.beta_deg*pi/180;
        if k == 1, M{k}.sens_x = SOFF + M{k}.e_apex;                       % 平面：x = σ
        else,      M{k}.sens_x = (SOFF*cos(bn) + M{k}.e_apex)/cos(bet);    % 錐面：ρ = σ/cosβ
        end
        M{k}.sens_y = 0;
    end
    % P1 的另一個候選貼附面 = **底錐面**（FACE_lower='cone'）：在錐面正下方 φ=180°，
    % 那正是展開的剪開線 → 圖上出現在扇形的**兩個外緣**（同一物理點，被剪成兩半）。
    b1 = M{1}.beta_deg*pi/180;   r1 = (SOFF*cos(bn) + M{1}.e_apex)/cos(b1);
    p1 = b1 + (pi/2)*sin(b1);
    M{1}.sens2_x = r1*cos(p1)*[1 1];   M{1}.sens2_y = r1*sin(p1)*[1 -1];

    A_ref = [105.40 124.22];                       % 真實側表面積（trapz，s=1~14）
    B_ref = [11.33 10.84];                         % ΔΦ/A（鋼面）[mT]
    for k = 1:2
        m = M{k};  w = m.area;  v = m.val;  ok = isfinite(v);
        mu = sum(v(ok).*w(ok))/sum(w(ok));
        fprintf(['[P%d] δ=%.2fmm｜格 %d（NaN %d）｜展開面積 %.2f mm²（真實 %.2f, 差 %+.2f%%）｜' ...
                 'b·n̂ 面積加權平均 %.2f mT（鋼面 ΔΦ/A = %.2f → 比值 %.2f）｜全距 %.1f ~ %.1f mT\n'], ...
                 k, m.delta, numel(v), sum(~ok), sum(w), A_ref(k), (sum(w)/A_ref(k)-1)*100, ...
                 mu, B_ref(k), B_ref(k)/mu, min(v(ok)), max(v(ok)));
    end

    % ---- 相對模式：除以「該站位的面積加權環平均」→ 只留方位結構 ----
    if strcmpi(MODE,'rel')
        for k = 1:2
            m = M{k};  [st, ~] = station_of(m);
            for z = unique(st).'
                q = st == z;
                m.val(q) = m.val(q) / (sum(m.val(q).*m.area(q))/sum(m.area(q)));
            end
            M{k} = m;
            fprintf('[P%d] 相對值 min %.3f / max %.3f（>1.6 的格 %.1f%% 會飽和）\n', ...
                    k, min(m.val), max(m.val), 100*mean(m.val>1.6));
        end
        % 比值圖用**發散色階**（中心 = 1）：turbo 中段是綠→黃、對比最弱，而 1.0 正好
        % 落在綠色平台 → ±30% 的起伏看不出來（踩過）。發散色階在中心對比最強。
        CLIM = [0.6 1.4];   CTK = 0.6:0.2:1.4;   CMAP = cmap_div(256);
        CLAB = '$\mathbf{b\cdot\hat{n}\;/\;\langle b\cdot\hat{n}\rangle_{ring}}$';
        SFX  = '_rel';
        for k = 1:2
            fprintf('[P%d] 超出色階（<0.6 或 >1.4）的格 %.1f%%\n', k, ...
                    100*mean(M{k}.val < CLIM(1) | M{k}.val > CLIM(2)));
        end
    else
        hi   = max(cellfun(@(m) prctile(m.val(isfinite(m.val)),99.5), M));
        CLIM = [0, ceil(hi/10)*10];   CTK = linspace(CLIM(1),CLIM(2),5);  CMAP = turbo(256);
        CLAB = '$\mathbf{b\cdot\hat{n}\;(mT)}$';
        SFX  = '';
    end

    % 框由 tick 反推：3 個等距 tick + 兩端留白 = 間距（figure-style 2D）
    sx = fit_step(max(cellfun(@(m) max(m.VX), M))/4);      XL = [0, 4*sx];
    sy = fit_step(max(cellfun(@(m) max(abs(m.VY)), M))/2); YL = [-2*sy, 2*sy];
    fprintf('共用 clim = [%g %g] ; XL = [%g %g] ; YL = [%g %g]\n', CLIM, XL, YL);

    for k = 1:2, render(M{k}, k, CLIM, CTK, CLAB, CMAP, SFX, XL, YL, [sx sy], figdir); end
end

% ============================================================================
function C = cmap_div(n)
% 發散色階（coolwarm 風格）：藍 — 淺灰 — 紅，中心對比最強
    a = [0.230 0.299 0.754];  b = [0.865 0.865 0.865];  c = [0.706 0.016 0.150];
    h = floor(n/2);
    t1 = linspace(0,1,h).';    t2 = linspace(0,1,n-h).';
    C  = [a + t1.*(b-a); b + t2.*(c-b)];
end

% ============================================================================
function [st, psi] = station_of(m)
% 由展開座標還原每格的站位 s 與展開角 psi
    bet = m.beta_deg*pi/180;
    VX  = reshape(m.VX,4,[]).';   VY = reshape(m.VY,4,[]).';
    psi = mean(atan2(VY,VX),2);   rho = mean(hypot(VX,VY),2);
    isf = abs(psi) <= bet + 1e-9;                       % 削平面 vs 錐面
    sg  = rho*cos(bet);   sg(isf) = rho(isf).*cos(psi(isf));
    % ⚠ 站位中心落在 DS 的半格上 → 直接 round(σ/DS) 會踩在 .5 的分界、同一站被切成兩組
    %   （相對圖會出現縱向條紋）。先減半格再取整，落點才在整數上、穩定。
    st  = round((sg - m.e_apex - m.DS/2)/m.DS);
end

% ============================================================================
function render(m, k, CLIM, CTK, CLAB, CMAP, SFX, XL, YL, stp, figdir)
    FS = 36;
    fig = figure('Color','w','Position',[60 40 1500 1000]);  ax = axes(fig);  hold(ax,'on');

    patch(ax, 'Faces',m.F, 'Vertices',[m.VX m.VY], 'FaceVertexCData',m.val, ...
          'FaceColor','flat', 'EdgeColor','none');
    colormap(ax, CMAP);   clim(ax, CLIM);

    rr = [min(m.rho) max(m.rho)];                             % 中心線 = sensor 方位
    plot(ax, rr, [0 0], '--', 'Color','k', 'LineWidth',3.0);
    if k == 1                                                 % P1 切邊
        for sg = [-1 1]
            plot(ax, rr*cos(m.psi_edge), sg*rr*sin(m.psi_edge), '--', ...
                 'Color',[0.30 0.30 0.30], 'LineWidth',2.2);
        end
    end

    % ---- sensor 安裝位置（SOFF 4.572mm）----
    %   圓 = 現用貼附面（P1 削平面中心 / P2 錐面正上方，皆在中心線上）
    %   方 = P1 的另一候選「底錐面」：在剪開線上 → 兩個外緣各一（同一物理點）
    if isfield(m,'sens2_x')
        plot(ax, m.sens2_x, m.sens2_y, 's', 'MarkerSize',22, 'LineWidth',4, ...
             'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'Clipping','off');
    end
    plot(ax, m.sens_x, m.sens_y, 'o', 'MarkerSize',20, 'LineWidth',4, ...
         'MarkerFaceColor','w', 'MarkerEdgeColor','k');

    box(ax,'on');  grid(ax,'off');  daspect(ax,[1 1 1]);
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',2.5,'TickLength',[.012 .012],'TickDir','out');
    xlim(ax, XL);  ylim(ax, YL);
    set(ax,'XTick', XL(1)+(1:3)*stp(1), 'YTick', (-1:1)*stp(2));

    yoff = YL(1) - 0.025*diff(YL);                 % 水平軸端點只標數字、不畫 tick
    for j = 1:2
        text(ax, XL(j), yoff, sprintf('%g',XL(j)), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end

    cb = colorbar(ax);  cb.FontSize = FS;  cb.FontWeight = 'bold';
    cb.Label.Interpreter = 'latex';                % Interpreter 要先於 String
    cb.Label.String = CLAB;
    cb.Label.FontSize = FS;
    cb.Ticks = CTK;

    xlabel(ax, '$\mathbf{x\;(mm)}$', 'Interpreter','latex', 'FontSize',FS);
    ylabel(ax, '$\mathbf{y\;(mm)}$', 'Interpreter','latex', 'FontSize',FS);
    ax.Toolbar.Visible = 'off';  hold(ax,'off');

    out = fullfile(figdir, sprintf('surface_map_P%d%s_maxwell.png', k, SFX));
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function m = build_map(k, here, DELTA)
    DS = 0.20;  DL = 0.20;  S0 = 1.0;  S1 = 14.0;

    MW = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(MW,'function'));
    cnst = model_config('long2016_hexapole_halfcut','tip40um');
    G    = load(fullfile(here,'data','flux_vs_s_maxwell.mat'),'GEO');
    R0   = G.GEO(k).R0;  sl = G.GEO(k).sl;          % R(s) = R0 + sl·s  [mm]
    bet  = atan(sl);     e_ap = R0/sl;              % 半錐角、虛擬錐頂偏移
    low  = logical(cnst.pole_is_lower(k));

    % ---- 局部框：â 沿極軸(尖→根)、û 為含 ẑ 的上方向、v̂ = â×û ----
    T  = [cnst.pole_tip_x(k); cnst.pole_tip_y(k); cnst.pole_tip_z_wp(k)]*1e3;
    ah = cnst.pole_axis(:,k);  ah = ah/norm(ah);
    uh = [0;0;1] - dot([0;0;1],ah)*ah;  uh = uh/norm(uh);
    vh = cross(ah, uh);
    Rl = [ah, uh, vh].';

    fprintf('P%d: 讀 %s ...\n', k, cnst.fld_files_voltage{k});
    d  = import_maxwell_fld(fullfile(cnst.fld_dir, cnst.fld_files_voltage{k}));
    Pw = [d.x*1e3, d.y*1e3, d.z*1e3 - cnst.SPH_OFST*1e3];
    Bw = cnst.s_source(k) * [d.bx, d.by, d.bz] * 1e3;              % all-source, mT
    clear d;
    Pl = (Pw - T.')*Rl.';     Bl = Bw*Rl.';       clear Pw Bw;

    t = Pl(:,1);  u = Pl(:,2);  v = Pl(:,3);  r = hypot(u,v);  Rt = R0 + sl*t;
    near = t > 0.5 & t < 14.6 & r < Rt + 0.9;
    air  = r > Rt + 0.05;
    if low, near = near & u < 0.9;   air = air | (u > 0.05);  end
    keep = near & air & vecnorm(Bl,2,2) < 400;
    fprintf('   格點 %d → 殼層空氣側 %d\n', numel(t), sum(keep));
    Ps = Pl(keep,:);  Bs = Bl(keep,:);
    clear Pl Bl t u v r Rt near air keep;

    Ix = scatteredInterpolant(Ps, Bs(:,1), 'linear','none');
    Iy = scatteredInterpolant(Ps, Bs(:,2), 'linear','none');
    Iz = scatteredInterpolant(Ps, Bs(:,3), 'linear','none');

    sed = S0:DS:S1;   nst = numel(sed)-1;
    Q = zeros(0,3);  N = zeros(0,3);  AR = zeros(0,1);
    PS = zeros(0,2);  SG = zeros(0,2);  FLT = false(0,1);
    for i = 1:nst
        s1 = sed(i);  s2 = sed(i+1);  sc = (s1+s2)/2;
        Rc = R0 + sl*sc;   sg1 = s1 + e_ap;  sg2 = s2 + e_ap;  sgc = sc + e_ap;
        if low                                                  % ---- 削平面 ----
            nf = 2*max(1, round(Rc/DL)) + 1;                    % 奇數 → ψ=0 為格心
            ve = linspace(-Rc, Rc, nf+1);
            for j = 1:nf
                vc = (ve(j)+ve(j+1))/2;
                Q(end+1,:)  = [sc, DELTA, vc];                  %#ok<AGROW>  外法線 = +û
                N(end+1,:)  = [0, 1, 0];                        %#ok<AGROW>
                AR(end+1,1) = (ve(j+1)-ve(j))*DS;               %#ok<AGROW>
                PS(end+1,:) = [atan(ve(j)/sgc), atan(ve(j+1)/sgc)];  %#ok<AGROW>
                SG(end+1,:) = [sg1, sg2];                       %#ok<AGROW>
                FLT(end+1,1)= true;                             %#ok<AGROW>
            end
            phi0 = pi/2;  phi1 = 3*pi/2;                        % 半錐（鋼側 u<0）
            nc   = 2*max(2, round((pi*Rc/2)/DL));               % **偶數**：剪開線 φ=π 落在格界
        else
            phi0 = -pi;   phi1 = pi;                            % 完整錐（剪開線在 ±π）
            nc   = 2*max(2, round((pi*Rc)/DL)) + 1;             % **奇數**：ψ=0 為格心
        end
        pe = linspace(phi0, phi1, nc+1);
        for j = 1:nc
            pc = (pe(j)+pe(j+1))/2;
            er = [0; cos(pc); sin(pc)];
            nn = cos(bet)*er - sin(bet)*[1;0;0];
            Q(end+1,:)  = ([sc;0;0] + Rc*er + DELTA*nn).';      %#ok<AGROW>
            N(end+1,:)  = nn.';                                 %#ok<AGROW>
            AR(end+1,1) = Rc*(pe(j+1)-pe(j))*DS/cos(bet);       %#ok<AGROW>
            PS(end+1,:) = [dev_psi(pe(j),bet,low,pc), dev_psi(pe(j+1),bet,low,pc)];  %#ok<AGROW>
            SG(end+1,:) = [sg1, sg2];                           %#ok<AGROW>
            FLT(end+1,1)= false;                                %#ok<AGROW>
        end
    end

    B   = [Ix(Q), Iy(Q), Iz(Q)];
    val = sum(B.*N, 2);                                          % b·n̂ [mT]
    fprintf('   查詢 %d 格，NaN %d\n', size(Q,1), sum(~isfinite(val)));

    % ---- 展開座標：flat ρ = σ/cosψ（等效 x=σ）；cone ρ = σ/cosβ ----
    nQ = size(Q,1);   VX = zeros(nQ,4);  VY = zeros(nQ,4);  RH = zeros(nQ,1);
    for q = 1:nQ
        ps = [PS(q,1) PS(q,2) PS(q,2) PS(q,1)];
        sg = [SG(q,1) SG(q,1) SG(q,2) SG(q,2)];
        if FLT(q), rho = sg./cos(ps);  else, rho = sg/cos(bet); end
        VX(q,:) = rho.*cos(ps);   VY(q,:) = rho.*sin(ps);   RH(q) = mean(rho);
    end

    m = struct('val',val, 'area',AR, 'rho',RH, 'psi_edge',bet, ...
               'VX', reshape(VX.',[],1), 'VY', reshape(VY.',[],1), ...
               'F',  reshape(1:4*nQ, 4, nQ).', ...
               'DS',DS, 'DL',DL, 'delta',DELTA, 'R0',R0, 'slope',sl, ...
               'beta_deg',bet*180/pi, 'e_apex',e_ap, 'pole',k, 'SRC','maxwell', ...
               'note', ['Isometric development. cone: rho=sigma/cos(beta), psi=phi*sin(beta); ' ...
                        'flat: x=sigma, y=sigma*tan(psi). Cells laid on STEEL surface; ' ...
                        'b.n sampled at delta=0.15mm outward, air-side 3D linear interp. ' ...
                        'Cut line for development = cone bottom (phi=180deg).']);
end

% ============================================================================
function psi = dev_psi(phi, bet, low, phic)
% 錐面方位 φ → 展開角 ψ。下極半錐（φ∈[90°,270°]）接在削平面（|ψ|≤β）外側。
%   ⚠ 剪開線 φ=180° 在展開圖上是**雙值的**（同一條母線出現在扇形兩個外緣）。
%     分支必須由**該格的中點 phic** 決定，不可由角點自己的 φ 決定 ——
%     否則跨在剪開線兩側的那一格，兩個角點會被分到相反外緣，畫出橫跨整個扇形的
%     巨大四邊形（每站位一個，疊成垂直長條把真實結構蓋掉）。2026-08-07 踩過。
    if low
        d = phi - pi/2;                          % 0..π，自 +v 切邊起算
        if phic <= pi,  psi =  bet + d*sin(bet);          % +v 側
        else,           psi = -bet - (pi - d)*sin(bet);   % -v 側
        end
    else
        psi = phi*sin(bet);
    end
end

% ============================================================================
function s = fit_step(x)
% 最小的 nice 步長 s ≥ x（→ 框寬 = n·s 必涵蓋資料，且 tick 等距、兩端留白 = s）
    kk = floor(log10(max(x, realmin)));  cand = [1 2 2.5 4 5 10];
    for e = [kk kk+1]
        for c = cand
            s = c*10^e;  if s >= x - 1e-12, return; end
        end
    end
end
