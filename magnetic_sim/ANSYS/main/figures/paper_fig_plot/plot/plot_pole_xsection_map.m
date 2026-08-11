function plot_pole_xsection_map(S_MM, CELL, MODE, force)
% plot_pole_xsection_map -- 上下磁極截面（垂直極軸）的磁通密度顏色圖（Maxwell，合併雙 panel）
% =========================================================================
%   在**垂直極軸**的平面上切一個截面，平面內鋪 CELL x CELL 的直角小格，**格心取值**，
%   畫 b·â 的顏色圖（â 朝極尖為正，與 plot_flux_vs_s 同號誌）。左 panel = 上極 P2、
%   右 panel = 下極 P1（與 plot_circuit_side 的上左／下右一致），**共用絕對色階**。
%     P1 下極（半切）-> 半圓截面（削平面通過極軸、鋼件在 v <= 0 側）
%     P2 上極（完整）-> 整圓截面
%
%   面內座標：**u = 水平**（切向，垂直子午面）、**v = 垂直**（面內含 +z 的方向）。
%   下極的削平面因此落在 v = 0、鋼件在下半 -> 圖上一眼看出半切。
%
%   站位 s = S_MM（預設 4.572 mm）**沿極軸自極尖點**量。
%   ⚠ 4.572 mm 是 sensor SOFF 的**斜面**距離；若要對齊 sensor 的實際軸向站位，
%     應改用 4.4832 mm（= (SOFF - t_tan)·cos(beta) + t_tan）。本圖採沿極軸解讀。
%   截面半徑 R(s) = (e_v + s)·tan(beta)，e_v = r_f/sin(beta) - r_f（虛擬錐頂外移）。
%   幾何 = CAD STEP 實測 per-pole（下 tan(beta)=0.20330 / 上 0.19463，見 rule
%   ansys-cad-alignment）；R0 = e_v·tan(beta) 由 POLE_TIP_R 導出，不寫死。
%
%   資料源 = **voltage 粗格** .fld（B_voltage_p<k>.fld，步距 0.1 mm）。WP 細格（0.02 mm）
%   框只有 x,y ∈ ±0.6 mm，到不了 s ~ 4.5 mm，故不可用。激發 = **自激發**。
%
%   🔒 **內插源限定鋼件內**（使用者拍板 2026-08-10）：先用錐體幾何把空氣格點剔除，
%   只用鋼內格點建 scatteredInterpolant。否則邊界格心雖在鋼內，trilinear 的模板會
%   伸出鋼外抓到空氣值 —— 實測會把 P1 邊界格由 246 拉到 124 mT、P2 由 156 拉到 22 mT，
%   在圖上長出一圈**完全虛假**的「邊緣衰減」結構。凸包外的格用 nearest 補（截面內部
%   b·â 均勻到 ±1%，nearest 無害），並回報補值格數。
%
%   ⚠ 這是**內插**：截面不落在源格線上，上極截面又傾斜 36.59°。CELL 小於 0.1 mm
%     （源格距）只是把同一批資料畫得更平滑，不含額外資訊。
%
%   小格相位：下極 v 方向的格**邊**對齊削平面（v = 0）-> 沒有任何一格跨過削平面；
%   上極對稱於極軸。u 方向兩極皆以極軸為格心。邊界格只要**格心**落在鋼內就整格保留，
%   真實截面輪廓另以黑線疊上。
%
%   風格①粗體框圖（FS 36）；colorbar 寬度 = CBW_RATIO·figW（0.009，與 circuit_side 共用）。
%
%   MODE = 'rel'（**預設**）畫 b·â/⟨b·â⟩（各 panel 除以自己的平均、發散色階、中心 1）
%     -> 把截面內 ~2% 的結構撐滿色階。⚠ 顏色不代表物理值、也不可跨 panel 比大小
%     （P1 的 1.0 = 223 mT、P2 = 141 mT）。
%   MODE = 'abs' 畫 b·â [mT]，兩 panel 共用絕對色階 -> 呈現「下極比上極高 58%」。
%     ⚠ 該版使用者已棄（2026-08-10），保留參數但不再是預設產物。
%
%   用法：plot_pole_xsection_map                          % rel、s = 4.572 mm、CELL = 0.01 mm
%         plot_pole_xsection_map(4.572, 0.01, 'abs')      % 絕對色階版
%         plot_pole_xsection_map(4.572, 0.01, 'rel', true) % 強制重讀 .fld（2 GB x 2，約 90 s）
% =========================================================================
    clc;
    if nargin < 1 || isempty(S_MM), S_MM = 4.572; end
    if nargin < 2 || isempty(CELL), CELL = 0.01;  end   % [MODIFIED 2026-08-10] 加密 0.02→0.01
    if nargin < 3 || isempty(MODE), MODE = 'rel'; end   % [MODIFIED 2026-08-10] 預設改 rel（abs 版使用者已棄）
    if nargin < 4 || isempty(force), force = false; end

    here   = fileparts(fileparts(mfilename('fullpath')));      % → paper_fig_plot/
    figdir = fullfile(fileparts(here), 'paper_fig', 'Section3_A');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    MW = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    addpath(fullfile(MW,'function'));
    cnst = model_config('long2016_hexapole_halfcut','tip40um');

    M = cell(1,2);
    for k = 1:2
        src  = load_src(k, S_MM, here, cnst, force);           % 截面附近的源格點（快取）
        M{k} = build_map(k, S_MM, CELL, src, cnst);            % 鋼內限定內插 + 格心取值
        m = M{k};   q = ~m.extrap;                     % q = 內插格（凸包內）
        fprintf(['[P%d] s = %.3f mm｜R = %.4f mm（beta %.3f deg）｜格 %d 個（%.3f mm）｜' ...
                 'Φ = %.4f uWb｜全區平均 %.1f mT（±%.2f%%）\n' ...
                 '      內插 %d 格：%.1f ~ %.1f mT｜外插 %d 格（%.1f%%）：%.1f ~ %.1f mT\n'], ...
                 k, S_MM, m.R, m.beta_deg, numel(m.val), CELL, sum(m.val)*CELL^2*1e-3, ...
                 mean(m.val), 100*max(abs(m.val-mean(m.val)))/mean(m.val), ...
                 sum(q), min(m.val(q)), max(m.val(q)), ...
                 m.n_extrap, 100*m.n_extrap/numel(m.val), min(m.val(~q)), max(m.val(~q)));
        fprintf('      標註點（最外圈真實格點 %d 個平均）：u=%+.3f v=%+.3f -> %.2f mT\n', ...
                 m.mark.n, m.mark.u, m.mark.v, m.mark.val);
    end

    % ---- 共用座標尺度 ----
    ext  = max(cellfun(@(m) m.R, M)) + CELL/2;
    stp  = fit_step(ext/2);                            % 3 個等距 tick + 兩端留白 = 間距
    LIM  = [-2*stp, 2*stp];

    % ---- 共用色階（兩模式都共用，只是量不同）----
    if strcmpi(MODE, 'rel')
        for k = 1:2                                    % 除以自身平均（全截面，不丟格）
            mu = mean(M{k}.val);
            M{k}.val = M{k}.val/mu;   M{k}.mark.val = M{k}.mark.val/mu;
        end
        % 序列色階（白 = 最小）-> clim 下限 = 資料最小，**上限取 98 百分位**（不是最大值）。
        %   取最大值會被 P1 削平面那條 ~1.02 的窄帶把上限撐開，紅色只剩最頂端一小段、
        %   等於整條紅半段浪費掉（使用者：「紅色都被壓掉」）。收到 98% 後超出者飽和成紅，
        %   換取 0.99~1.01 這個資料最密的區間佔滿整條色階。
        allv = cell2mat(cellfun(@(m) m.val(:), M, 'UniformOutput', false).');
        lo0 = min(allv);   hi0 = prctile(allv, 98);
        %   刻度步長取 0.001 的整數倍（不強求 0.005 的「nice」值 —— 那會把上限又推回去，
        %   反而白做；規則要求的是「奇數個 + 等距」，3 位小數的標籤照樣乾淨）。
        lo  = floor(lo0/0.001)*0.001;
        stq = ceil((hi0 - lo)/4/0.001)*0.001;
        CLIM = [lo, lo + 4*stq];   CTK = lo + (0:4)*stq;
        fprintf('色階上限取 98 百分位 %.4f（資料最大 %.4f）-> 飽和格 %.1f%%\n', ...
                hi0, max(allv), 100*mean(allv > CLIM(2)));
        CLAB = '$\mathbf{b\cdot\hat{a}\;/\;\langle b\cdot\hat{a}\rangle}$';
        CMAP = cmap_wbr(256);   SFX = '_rel';
        fprintf('共用 clim = [%.3f %.3f]（資料 %.4f ~ %.4f）; 軸 = [%g %g] mm（tick 間距 %g）\n', ...
                CLIM, lo0, hi0, LIM, stp);
    else
        hi   = max(cellfun(@(m) max(m.val), M));
        CLIM = [0, ceil(hi/50)*50];   CTK = linspace(CLIM(1), CLIM(2), 5);
        CLAB = '$\mathbf{b\cdot\hat{a}\;(mT)}$';
        CMAP = turbo(256);   SFX = '';
        fprintf('共用 clim = [%g %g] mT ; 軸 = [%g %g] mm（tick 間距 %g）｜下/上 = %.1f%%\n', ...
                CLIM, LIM, stp, 100*(mean(M{1}.val)/mean(M{2}.val) - 1));
    end

    render_merged(M{2}, M{1}, S_MM, CELL, CLIM, CTK, CLAB, CMAP, SFX, LIM, stp, figdir);
end

% ============================================================================
function s = fit_step(x)
% 取「>= x 的最小 nice 步長」（figure-style：等距、疏密適中）。
% 用 >= 而非最近值 —— 步長取小了會把資料切出框外（踩過：R=0.846 配 stp=0.4）。
    n = [1 2 2.5 3 4 5 10];
    e = floor(log10(x));   c = x/10^e;
    i = find(n >= c - 1e-12, 1);   s = n(i)*10^e;
end

% ============================================================================
function [T, ah, uh, vh, R, bet] = pole_geom(k, S_MM, cnst)
% 極尖 T、極軸 ah（尖→根）、面內水平 uh、面內垂直 vh（含 +z）、站位 s 的截面半徑 R  [mm]
    T  = [cnst.pole_tip_x(k); cnst.pole_tip_y(k); cnst.pole_tip_z_wp(k)]*1e3;
    ah = cnst.pole_axis(:,k);  ah = ah/norm(ah);
    vh = [0;0;1] - dot([0;0;1],ah)*ah;   vh = vh/norm(vh);     % 面內垂直（含 +z）
    uh = cross(ah, vh);                                        % 面內水平（切向）
    low = logical(cnst.pole_is_lower(k));
    sl  = cnst.pole_cone_slope(2 - low);                       % 下極 = (1)、上極 = (2)
    bet = atan(sl);
    r_f = cnst.POLE_TIP_R*1e3;
    e_v = r_f/sin(bet) - r_f;                                  % 虛擬錐頂外移
    R   = (e_v + S_MM)*sl;
end

% ============================================================================
function src = load_src(k, S_MM, here, cnst, force)
% 截面附近的源格點（局部座標 + 鋼件遮罩），快取成小 .mat；改 CELL 不必重讀 2 GB .fld。
    cf = fullfile(here, 'data', sprintf('xsec_src_P%d_s%03d_maxwell.mat', k, round(S_MM*100)));
    if exist(cf,'file') && ~force
        src = load(cf);
        fprintf('P%d: 由快取載入 %s（鋼內 %d 點 / 框內 %d 點）\n', ...
                k, cf, sum(src.steel), numel(src.steel));
        return
    end

    [T, ah, uh, vh, R] = pole_geom(k, S_MM, cnst);
    C  = T + S_MM*ah;                                 % 截面中心（WP frame, mm）
    HB = R + 0.35;                                    % 全域框半邊（含內插餘裕）

    fprintf('P%d: 讀 %s（約 2 GB，數十秒）...\n', k, cnst.fld_files_voltage{k});
    d  = import_maxwell_fld(fullfile(cnst.fld_dir, cnst.fld_files_voltage{k}));
    X  = d.x*1e3;  Y = d.y*1e3;  Z = d.z*1e3 - cnst.SPH_OFST*1e3;   % → WP frame, mm
    in = abs(X-C(1)) <= HB & abs(Y-C(2)) <= HB & abs(Z-C(3)) <= HB;
    P  = [X(in), Y(in), Z(in)];
    B  = cnst.s_source(k) * [d.bx(in), d.by(in), d.bz(in)] * 1e3;   % all-source, mT
    clear d X Y Z in;

    % ---- 全域 → 局部座標（t 沿極軸自極尖、u 水平、v 垂直）----
    Q = P - T.';
    src.t = Q*ah;   src.u = Q*uh;   src.v = Q*vh;   src.B = B;

    % ---- 鋼件遮罩：r <= R(t)，下極再加 v <= 0 ----
    low = logical(cnst.pole_is_lower(k));
    sl  = cnst.pole_cone_slope(2 - low);   bet = atan(sl);
    r_f = cnst.POLE_TIP_R*1e3;   e_v = r_f/sin(bet) - r_f;
    Rt  = (e_v + src.t)*sl;
    src.steel = hypot(src.u, src.v) <= Rt & src.t > 0;
    if low, src.steel = src.steel & src.v <= 0; end

    bn = vecnorm(B, 2, 2);
    fprintf(['P%d: 框內 %d 點 → 鋼內 %d（‖b‖ 平均 %.1f mT）／鋼外 %d（%.1f mT）' ...
             '— 比 %.1fx，遮罩合理\n'], k, numel(src.t), sum(src.steel), ...
             mean(bn(src.steel)), sum(~src.steel), mean(bn(~src.steel)), ...
             mean(bn(src.steel))/mean(bn(~src.steel)));

    src.C = C;  src.R = R;  src.S_MM = S_MM;
    save(cf, '-struct', 'src');
    fprintf('P%d: 已存 %s\n', k, cf);
end

% ============================================================================
function m = build_map(k, S_MM, CELL, src, cnst)
% 截面內鋪小格、格心取值（**只用鋼內源格點**內插）
    [~, ~, ~, ~, R, bet] = pole_geom(k, S_MM, cnst);
    low = logical(cnst.pole_is_lower(k));

    % 小格中心（u = 水平、v = 垂直）
    n = ceil(R/CELL) + 1;
    uc = CELL*(-n:n);                          % 對稱於極軸（含 0）
    if low
        vc = -CELL*((1:n) - 0.5);              % 格邊對齊削平面 v = 0，鋼在下方
    else
        vc = CELL*(-n:n);
    end
    [UU, VV] = meshgrid(uc, vc);
    keep = (UU.^2 + VV.^2) <= R^2;             % 格心在鋼內
    UU = UU(keep);  VV = VV(keep);

    % ---- 鋼內限定內插（局部座標；linear 為主、凸包外用 nearest 補）----
    %   🔒 **先扣掉軸向線性項再內插**（使用者回報 P2 頂端突兀紅塊 → 2026-08-10 診斷）：
    %   b·â 沿極軸的梯度是 **17%/mm**（單獨解釋源雲 R²=0.97~0.98），跨一個源格距 0.1mm
    %   就變化 1.7% —— **比整個面內結構（~1.1%）還大**。3D 內插的四面體只要跨到不同 t，
    %   就把軸向梯度漏進面內圖，造成隨四面體走向而變的斑塊。做法 = 在薄板 |t-S|<=TSLAB
    %   內對 t 做線性迴歸，把每點投影回 t = S（dt=0 時修正為 0，故平面上的絕對值不變），
    %   再做 3D 內插；殘差場幾乎與 t 無關 -> 內插對四面體走向不敏感。
    %   再剝掉**錐面最外一層**源格點（離錐面 < SKIN）：那層落在鋼/空氣界面上，Maxwell
    %   匯出值被界面平均掉（實測該層 ‖b‖ min 只有 11~12 mT，往內每層 min 都 >109 mT）。
    %   🔒 **SKIN 只套錐面、不套下極的削平面**（2026-08-10 實測）：削平面**恰好與源格
    %   平面重合**（v = 0.000），格點就落在平面上、值完全乾淨（v=0 層 225.15~226.97 mT、
    %   std 0.40，連離錐面 0.06mm 的點都正常）；錐面是斜曲面才會讓格點落進混合元素。
    %   保留 v=0 層 -> P1 的資料一路涵蓋到削平面，最關鍵的位置不必外插。
    SKIN = 0.05;                                                    % [mm]
    sl2  = cnst.pole_cone_slope(2 - low);   b2 = atan(sl2);
    rf2  = cnst.POLE_TIP_R*1e3;   ev2 = rf2/sin(b2) - rf2;
    dcone = (ev2 + src.t)*sl2 - hypot(src.u, src.v);                % 離錐側面（削平面不計）
    TSLAB = 0.35;                                                   % 迴歸/內插薄板半厚 [mm]
    [~, ah] = pole_geom(k, S_MM, cnst);
    st = src.steel & dcone > SKIN & abs(src.t - S_MM) <= TSLAB;
    y  = -src.B(st, :)*ah;                      % 純量 b·â，â = -ah（朝極尖為正）[mT]
    dt = src.t(st) - S_MM;
    cf = [ones(size(dt)), dt] \ y;              % 軸向線性迴歸
    yd = y - cf(2)*dt;                          % 投影回 t = S（平面上 dt=0，絕對值不變）
    r2a = 1 - sum((y - [ones(size(dt)) dt]*cf).^2) / sum((y - mean(y)).^2);
    vfr = sum((yd - mean(yd)).^2) / sum((y - mean(y)).^2);
    fprintf('   P%d 去趨勢：軸向斜率 %+.2f mT/mm（解釋 R^2 = %.3f）→ 殘差變異剩 %.1f%%\n', ...
            k, cf(2), r2a, 100*vfr);

    % **3D** 內插（保留 t 維度）+ 凸包外用 **linear 外插**補滿到真實截面邊界。
    % ⚠ 不可改成「投影到 (u,v) 做 2D」：P1 的極軸剛好是 +x、各 t 層投影到完全相同的
    %   (u,v) 格上（重複點被平均、乾淨），但 **P2 的極軸是斜的**，投影後點位散開，2D
    %   三角化會把不同 t 的點連起來 -> 整片**水平條紋**（踩過，2026-08-10）。
    % ⚠ 外插用 linear 不用 nearest：nearest 會複製出徑向條紋，在 ±1% 尺度上像真結構。
    Pt = [src.t(st), src.u(st), src.v(st)];
    Qq = [S_MM*ones(numel(UU),1), UU, VV];
    % 內插用 'natural'（C1 連續）而非 'linear'：linear 在凸包邊界會留下階梯狀方塊補丁。
    Fe = scatteredInterpolant(Pt, yd, 'natural', 'linear');   % 內插 + 線性外插（保險）
    val = Fe(Qq);

    % ---- 外緣：自「圓形安全半徑」沿徑向線性外推補滿到真實壁面 ----
    %   最外一圈（厚 ~0.06~0.1mm）必然落在凸包外（SKIN 0.05 + 格距 0.1 的必然結果，
    %   加厚薄板也只從 14.2% 降到 12.1%）。直接吃 scatteredInterpolant 的外插會照著
    %   **鋸齒狀的凸包**起算 -> P2 頂端出現階梯方塊（踩過）。改為：取一個整圈都在凸包內
    %   的圓 r = R - RIM，沿每個方位量該處的值與徑向斜率，線性延伸到 r。起點是圓、
    %   在方位上連續 -> 邊緣平滑。⚠ r > R-RIM 的環帶是**延伸值、非原始資料**。
    RIM = 0.12;   HR = 0.05;                     % 安全內縮 / 求徑向斜率的差分步長 [mm]
    rr  = hypot(UU, VV);   th = atan2(VV, UU);
    extrap = rr > (R - RIM);
    if any(extrap)
        ce = cos(th(extrap));   se = sin(th(extrap));   ne = sum(extrap);
        f1 = Fe([S_MM*ones(ne,1), (R-RIM)   *ce, (R-RIM)   *se]);
        f2 = Fe([S_MM*ones(ne,1), (R-RIM-HR)*ce, (R-RIM-HR)*se]);
        val(extrap) = f1 + (rr(extrap) - (R-RIM)) .* ((f1-f2)/HR);
    end

    % [ADDED 2026-08-10] 標註值 = **最外圈真實格點**（不是端點的徑向延伸值）。
    %   取「近子午線（|u|<=0.06）且 |v| > 0.85R」那一小群源格點**取平均**，不是取單點：
    %   單點會被 ±0.3~0.4 mT 的點間雜訊帶偏（P2 最外那點 v=0.869 離錐壁其實比鄰點遠、
    %   值偏高 0.5%，標出來會與該處顯示的顏色矛盾）。平均後 std 只剩 0.35~0.42 mT。
    uu = src.u(st);   vv = src.v(st);
    w = abs(uu) <= 0.06 & abs(vv) > 0.85*R & (vv < 0) == low;
    if ~any(w), w = abs(uu) <= 0.06; end                    % 保險
    m.mark = struct('u', mean(uu(w)), 'v', mean(vv(w)), 'val', mean(yd(w)), 'n', sum(w));

    m.val  = val;
    m.uc   = UU;   m.vc = VV;   m.CELL = CELL;
    m.R    = R;    m.beta_deg = bet*180/pi;
    m.low  = low;  m.S_MM = S_MM;   m.pole = k;
    m.extrap = extrap;   m.n_extrap = sum(extrap);   m.SKIN = SKIN;
end

% ============================================================================
function C = cmap_wbr(n)
% 序列色階（使用者拍板 2026-08-10 定案）：**白（最小）→ 藍 → 紅（最大）**。
% 顏色值沿用 coolwarm 那三個，只是把白從中段移到最小端 —— 發散版（白=平均）整體太淡、
% 分辨度不足，故改此版換取對比。⚠ 白代表**資料最小值**，不是「平均」或「1.0」。
    a = [0.865 0.865 0.865];      % 白（coolwarm 中段色）= 最小
    b = [0.230 0.299 0.754];      % 藍（coolwarm 低端色）= 中段
    c = [0.706 0.016 0.150];      % 紅（coolwarm 高端色）= 最大
    h = floor(n/2);
    C = [a + linspace(0,1,h).'.*(b-a); b + linspace(0,1,n-h).'.*(c-b)];
end

% ============================================================================
function render_merged(mL, mR, S_MM, CELL, CLIM, CTK, CLAB, CMAP, SFX, LIM, stp, figdir)
% 合併雙 panel（左 = 上極 P2、右 = 下極 P1）+ 單一共用 colorbar
    FS = 36;                                    % 字級 36（figure-style 統一值）
    % [MODIFIED 2026-08-10] panel 縮小（820→620）：兩 panel 併排時圖寬 2175px 會讓 36pt
    %   的數字在同樣版面寬度下看起來偏小；縮 panel 讓字/圖比例貼近單 panel 的論文圖。
    H = 620; y0 = 150; leftm = 145; midgap = 105; cbgap = 22; cblab = 165; rightm = 30;
    CBW_RATIO = 0.009;                          % colorbar 寬度佔比（與 circuit_side 共用同值）

    w    = H;                                   % 兩軸等範圍 + daspect 1:1 → 正方
    base = leftm + w + midgap + w + cbgap + cblab + rightm;
    cbw  = CBW_RATIO*base/(1-CBW_RATIO);
    figW = base + cbw;   figH = y0 + H + 70;
    x1 = leftm;   x2 = x1 + w + midgap;

    fig = figure('Color','w','Units','pixels','Position',[20 40 figW figH]);
    ax1 = axes(fig,'Units','pixels');  panel(ax1, mL, CLIM, CMAP, LIM, stp, FS, true);
    ax1.Units = 'pixels';  ax1.Position = [x1 y0 w H];
    ax2 = axes(fig,'Units','pixels');  panel(ax2, mR, CLIM, CMAP, LIM, stp, FS, false);
    ax2.Units = 'pixels';  ax2.Position = [x2 y0 w H];

    cb = colorbar(ax2,'Units','pixels');
    cb.Position = [x2+w+cbgap y0 cbw H];   ax2.Position = [x2 y0 w H];
    cb.FontSize = FS;  cb.FontWeight = 'bold';
    cb.Label.Interpreter = 'latex';             % Interpreter 要先於 String
    cb.Label.String = CLAB;
    cb.Label.FontSize = FS;
    cb.Ticks = CTK;
    if all(abs(CTK-1) < 0.5)                    % rel 模式：統一 3 位小數（否則 1 會印成 "1"）
                                                %   （abs 模式的刻度是 0~250，不會誤觸）
        cb.TickLabels = compose('%.3f', CTK);
    end

    out = fullfile(figdir, sprintf('xsec_bax_P1P2_s%03d_c%03d_maxwell%s.png', ...
                   round(S_MM*100), round(CELL*1000), SFX));
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);
end

% ============================================================================
function panel(ax, m, CLIM, CMAP, LIM, stp, FS, with_ylabel)
    LWBOX = 4.0;                                % 座標框線寬（與 3D box 圖同值）
    hold(ax,'on');
    h  = m.CELL/2;
    VX = [m.uc-h, m.uc+h, m.uc+h, m.uc-h].';
    VY = [m.vc-h, m.vc-h, m.vc+h, m.vc+h].';
    F  = reshape(1:numel(VX), 4, []).';
    ec = 'none';  if m.CELL >= 0.05, ec = [0.35 0.35 0.35]; end   % 細格不畫格線
    patch(ax, 'Faces',F, 'Vertices',[VX(:) VY(:)], 'FaceVertexCData',m.val, ...
          'FaceColor','flat', 'EdgeColor',ec, 'LineWidth',0.25);
    colormap(ax, CMAP);   clim(ax, CLIM);

    % [MODIFIED 2026-08-10] 使用者拍板：拿掉磁極輪廓黑線（小格自己就描出截面形狀）

    % [MODIFIED 2026-08-10] 使用者拍板：圖上不畫端點標註（黑點 + 引線 + 數值）。
    %   m.mark（最外圈真實格點平均）仍在 build_map 算出並印到 console，供引用查核。
    box(ax,'on');  grid(ax,'off');  daspect(ax,[1 1 1]);
    % [MODIFIED 2026-08-10] TickDir 由 out 改 **in**：兩 panel 併排時，`box on` 會把左圖
    %   右框的鏡像刻度往外戳進兩圖之間的間隙，跟右圖左框的刻度一起夾住右圖的 y 標籤，
    %   看起來像「多一個 tick」。朝內就不會侵入間隙。
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',LWBOX, ...   % 框線加粗（使用者拍板 2026-08-10）
           'TickLength',[.012 .012],'TickDir','in');
    xlim(ax, LIM);  ylim(ax, LIM);
    set(ax,'XTick', (-1:1)*stp, 'YTick', (-1:1)*stp);

    yoff = LIM(1) - 0.025*diff(LIM);            % 水平軸端點只標數字、不畫 tick
    for j = 1:2
        text(ax, LIM(j), yoff, sprintf('%g',LIM(j)), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontSize',FS,'FontWeight','bold','Clipping','off');
    end

    % [MODIFIED 2026-08-10] 使用者拍板：拿掉 P1 / P2 字樣
    xlabel(ax, '$\mathbf{u\;(mm)}$', 'Interpreter','latex', 'FontSize',FS);
    if with_ylabel
        ylabel(ax, '$\mathbf{v\;(mm)}$', 'Interpreter','latex', 'FontSize',FS);
    else
        % [MODIFIED 2026-08-10] 右 panel 不印 y 刻度數字：兩 panel 的 y 尺度完全相同，
        %   而那組數字夾在兩個框之間，會被左右兩邊的框內刻度夾住、看起來像多一個 tick。
        %   共用軸的標準做法 = 只由左 panel 標 y。刻度線本身保留。
        set(ax, 'YTickLabel', []);
    end
    ax.Toolbar.Visible = 'off';  hold(ax,'off');
end
