function [P, Bstack, info] = conv_design_ws(Nr, R, data)
% conv_design_ws -- 工作空間（球）：決定取樣點位置 + Delaunay 重心內插取場（APDL 分支）
% =========================================================================
%   [ADDED 2026-09-02 使用者指示] 與 Maxwell 分支的 conv_design_ws 同名同角色，
%   差別只在**取場引擎**：
%       Maxwell：規則格點 .fld -> 三線性內插
%       APDL   ：非結構節點雲 -> Delaunay 四面體剖分 + 重心內插（Algorithm 2）
%   產點方式（六軸殼層）與 frame / 單位慣例兩邊逐行相同，故兩支可直接對照。
%
%   介面（使用者指定）
%       [P, Bstack, info] = conv_design_ws(Nr, R, data)
%     Nr    殼層數（每根軸每個方向的等分數）；N_tot = 1 + 6*Nr
%     R     取樣球半徑 [m]
%     data  struct：
%           .cfg     model_config 輸出（必要；提供 SPH_OFST / R_act / N_I / 濾鐵幾何）
%           .x .y .z N x 1  節點座標（**APDL/measure frame**, [m]）
%           .B       N x 3 x Ns  各激發的 B（[T]，measure frame）
%           .margin  選填，篩選球外擴量 [m]，預設 100e-6
%
%   輸出
%     P       Nkeep x 3      取樣點 p_a（**actuator frame**, [m]）
%     Bstack  3*Nkeep x Ns   b(p_a)，逐點堆疊 [bx;by;bz]（**mT**, actuator frame）
%     info    struct         .N .Nr .radii .keep .npts_kept .n_iron .n_nan .B .frame
%                            .Nt .Np_tri .dV（三角化的規模與體積自檢）
%
%   ── 六步流程（對應 Maxwell 版的同名步驟）──────────────────────
%     1. 決定設計    Nr -> 六軸 x Nr 層 + 中心
%     2. 建查詢點    actuator -> measure（P_meas = p_a * R_act）-> z + SPH_OFST
%     3. 濾鐵        filter_iron_nodes(x, y, z, cfg)，吃全域座標
%     4. 三角化      節點雲篩 r <= R + margin（WP 框）-> Bowyer-Watson -> T
%                    ⚠ 只跟節點位置有關 -> **persistent 快取**，同一份資料只建一次
%     5. 定位 + 內插 每個查詢點解重心座標 lam -> b = sum_i lam_i * B(v_i)
%     6. 最後轉換    B [T] -> mT；平移 z - SPH_OFST，再旋轉 cfg.R_act -> actuator
%
%   ⚠ 這是**引擎**（無對應圖），勿當孤兒刪除。
%   ⚠ 三角化成本高（34k 點約 163 s），但只跟節點位置有關；掃階梯時 Nr 改變**不會**
%     重建，只重跑第 5 步。
% =========================================================================
    MARGIN_DEF = 100e-6;                       % 使用者拍板 2026-09-02

    assert(isstruct(data) && isfield(data,'cfg'), 'conv_design_ws:data', ...
           'data 必須是含 .cfg/.x/.y/.z/.B 的 struct');
    cfg    = data.cfg;
    margin = getdef_(data, 'margin', MARGIN_DEF);
    Ns     = size(data.B, 3);

    %% ==== Step 1+2：六軸殼層產點 -> measure frame ========================
    [Pm, info] = axes_shells_(Nr, R, cfg.R_act);
    x = Pm(:,1);   y = Pm(:,2);   z = Pm(:,3) + cfg.SPH_OFST;   % -> APDL 框
    info.N = numel(x);

    %% ==== Step 3：濾鐵（吃全域座標，零轉換）==============================
    keep = filter_iron_nodes(x, y, z, cfg);   keep = keep(:);
    info.n_iron = nnz(~keep);

    %% ==== Step 4：三角化（persistent 快取）===============================
    R_keep = R + margin;
    [T, Pt, Ainv, tinfo] = tri_cache_(data, cfg, R_keep);
    info.Nt = tinfo.Nt;   info.Np_tri = tinfo.Np;   info.dV = tinfo.dV;

    %% ==== Step 5：定位 + 重心內插 =======================================
    zw = z - cfg.SPH_OFST;                                  % 查詢點回 WP 框
    pq = [x(keep), y(keep), zw(keep)];
    [tetID, lam] = locate_bf_(pq, Ainv);
    Bsub = data.B(tinfo.idx, :, :);                         % **與 Pt 同一套索引**
    Bq   = interp_bary_(T, tetID, lam, Bsub, tinfo.Np);     % Nq x 3 x Ns [T]

    nan_q = isnan(tetID);
    info.n_nan = nnz(nan_q);
    if info.n_nan > 0
        warning('conv_design_ws:outside', ...
                ['%d 個取樣點落在三角化凸包外 -> B 為 NaN。' ...
                 'margin (%.0f um) 可能不足。'], info.n_nan, margin*1e6);
    end

    %% ==== Step 6：單位與 frame ==========================================
    B = 1e3 * Bq;                                           % T -> mT
    Pr = (cfg.R_act * [x(keep), y(keep), zw(keep)].').';    % -> actuator frame
    for m = 1:Ns
        B(:,:,m) = (cfg.R_act * B(:,:,m).').';              % B 只旋轉
    end

    nkeep  = size(Pr,1);
    P      = Pr;
    Bstack = zeros(3*nkeep, Ns);                            % 逐點 [bx;by;bz]
    for m = 1:Ns
        Bstack(:,m) = reshape(B(:,:,m).', [], 1);
    end
    info.keep = keep;   info.npts_kept = nkeep;   info.B = B;   info.frame = 'actuator';
    report_(info);
end

% ============================================================================
function [Pm, info] = axes_shells_(Nr, R, RACT)
% 六軸殼層產點 —— 與 Maxwell 分支的 axes_shells_ 逐行相同（Algorithm 1）。
%      1: d     <- R / N_a
%      2: N_tot <- 1 + 6*N_a
%      5:        中心點
%   6~13: 軸 j -> 正負 s -> 等分 i
%   回傳轉到 measure frame 的 Pm；info.P_act 保留 actuator frame。
    assert(~isempty(R) && isscalar(R) && R > 0, 'conv_design_ws:noR', '必須給球半徑 R');
    assert(~isempty(Nr) && Nr >= 1, 'conv_design_ws:noNr', '必須給殼層數 Nr');
    Na   = Nr;
    d    = R / Na;
    Ntot = 1 + 6*Na;
    pa   = zeros(Ntot, 3);
    E    = eye(3);
    cnt  = 1;   pa(cnt,:) = [0, 0, 0];
    for j = 1:3
        for s = [+1, -1]
            for i = 1:Na
                cnt = cnt + 1;
                pa(cnt,:) = s * i * d * E(j,:);
            end
        end
    end
    % actuator -> measure：R_act 的**列**是 measure 座標下的致動軸單位向量
    Pm   = pa * RACT;
    info = struct('Nr',Nr, 'R',R, 'radii',(1:Na)*d, 'P_act',pa, 'npts_design',Ntot);
end

% ============================================================================
function [T, P, Ainv, tinfo] = tri_cache_(data, cfg, R_keep)
% 三角化只跟**節點位置**有關 -> 同一份節點雲 + 同一個 R_keep 只建一次。
%   掃階梯（Nr 變）時直接命中快取；換模型/換 variant/換 R_keep 才重建。
    persistent C
    key = sprintf('%d|%.9g|%.9g|%.9g|%.9g', numel(data.x), R_keep, ...
                  data.x(1), data.y(end), data.z(round(end/2)));
    if ~isempty(C) && strcmp(C.key, key)
        T = C.T;   P = C.P;   Ainv = C.Ainv;   tinfo = C.tinfo;   return
    end
    zw  = data.z - cfg.SPH_OFST;
    idx = find(sqrt(data.x.^2 + data.y.^2 + zw.^2) <= R_keep);
    P   = [data.x(idx), data.y(idx), zw(idx)];
    fprintf('  [tri] 建 Delaunay：%d 個節點（r <= %.0f um）...\n', numel(idx), R_keep*1e6);
    t0 = tic;   T = bowyer_watson_(P);   t_tri = toc(t0);

    [~, Vhull] = convhull(P(:,1), P(:,2), P(:,3));
    Vtri = sum(tet_volume_(P, T));   dV = abs(Vtri - Vhull)/Vhull;
    if dV > 1e-8
        warning('conv_design_ws:volume', '剖分體積與凸包差 %.2e -> 剖分可能不完整', dV);
    end
    Ainv  = tet_ainv_(P, T);
    tinfo = struct('idx',idx, 'Np',numel(idx), 'Nt',size(T,1), 'dV',dV, ...
                   'R_keep',R_keep, 't_tri',t_tri);
    fprintf(['  [tri] Nt = %d（Nt/Np = %.2f）｜%.1f s｜體積自檢 %.1e' newline], ...
            tinfo.Nt, tinfo.Nt/tinfo.Np, t_tri, dV);
    C = struct('key',key, 'T',T, 'P',P, 'Ainv',Ainv, 'tinfo',tinfo);
end

% ============================================================================
function T = bowyer_watson_(P)
% 教科書 Bowyer-Watson。超級四面體用使用者定義：
%     S_1..4 = c + M*L*{(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)}
%   c = WP 框原點、L = 點雲半徑、M = 1000。
%   ⚠ M 太小會讓外殼上正確的四面體「吃到」超級頂點而被誤挖，最後刪除超級四面體時
%     那塊就空了（M=20 實測缺 16 顆、體積短少 2.5e-5）。M=1000 起與內建逐顆相同。
    M   = 1000;
    n   = size(P,1);
    L   = max(vecnorm(P, 2, 2));
    sv  = M * L * [ 1  1  1;  1 -1 -1; -1  1 -1; -1 -1  1];
    Pa  = [P; sv];

    cap  = max(16*n, 64);
    T    = zeros(cap,4);   cc = zeros(cap,3);   r2 = zeros(cap,1);   live = false(cap,1);
    T(1,:) = [n+1 n+2 n+3 n+4];   live(1) = true;   nt = 1;
    [cc(1,:), r2(1)] = circum_rows_(Pa, T(1,:));
    RELTOL = 1e-10;                                  % 相對容差（絕對值會誤挖）

    for i = 1:n
        p   = Pa(i,:);
        act = find(live(1:nt));
        dd  = (cc(act,1)-p(1)).^2 + (cc(act,2)-p(2)).^2 + (cc(act,3)-p(3)).^2;
        bad = act(dd <= r2(act) .* (1 + RELTOL));
        if isempty(bad), continue; end

        Tb = T(bad,:);
        F  = [Tb(:,[1 2 3]); Tb(:,[1 2 4]); Tb(:,[1 3 4]); Tb(:,[2 3 4])];
        [u, ~, ic] = unique(sort(F,2), 'rows');
        Bf = u(accumarray(ic,1) == 1, :);            % 只出現一次 = 洞穴邊界

        live(bad) = false;
        nb = size(Bf,1);
        if nt + nb > cap
            k = find(live(1:nt));   m = numel(k);
            T(1:m,:) = T(k,:);   cc(1:m,:) = cc(k,:);   r2(1:m) = r2(k);
            live(:) = false;   live(1:m) = true;   nt = m;
            if nt + nb > cap
                cap = max(2*cap, nt + nb + 1024);
                T(cap,4)=0; cc(cap,3)=0; r2(cap,1)=0; live(cap,1)=false;
            end
        end
        rows = nt + (1:nb);
        T(rows,:) = [Bf, repmat(i, nb, 1)];
        [cc(rows,:), r2(rows)] = circum_rows_(Pa, T(rows,:));
        live(rows) = true;   nt = nt + nb;
    end
    T = T(live(1:nt), :);
    T(any(T > n, 2), :) = [];                        % 移除含超級頂點者
end

% ============================================================================
function Ainv = tet_ainv_(P, T, VMIN)
% 每顆四面體的重心座標反矩陣：A = [x1..x4; y1..y4; z1..z4; 1 1 1 1]，lam = A \ [p;1]。
%   體積過小的退化 sliver 填 NaN -> 定位時永不被選中（NaN >= -tol 為 false）。
    if nargin < 3 || isempty(VMIN), VMIN = 1e-6 * median(tet_volume_(P, T)); end
    nt = size(T,1);
    A  = ones(4,4,nt);
    for i = 1:4, A(1:3,i,:) = permute(P(T(:,i),:), [2 3 1]); end
    Ainv = pageinv(A);
    Ainv(:,:,tet_volume_(P,T) < VMIN) = NaN;
end

% ============================================================================
function [tetID, lam] = locate_bf_(pq, Ainv, TOL)
% 暴力定位：對每個查詢點掃全部四面體，取 lam 四個都 >= -TOL 者。
%   多顆同時符合（點落在面/邊上）-> 取 min(lam) 最大者，結果有確定性。
%   找不到 = 查詢點在凸包外 -> NaN（margin 不足的警報）。
    if nargin < 3 || isempty(TOL), TOL = 1e-12; end
    nq = size(pq,1);
    tetID = nan(nq,1);   lam = nan(nq,4);
    for q = 1:nq
        Lq = reshape(pagemtimes(Ainv, [pq(q,:).'; 1]), 4, []);
        [best, j] = max(min(Lq, [], 1));
        if best >= -TOL, tetID(q) = j;   lam(q,:) = Lq(:,j).'; end
    end
end

% ============================================================================
function Bq = interp_bary_(T, tetID, lam, Bsub, Np_expect)
% Algorithm 2 第 5~7 行：Bq(q,:,s) = sum_i lam(q,i) * Bsub(T(tetID(q),i), :, s)
%   ⚠ Bsub 必須已用 tinfo.idx 切成與三角化的 P 同一套索引。
    Np = size(Bsub,1);   Ns = size(Bsub,3);   nq = numel(tetID);
    assert(Np == Np_expect, ['Bsub 有 %d 列，但三角化的 P 有 %d 列 ' ...
           '-> 沒有切成同一套索引'], Np, Np_expect);
    Bq = nan(nq, 3, Ns);
    ok = ~isnan(tetID);
    if ~any(ok), return; end
    V  = T(tetID(ok), :);   nk = size(V,1);
    Bv = reshape(Bsub(V(:),:,:), nk, 4, 3, Ns);
    Bq(ok,:,:) = reshape(sum(lam(ok,:) .* Bv, 2), nk, 3, Ns);
end

% ============================================================================
function v = tet_volume_(P, T)
    ba = P(T(:,2),:) - P(T(:,1),:);
    ca = P(T(:,3),:) - P(T(:,1),:);
    da = P(T(:,4),:) - P(T(:,1),:);
    v  = abs(sum(ba .* cross(ca, da, 2), 2)) / 6;
end

% ============================================================================
function [cc, r2] = circum_rows_(P, T)
% 四面體外接球心與半徑平方（對 T 每列向量化，Cramer 展開）
    a = P(T(:,1),:);  b = P(T(:,2),:);  c = P(T(:,3),:);  d = P(T(:,4),:);
    ba = b - a;   ca = c - a;   da = d - a;
    n1 = sum(ba.^2,2);   n2 = sum(ca.^2,2);   n3 = sum(da.^2,2);
    x1 = cross(ca, da, 2);   x2 = cross(da, ba, 2);   x3 = cross(ba, ca, 2);
    det = sum(ba .* x1, 2);
    cc  = a + 0.5 * (n1.*x1 + n2.*x2 + n3.*x3) ./ det;
    r2  = sum((cc - a).^2, 2);
end

% ============================================================================
function report_(info)
    fprintf('  [axsh] R=%.0f um：Nr=%d -> 六軸 x %d 層 + 中心，產 %d 點；r = %.1f ... %.1f um\n', ...
            info.R*1e6, info.Nr, info.Nr, info.N, info.radii(1)*1e6, info.radii(end)*1e6);
    fprintf('  [sample] %d 點 -> 保留 %d（鐵件 %d、凸包外 %d）；frame=%s\n', ...
            info.N, info.npts_kept, info.n_iron, info.n_nan, info.frame);
end

% ============================================================================
function v = getdef_(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
