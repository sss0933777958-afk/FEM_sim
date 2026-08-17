function [tri, npts] = conv_design(R, cfg, F, USE_BIAS, here, l0, sg)
% conv_design -- 取某取樣半徑 R 的「收斂點設計」三元組 (N_r, N_phi, N_theta)
% =========================================================================
%   ⚠ 這是**引擎**（無對應圖），被 plot_err_hist_shell.m / plot_err_hist_overlay.m
%     呼叫，勿當孤兒刪除。
%
%   判準 = 三者交集，各自持續 KWIN = 10 個設計：
%     ① l_hat  相對前一個設計的變化 < 0.5%
%     ② g_I    同一把尺
%     ③ K_I_bar 符合物理：off-diag 全負 + 對角全正 + 對角占優
%   （與 plot_full_vs_conv_vs_R.m 完全一致；理由見該檔 conv_fit 的註解。）
%
%   先查 data/full_vs_conv_vs_R_maxwell.mat（R-sweep 是 80:20:500）；
%   該 R 不在網格上（例如 150）才就地沿階梯搜一次。
%
%   R        : 取樣半徑 [m]
%   USE_BIAS : false = single、true = eighteen
%   here     : paper_fig_plot/ 的絕對路徑（用來找 data/ 快取）
% =========================================================================
    if nargin < 6 || isempty(l0), l0 = 0.5e-3; end
    % [ADDED 2026-08-15] sg：傳給 sphere_grid_sample 的 model 選項（.model/.geom）。
    %   空 = long2016（既有行為）。另可用 .ki_gate=false 放寬判準（見下）。
    if nargin < 7 || isempty(sg), sg = struct(); end
    SG = struct('frame','actuator');
    if isfield(sg,'model'), SG.model = sg.model; end
    if isfield(sg,'geom'),  SG.geom  = sg.geom;  end
    KI_GATE = true;   if isfield(sg,'ki_gate'), KI_GATE = logical(sg.ki_gate); end
    % [ADDED 2026-08-15] ki_req=false → **K̄_I 完全不參與判準**，只用 l_hat + g_I 收斂
    %   （使用者拍板，用於 zhi_peng：六極不等強，K̄_I 的物理結構條件本來就不適用）。
    %   ki_gate 只放寬「非對角全負」；ki_req 才是「這條判準要不要算」。
    KI_REQ  = true;   if isfield(sg,'ki_req'),  KI_REQ  = logical(sg.ki_req);  end
    is_long = ~isfield(sg,'model') || strcmp(sg.model,'long2016_hexapole_halfcut');
    TOL = 0.005;   KWIN = 10;   NDMAX = 150;

    % ---- 先查 R-sweep 快取（**只對 long2016 有效**：那顆快取是它專屬的）----
    df = fullfile(here, 'data', 'full_vs_conv_vs_R_maxwell.mat');
    if is_long && exist(df, 'file')
        D  = load(df);
        ia = find(abs(D.R_um - R*1e6) < 1e-6, 1);
        if ~isempty(ia)
            fld = 'tri_c1';  if USE_BIAS, fld = 'tri_c2'; end
            tri = D.(fld)(ia,:);
            if all(tri > 0)
                npts = prod(tri);
                fprintf('  [conv_design] R=%.0f um（快取）設計 (%d,%d,%d)、%d 點\n', ...
                        R*1e6, tri, npts);
                return
            end
        end
    end

    % ---- 快取沒有 → 就地搜尋 ----
    w = [1 3 3*pi];   t = [1 2 3];   TRI = zeros(NDMAX,3);
    for q = 1:NDMAX, TRI(q,:) = t;  [~,j] = min(t./w);  t(j) = t(j)+1;  end
    [ell, gv] = deal(nan(1,NDMAX));   ok = false(1,NDMAX);
    for q = 1:NDMAX
        o = SG;   o.NRPT = TRI(q,:);
        [x,y,z,B] = sphere_grid_sample(R, [], o);
        P = [x y z];   np = size(P,1);
        Bs = zeros(3*np, size(B,3));
        for j = 1:size(B,3), Bs(:,j) = reshape(B(:,:,j).', [], 1); end
        try
            [e, l] = fitting(P, Bs, cfg.Pc_base, l0, USE_BIAS);
            [K, gv(q)] = solve_current(l, e, cfg.Pc_base, P, Bs, F);
            ell(q) = l*1e6;
            od = K(~eye(6));   [~,am] = max(abs(K),[],2);
            % KI_GATE=false 時放寬「非對角全負」（六極不等強的設計，如 zhi_peng，
            %   該條恆不成立 → 判準永遠達不到）；對角全正 + 對角占優仍要求。
            ok(q) = all(diag(K) > 0) && isequal(am(:).', 1:6) && (~KI_GATE || all(od < 0));
        catch
        end
        ie = fstab(ell(1:q), TOL, KWIN);
        ig = fstab(gv(1:q),  TOL, KWIN);
        if KI_REQ, ik = ftrue(ok(1:q), KWIN); else, ik = 1; end   % ki_req=false → 只看 l_hat + g_I
        if ~isnan(ie) && ~isnan(ig) && ~isnan(ik)
            tri  = TRI(max([ie ig ik]), :);
            npts = prod(tri);
            fprintf('  [conv_design] R=%.0f um（搜尋）設計 (%d,%d,%d)、%d 點\n', ...
                    R*1e6, tri, npts);
            return
        end
    end
    error('conv_design:noconv', 'R=%.0f um 在 %d 個設計內未達判準', R*1e6, NDMAX);
end

% ============================================================================
function i0 = fstab(v, tol, K)
% 第一個 i，使其後連續 K 個「相對前一點」的變化率都 < tol。
    i0 = NaN;   if numel(v) < K+1, return; end
    rel = abs(diff(v) ./ v(1:end-1));
    for i = 1:numel(rel)-K+1
        if all(rel(i:i+K-1) < tol), i0 = i;  return; end
    end
end

% ============================================================================
function i0 = ftrue(v, K)
% 第一個 i，使 v(i..i+K-1) 全為 true。
    i0 = NaN;
    for i = 1:numel(v)-K+1
        if all(v(i:i+K-1)), i0 = i;  return; end
    end
end
