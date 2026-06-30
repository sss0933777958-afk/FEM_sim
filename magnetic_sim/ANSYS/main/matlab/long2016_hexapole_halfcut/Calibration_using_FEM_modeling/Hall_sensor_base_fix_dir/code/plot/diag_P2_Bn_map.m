function diag_P2_Bn_map(PREVIEW)
% DIAG_P2_BN_MAP  P2-under-P1 sensor 讀值 sign 診斷
%   問題：all-source 圖在 P2 sensor 上方那塊「往上/往外」磁通，造成模型 P2-under-P1 讀正；
%         實驗是負（磁通匯進 P2）。本診斷沿 P2 cone 兩側 flank、不同沿面距尖 SOFF、不同
%         air-gap，掃 all-source mean(B·n+) 的正負，找出「匯進極(負)」區在哪、現有 sensor
%         （flank A, SOFF=4.572mm, AIR=0.41mm）落在正區或負區。
%
%   資料源：baseline coil1/standard（P1 激發），真實 FEM 節點（球內平均，非內插）。
%   all-source：P1 為下極 → exc_sign=-1，故 val = -mean(B_raw·n+)。
%
%   PREVIEW=true 輸出暫存 png；false 存 figures/P2_Bn_sign_map_standard.png。

    if nargin < 1 || isempty(PREVIEW), PREVIEW = false; end

    %% ---- paths ----
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants, import_ansys_data
    cnst = mt_constants();
    rr   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

    %% ---- 載 baseline coil1 (P1 激發) ----
    d = import_ansys_data(fullfile(rr,'coil1','standard'),'all','coil1');
    X = [d.x, d.y, d.z - cnst.SPH_OFST];                 % 節點座標(WP 框) [m]
    B = [d.bx, d.by, d.bz];                              % raw FEM B [T]
    fprintf('載入 coil1/standard：%d 節點，|B|max=%.4f T\n', size(X,1), max(vecnorm(B,2,2)));

    exc_sign = -1;                                       % all-source：P1 下極翻號

    %% ---- P2 cone 幾何 (i=2，上極) ----
    th     = cnst.pole_angles(2)*pi/180;                 % P2 方位角
    inc_up = cnst.upper_incline;                         % 上極軸傾角 ≈36.59° [rad]
    beta   = atan2(cnst.POLE_R, cnst.POLE_CONE_LEN);     % 半錐角 ≈11.31° [rad]
    dirv   = @(e,a) [cos(e)*cos(a); cos(e)*sin(a); sin(e)];
    tip    = [cnst.R_norm_xy*cos(th); cnst.R_norm_xy*sin(th); +cnst.R_norm_z];

    % 兩側 flank：A=上側(現用，外法線往上/外)、B=下側(WP-facing，外法線往下/外)
    flank(1).name='A 上側(現用)'; flank(1).slant=dirv(inc_up+beta,th);      flank(1).n=dirv(inc_up+beta+pi/2,th);
    flank(2).name='B 下側(WP面)'; flank(2).slant=dirv(inc_up-beta,th);      flank(2).n=dirv(inc_up-beta-pi/2,th);

    %% ---- 掃描參數 ----
    SOFF = (0.5:0.5:12)*1e-3;        % 沿面距尖 [m]
    AIRs = [0.2 0.41 0.8]*1e-3;      % air-gap [m]
    RS   = 0.3e-3;                   % 球取樣半徑 [m]（粗網格穩定 sign 用，比 sensor 盤略大）
    KNN  = 5;                        % 球內無點時兜底：取最近 KNN 個節點
    cur  = struct('soff',4.572e-3,'air',0.41e-3,'flank',1);  % 現有 sensor

    % 區域取點：球內節點；空則取最近 KNN（回傳 B·n+ 平均 + 命中數）
    sample = @(pos,nn) local_mean_Bn(X, B, pos, nn, RS, KNN);

    %% ---- 計算 ----
    val = nan(2,numel(AIRs),numel(SOFF));   % signed all-source mean(B·n+)
    cnt = zeros(2,numel(AIRs),numel(SOFF)); % 取樣節點數
    for f = 1:2
        sl = flank(f).slant; nn = flank(f).n;
        for ia = 1:numel(AIRs)
            for is = 1:numel(SOFF)
                pos = tip + SOFF(is)*sl + AIRs(ia)*nn;        % 候選 sensor 中心
                [mb, nc] = sample(pos, nn);
                val(f,ia,is) = exc_sign * mb;  cnt(f,ia,is) = nc;
            end
        end
    end

    % 現有 sensor 點的值（同法）
    pc = tip + cur.soff*flank(cur.flank).slant + cur.air*flank(cur.flank).n;
    [mbc, ncc] = sample(pc, flank(cur.flank).n);
    valc = exc_sign * mbc;
    fprintf('\n現有 sensor (flank A, SOFF=4.572mm, AIR=0.41mm)：all-source mean(B·n+)=%+.3e T (%s)，取樣 %d 節點\n', ...
            valc, ternary(valc<0,'負=匯進P2','正=回流'), ncc);

    %% ---- signed-log 變換（保留正負、壓縮量級）----
    slog = @(v) sign(v).*log10(max(abs(v),1e-9)/1e-9);  % 0→0；正負分明

    %% ---- 圖 ----
    fig = figure('Position',[60 60 1100 820],'Color','w');
    cols = lines(numel(AIRs));
    for f = 1:2
        ax = subplot(2,1,f); hold(ax,'on'); grid(ax,'on');
        for ia = 1:numel(AIRs)
            y = slog(squeeze(val(f,ia,:)));
            plot(SOFF*1e3, y, '-o','Color',cols(ia,:),'LineWidth',1.6,'MarkerSize',4, ...
                 'DisplayName',sprintf('gap=%.2fmm',AIRs(ia)*1e3));
        end
        yline(ax,0,'k-','LineWidth',1.2,'HandleVisibility','off');     % sign 翻轉線
        if f==cur.flank
            plot(cur.soff*1e3, slog(valc),'p','MarkerSize',16,'MarkerFaceColor',[1 .84 0], ...
                 'MarkerEdgeColor','k','LineWidth',1,'DisplayName','現有 sensor');
        end
        xlabel('沿錐面距極尖 SOFF [mm]'); ylabel('signed log_{10}|B·n_+|  (>0 回流, <0 匯進P2)');
        title(sprintf('P2 cone flank %s ｜ all-source mean(B·n+) sign (P1 激發, baseline)', flank(f).name),'Interpreter','none');
        legend(ax,'Location','best'); xlim([0 12]);
    end

    %% ---- 輸出 ----
    if PREVIEW
        out = fullfile(tempdir,'P2_Bn_sign_map_preview.png');
    else
        out = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
               'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\figures\' ...
               'P2_Bn_sign_map_standard.png'];
    end
    exportgraphics(fig, out, 'Resolution', 200);
    fprintf('saved: %s\n', out);

    %% ---- 文字摘要：每 flank 在 gap=0.41mm 的 sign 翻轉位置 ----
    ia0 = find(abs(AIRs-0.41e-3)<1e-9,1);
    for f = 1:2
        y = squeeze(val(f,ia0,:));
        neg = SOFF(y<0)*1e3; pos = SOFF(y>0)*1e3;
        fprintf('flank %s (gap0.41): 負(匯進)區 SOFF=%s mm ; 正(回流)區 SOFF=%s mm\n', ...
                flank(f).name, mat2str(round(neg,1)), mat2str(round(pos,1)));
    end
end

function o = ternary(c,a,b), if c, o=a; else, o=b; end, end

function [mb, nc] = local_mean_Bn(X, B, pos, nn, RS, KNN)
% 球內(半徑 RS)真實節點的 mean(B·n+)；球內無點 → 取最近 KNN 個。
    dr  = vecnorm(X - pos.', 2, 2);
    sel = dr <= RS;
    if ~any(sel)
        [~, idx] = mink(dr, KNN); sel = idx;   % 兜底：最近 KNN
    end
    Bn = B(sel,:) * nn;  mb = mean(Bn);  nc = numel(Bn);
end
