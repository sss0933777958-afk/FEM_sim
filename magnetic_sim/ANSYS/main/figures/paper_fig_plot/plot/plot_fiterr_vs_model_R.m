function plot_fiterr_vs_model_R(force, USE_BIAS, SRC)
%PLOT_FITERR_VS_MODEL_R  「用大範圍擬合的模型去描述 R≤150µm 資料」的誤差曲線。
%   plot_fiterr_vs_model_R()                    % APDL graded、有快取就用
%   plot_fiterr_vs_model_R(true)                % 強制重算（每個 R 重跑一次 18 參數擬合，很慢）
%   plot_fiterr_vs_model_R(true, true,'maxwell')% 改用 Maxwell WP 細格（0.02mm）
%   SRC = 'apdl'（預設，graded .dat）| 'maxwell'（B_p<k>.fld，±0.6mm 框涵蓋 R≤500µm）
%
%   橫軸 R_model = 300:20:500 µm：用 r ≤ R_model 的資料擬合 (l_hat, e)，並在同一顆球上解 G。
%   縱軸 = 把**整顆模型（l_hat, e, G 全部固定）**拿去描述 r ≤ 150 µm 資料的 RMSPE [%]
%          RMSPE = sqrt(Σ‖S₁₅₀·G − b₁₅₀‖² / Σ‖b₁₅₀‖²)·100
%   另存一條對照：幾何(l_hat,e)沿用、但 G 在 150 球上重解（只測幾何、不測振幅）。
%
%   ⚠ R_model = 150 的自身值（in-sample）另存於 rec.rmspe_self，供參考不畫。
%   風格①粗體框圖。輸出 figures/paper_fig/Section2_E/fiterr_vs_model_R_apdl_<model>.png

    if nargin < 1 || isempty(force),    force = false; end
    if nargin < 2 || isempty(USE_BIAS), USE_BIAS = true; end
    if nargin < 3 || isempty(SRC),      SRC = 'apdl'; end
    HERE   = fileparts(fileparts(mfilename('fullpath')));
    FIGDIR = fullfile(fileparts(HERE), 'paper_fig', 'Section2_E');
    if ~exist(FIGDIR,'dir'), mkdir(FIGDIR); end
    mstr  = 'single'; if USE_BIAS, mstr = 'eighteen'; end
    CACHE = fullfile(HERE, 'data', sprintf('fiterr_vs_model_R_%s_%s.mat', lower(SRC), mstr));
    FS = 36;  R_IN = 150e-6;  l0 = 0.5e-3;
    Rlist = (300:20:500)*1e-6;

    if ~force && exist(CACHE,'file')
        S = load(CACHE);  fprintf('cache hit: %s\n', CACHE);
    else
        solver_path(SRC);                                       % [MODIFIED] 掛對應分支、移除另一分支
        cfg = model_config('long2016_hexapole_halfcut','tip40um');
        if strcmpi(SRC,'maxwell')
            raw = extract_maxwell_data(cfg,'all',cfg.default_variant);
        else
            raw = extract_ansys_data(cfg,'all',cfg.default_variant);
        end
        ad  = build_actuator_data(raw,cfg);
        F = zeros(6,cfg.N_I);  for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j),j) = 1; end
        [P0,B0,n0] = cfg.select_ball(ad, R_IN);                 % 評估用的內圈資料
        fprintf('內圈 R<=%.0fum: %d 點\n', R_IN*1e6, n0);

        nR = numel(Rlist);
        S = struct('Rum',Rlist*1e6, 'rmspe',zeros(1,nR), 'rmspe_refitG',zeros(1,nR), ...
                   'ell',zeros(1,nR), 'gI',zeros(1,nR), 'npts',zeros(1,nR), 'USE_BIAS',USE_BIAS);
        for i = 1:nR
            [P,B,np] = cfg.select_ball(ad, Rlist(i));
            [e,l_hat] = fitting(P, B, cfg.Pc_base, l0, USE_BIAS);          % 在大球上擬合幾何
            [~, gI, G] = solve_current(l_hat, e, cfg.Pc_base, P, B, F);    % 同一顆球解 G
            % ① 整顆模型固定（l_hat, e, G）→ 在內圈算殘差
            S150 = kern(P0, l_hat, mkPc(e, cfg.Pc_base));
            S.rmspe(i) = sqrt(sum((S150*G - B0).^2,'all')/sum(B0.^2,'all'))*100;
            % ② 幾何沿用、G 在內圈重解（只測幾何）
            [~,~,~,rm2] = solve_current(l_hat, e, cfg.Pc_base, P0, B0, F);
            S.rmspe_refitG(i) = rm2.RMSPE;
            S.ell(i) = l_hat*1e6;  S.gI(i) = gI;  S.npts(i) = np;
            fprintf('R_model=%3.0fum  n=%7d  l=%6.1fum  gI=%7.4f | RMSPE(150) 固定G %6.3f%%  重解G %6.3f%%\n', ...
                    Rlist(i)*1e6, np, l_hat*1e6, gI, S.rmspe(i), S.rmspe_refitG(i));
        end
        [e0,l0h] = fitting(P0,B0,cfg.Pc_base,l0,USE_BIAS);                  % in-sample 參考值
        [~,~,~,rm0] = solve_current(l0h,e0,cfg.Pc_base,P0,B0,F);
        S.rmspe_self = rm0.RMSPE;  S.ell_self = l0h*1e6;
        fprintf('in-sample (R_model=150): l=%.1fum  RMSPE=%.3f%%\n', S.ell_self, S.rmspe_self);
        save(CACHE,'-struct','S');  fprintf('saved %s\n',CACHE);
    end

    %% ---- 畫圖（風格①）----
    % [MODIFIED 2026-08-06] 縱軸標題太長會被裁掉 → 拆兩行 + 加高畫布（字級維持 36）
    fig = figure('Color','w','Position',[80 60 1250 1000]);
    ax  = axes(fig);  hold(ax,'on');
    plot(ax, S.Rum, S.rmspe, '-o', 'Color',[0.00 0.45 0.74], 'LineWidth',3.5, ...
         'MarkerSize',10, 'MarkerFaceColor',[0.00 0.45 0.74]);
    hold(ax,'off');
    xt = 300:50:500;                                    % 5 個等距 tick、含起訖點
    [yl, yt] = ticks_y(min(S.rmspe), max(S.rmspe));
    xlim(ax,[300 500]);  ylim(ax,yl);
    set(ax,'XTick',xt,'YTick',yt);
    set(ax,'FontSize',FS,'FontWeight','bold','LineWidth',3.0,'TickLength',[.018 .018]);
    box(ax,'on');  grid(ax,'off');
    xlabel(ax,'$\mathbf{Model\;sampling\;range\;(micro\;meter)}$','Interpreter','latex','FontSize',FS);
    ylabel(ax,{'$\mathbf{Fitting\;error\;within}$', ...
               '$\mathbf{R\le150\;micro\;meter\;(\%)}$'},'Interpreter','latex','FontSize',FS);
    ax.Toolbar.Visible = 'off';

    out = fullfile(FIGDIR, sprintf('fiterr_vs_model_R_%s_%s.png', lower(SRC), mstr));
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('wrote %s\n', out);  close(fig);
    fprintf('R_model 300→500: RMSPE(150) %.3f%% → %.3f%%（in-sample %.3f%%）\n', ...
            S.rmspe(1), S.rmspe(end), S.rmspe_self);
end

% ---- 依 SRC 掛上對應分支、並移除另一分支（兩包有同名函式，只 addpath 會被遮蔽）----
function CAL = solver_path(SRC)
    APDL = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    MW   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\Maxwell';
    switch lower(SRC)
        case 'apdl',    CAL = APDL;  OTHER = MW;
        case 'maxwell', CAL = MW;    OTHER = APDL;
        otherwise, error('SRC 必為 ''apdl'' | ''maxwell''');
    end
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(fullfile(OTHER,'function'));  rmpath(fullfile(OTHER,'common_path'));
    warning('on','MATLAB:rmpath:DirNotFound');
    addpath(fullfile(CAL,'function'));   addpath(fullfile(CAL,'common_path'));
    fprintf('[SRC=%s] model_config -> %s\n', SRC, which('model_config'));   % 自證用哪一份 code
end

% ---- 電荷格 Pc = Pc_base + E(e)（含 e6z 約束；與 solve_current 一致）----
function Pc = mkPc(e17, Pc_base)
    if isempty(e17) || all(e17(:)==0), Pc = Pc_base; return; end
    E = zeros(3,6);
    E(:,1)=e17(1:3);  E(:,2)=e17(4:6);  E(:,3)=e17(7:9);
    E(:,4)=e17(10:12); E(:,5)=e17(13:15);
    E(1,6)=e17(16);   E(2,6)=e17(17);
    E(3,6)=e17(1)-e17(4)+e17(8)-e17(11)+e17(15);
    Pc = Pc_base + E;
end

% ---- 無因次 Coulomb kernel（3Np×6；與 solve_current 的 build_S 同式）----
function S = kern(P, l_hat, Pc)
    Np = size(P,1);  pbar = P/l_hat;  S = zeros(3*Np,6);
    for k = 1:6
        d  = pbar - Pc(:,k).';
        r3 = sum(d.^2,2).^1.5;
        S(:,k) = reshape((d./r3).', 3*Np, 1);
    end
end

% ---- 縱軸：奇數個等距 tick、首末不標、兩端留白 = 間距 ----
function [lim, tk] = ticks_y(ymin, ymax)
    nice = [1 2 2.5 3 4 5 10];
    for e = -4:3
        for q = nice
            s = q*10^e;
            lo = floor(ymin/s)*s;
            if lo + 4*s >= ymax*1.0001
                tk = lo + [s 2*s 3*s];  lim = [lo, lo+4*s];  return
            end
        end
    end
    s = (ymax-ymin)/3;  tk = ymin + [s 2*s 3*s];  lim = [ymin-s, ymax+s];
end
