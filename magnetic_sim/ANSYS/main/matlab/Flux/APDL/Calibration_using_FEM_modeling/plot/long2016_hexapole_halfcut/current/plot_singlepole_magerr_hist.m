%% plot_singlepole_magerr_hist.m -- single-pole BIAS field VECTOR error histogram (3 shapes overlaid)
% =========================================================================
%  no_fix_l 版（bias）：對照 current_base/code/plot/plot_singlepole_magerr_hist.m（無-bias）。
%  橫軸 = 逐節點向量誤差 e_i = |B_model - B_FEM| / |B_FEM| * 100 %。
%  縱軸 = count。filled / halfcut / tipcut 三形狀疊圖（reference 風格：半透明 + mean 虛線）。
%  ⚠ 相減在「同一座標系」：B_model 用 pc_hat=[1,e_y,e_z] 在 global +x 框算；B_FEM=+BmT(SGN=+1)
%    也在同一 global +x 框（load_singlepole 不旋轉、cnst.SPH_OFST=0）→ 逐分量相減合法。
%  B_model = gI*I*(p/ell - pc_hat)/|p/ell - pc_hat|^3。取樣 R<=150µm。
%  REUSES：current_base 的 load_singlepole；fit 結果讀 data/singlepole_bias_fit.mat。風格① 粗體框圖。
%  對照無-bias：bias 讓 halfcut 分布明顯左移（誤差下降）。
%  輸出 -> 本組 figures/singlepole_magerr_hist_R150.png。
% =========================================================================

clear; clc;

%% ---- config ----
R_fit  = 150e-6;      % 取樣半徑 [m]
I      = 1;           % drive current [A]

%% ---- paths ----
TREE     = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\Calibration_using_FEM_modeling\current_base';
FIX_MFUN = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\Calibration_using_FEM_modeling\current_base\code\main_function';  % load_singlepole
% [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live config。
CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
addpath(FIX_MFUN);                                                                              % load_singlepole (reuse)
DATA = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil1\singlepole';
figdir = fullfile(TREE,'figures','eighteen_param');  if ~exist(figdir,'dir'); mkdir(figdir); end

cnst = model_config('long2016_hexapole_halfcut','tip40um');
cnst.SPH_OFST = 0;    % +x mesh: WP at ORIGIN (no z-shift) -> same global +x frame as B_model

%% ---- load bias fit results ----
S = load(fullfile(TREE,'data','singlepole_bias_fit.mat'));   % SHAPES, ell_um, ey, ez, gI
SHAPES = S.SHAPES;  ns = numel(SHAPES);

%% ---- per-shape VECTOR error (bias model, same-frame subtraction) ----
E = cell(1,ns);  mu = zeros(1,ns);  sd = zeros(1,ns);  Npt = zeros(1,ns);
for s = 1:ns
    [Pw, BmT] = load_singlepole(fullfile(DATA, SHAPES{s}), cnst, 'tip');
    in = vecnorm(Pw,2,2) <= R_fit;  P = Pw(in,:);  Bfem = BmT(in,:);   % SGN=+1, global +x frame
    ell = S.ell_um(s)*1e-6;  pchat = [1, S.ey(s), S.ez(s)];           % bias charge位置
    D    = P/ell - pchat;
    nrm  = vecnorm(D,2,2);
    Bmod = S.gI(s) * I * D ./ (nrm.^3);                               % same global +x frame
    e = vecnorm(Bmod - Bfem, 2, 2) ./ vecnorm(Bfem, 2, 2) * 100;      % vector err % (>=0)
    E{s} = e;  mu(s) = mean(e);  sd(s) = std(e);  Npt(s) = numel(e);
    fprintf('%-8s | N=%d | mean=%.3f%% | std=%.3f%% | max=%.3f%%\n', SHAPES{s}, Npt(s), mu(s), sd(s), max(e));
end

%% ---- common edges（figure-style 規則：nb=180 均分資料合併全距）----
allE = cell2mat(E(:));
nb = 180;  edges = linspace(min(allE), max(allE), nb+1);   % 共用 edges、離散 180 bins

%% ---- plot: reference 風格（半透明無邊框 + mean 虛線）+ 風格① 粗體框 ----
col = [0.35 0.45 0.75; 0.80 0.45 0.45; 0.40 0.65 0.45];   % filled=blue, halfcut=red, tipcut=green (muted)
fig = figure('Color','w','Position',[100 100 860 600]);
ax = axes(fig);  hold(ax,'on');
h = gobjects(1,ns);
for s = 1:ns
    h(s) = histogram(ax, E{s}, edges, 'FaceColor',col(s,:), 'FaceAlpha',0.55, ...
                     'EdgeColor','w', 'LineWidth',0.3);   % 規則：FaceAlpha 0.55 + 白細邊
end
for s = 1:ns
    xl = xline(ax, mu(s), '--');  xl.Color = col(s,:)*0.75;  xl.LineWidth = 2.5;  xl.HandleVisibility = 'off';
end
hold(ax,'off');

set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box(ax,'on'); grid(ax,'off');
xlim(ax,[min(allE) max(allE)]);                    % 資料全距
yl = ylim(ax); ylim(ax,[0 yl(2)*1.20]);            % headroom：留白讓 legend 不擋 bar
xt = get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
yt = get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
xlabel(ax,'|B_{model} - B_{FEM}| / |B_{FEM}| (%)','FontWeight','bold');
ylabel(ax,'Count','FontWeight','bold');

lg = arrayfun(@(s) sprintf('%s  (mean=%.2f%%, std=%.2f%%)', SHAPES{s}, mu(s), sd(s)), ...
              1:ns, 'UniformOutput',false);
legend(ax, h, lg, 'Location','northeast', 'FontSize',13, 'FontWeight','bold', 'Box','on');

%% ---- output ----
outpng = fullfile(figdir, 'singlepole_magerr_hist_R150.png');
exportgraphics(fig, outpng, 'Resolution',1200);
fprintf('已輸出 %s\n', outpng);
