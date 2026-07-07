%% plot_singlepole_magerr_hist.m -- single-pole field VECTOR error histogram (3 shapes overlaid)
% =========================================================================
%  橫軸 = 逐節點向量誤差 e_i = |B_model - B_FEM| / |B_FEM| * 100 %  (含方向)
%  縱軸 = count。filled / halfcut / tipcut 三形狀疊圖（reference 風格：半透明無邊框 + mean 虛線）。
%  ⚠ 相減在「同一座標系」：B_model 用 dhat=[1 0 0] 在 global +x 框算；B_FEM=+BmT(SGN=+1)
%    也在同一 global +x 框（load_singlepole 不旋轉、cnst.SPH_OFST=0）→ 逐分量相減合法。
%  B_model = gI*I*(p/ell - dhat)/|p/ell - dhat|^3。取樣 R<=150µm。legend 標 mean 與 std。
%  REUSES：load_singlepole；fit 結果讀 data/singlepole_fit.mat。風格① 粗體框圖。
%  輸出 -> 本組 figures/singlepole_magerr_hist_R150.png。
% =========================================================================

clear; clc;

%% ---- config ----
R_fit  = 150e-6;      % 取樣半徑 [m]
I      = 1;           % drive current [A]
NB     = 80;          % bin 數（更細）
d_act  = [1;0;0];     % +x 極軸

%% ---- paths ----
TREE = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\Calibration_using_FEM_modeling\fix_dir';
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % mt_constants
addpath(fullfile(TREE,'code','main_function'));                                                 % load_singlepole
DATA = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil1\singlepole';
figdir = fullfile(TREE,'figures');  if ~exist(figdir,'dir'); mkdir(figdir); end

cnst = mt_constants();
cnst.SPH_OFST = 0;    % +x mesh: WP at ORIGIN (no z-shift) -> same global +x frame as B_model

%% ---- load fit results ----
S = load(fullfile(TREE,'data','singlepole_fit.mat'));   % SHAPES, ell_um, gI
SHAPES = S.SHAPES;  ns = numel(SHAPES);

%% ---- per-shape VECTOR error (same-frame subtraction) ----
E = cell(1,ns);  mu = zeros(1,ns);  sd = zeros(1,ns);  Npt = zeros(1,ns);
for s = 1:ns
    [Pw, BmT] = load_singlepole(fullfile(DATA, SHAPES{s}), cnst, 'tip');
    in = vecnorm(Pw,2,2) <= R_fit;  P = Pw(in,:);  Bfem = BmT(in,:);   % SGN=+1, global +x frame
    ell = S.ell_um(s)*1e-6;  gI = S.gI(s);
    D    = P/ell - d_act.';
    nrm  = vecnorm(D,2,2);
    Bmod = gI * I * D ./ (nrm.^3);                                     % same global +x frame
    e = vecnorm(Bmod - Bfem, 2, 2) ./ vecnorm(Bfem, 2, 2) * 100;      % vector err % (>=0)
    E{s} = e;  mu(s) = mean(e);  sd(s) = std(e);  Npt(s) = numel(e);
    fprintf('%-8s | N=%d | mean=%.3f%% | std=%.3f%% | max=%.3f%%\n', SHAPES{s}, Npt(s), mu(s), sd(s), max(e));
end

%% ---- common edges (0 .. robust upper) ----
allE = cell2mat(E(:));
BW = 0.05;  hi = prctile(allE, 99.5);  edges = 0:BW:ceil(hi/BW)*BW;   % 固定格寬 0.05%

%% ---- plot: reference 風格（半透明無邊框 + mean 虛線）+ 風格① 粗體框 ----
col = [0.35 0.45 0.75; 0.80 0.45 0.45; 0.40 0.65 0.45];   % filled=blue, halfcut=red, tipcut=green (muted)
fig = figure('Color','w','Position',[100 100 860 600]);
ax = axes(fig);  hold(ax,'on');
h = gobjects(1,ns);
for s = 1:ns
    h(s) = histogram(ax, E{s}, edges, 'FaceColor',col(s,:), 'FaceAlpha',0.45, ...
                     'EdgeColor','w', 'LineWidth',0.5);   % 半透明疊 + 白色細邊框
end
for s = 1:ns
    xl = xline(ax, mu(s), '--');  xl.Color = col(s,:)*0.75;  xl.LineWidth = 2.5;  xl.HandleVisibility = 'off';
end
hold(ax,'off');

set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
box(ax,'on'); grid(ax,'off');
xlim(ax,[0 hi]);
yl = ylim(ax); ylim(ax,[0 yl(2)*1.32]);            % headroom：留白讓 legend 不擋 bar
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
