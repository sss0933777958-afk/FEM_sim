%% plot_Dmodel_field_err_hist.m -- NTU: 18-param bias 場誤差直方圖
% =========================================================================
%  NTU full_assembly：18-param bias 電荷模型重建場 |B_model − B_FEM| 逐點逐激發向量差。
%  自足：載 NTU 6-coil → select_ball(WP,300µm) → all-source → bias fit → A·Dv 重建 → 殘差。
%  風格：選項① 粗體框、y=Percentage(%)、legend 標 mean + CV、histogram bars（nb=180）。
%  輸出 → voltage_base/figures/eighteen_param/Dmodel_field_err_hist.png。
% =========================================================================
clear; clc;

VARIANT = 'full_assembly';  R_select = 300e-6;  ell0 = 0.5e-3;
CAL = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
       'NTU_hexapole\Calibration_using_FEM_modeling'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');   % import_ansys_data
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');  % ansys_path
addpath(fullfile(CAL,'current_base','code','main_function'));   % NTU mt_constants/load/select_ball/fitting/build_S_matrix (LAST→NTU 版勝出)
model = 'NTU_hexapole';  apdl_to_paper_idx = [1,3,6,5,2,4];
figdir = fullfile(CAL,'voltage_base','figures','eighteen_param');  if ~exist(figdir,'dir'); mkdir(figdir); end

cnst = mt_constants();
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, 'all', VARIANT);
[P, Bstack, npts] = select_ball(D, R_select);
s_sink = ones(1,6);  for j=1:6, if ismember(apdl_to_paper_idx(j),[1 3 6]), s_sink(j)=-1; end; end
Bstack = (-Bstack).*s_sink;                                 % all-source（已 mT）

% ---- 18-param bias 場重建殘差 ----
[ellb,~,Pcb,~] = fitting(P, Bstack, D.Pc_base, ell0, true);
Ab = build_S_matrix(ellb, Pcb, P);  Dvb = (Ab.'*Ab)\(Ab.'*Bstack);
errb = vecnorm(reshape(Ab*Dvb - Bstack, 3, npts*6), 2, 1).';   % mT

mb=mean(errb); cvb=std(errb)/mb;
fprintf('18-param bias: mean=%.4f mT CV=%.3f (N=%d)\n', mb,cvb,numel(errb));

% ---- 直方圖（Percentage、選項① 粗體框、nb=180）-----------------------------
XMAX = prctile(errb, 99.5);  nb = 180;  edges = linspace(0, XMAX, nb+1);
fig = figure('Color','w','Position',[100 100 900 620]);  ax = axes(fig);  hold(ax,'on');
histogram(ax, errb, edges, 'Normalization','probability', 'FaceColor',[0.35 0.50 0.80], 'FaceAlpha',0.55, 'EdgeColor','w','LineWidth',0.3);
xline(ax, mb, '--', 'Color',[0.20 0.30 0.60], 'LineWidth',2);           % mean 虛線
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);  box(ax,'on'); grid(ax,'off');
xlim(ax,[0 XMAX]);
yl=ylim(ax); ylim(ax,[0 yl(2)*1.20]);                                   % headroom 給 legend
xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));
yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end)); set(ax,'YTickLabel', arrayfun(@(v)sprintf('%g',v*100), get(ax,'YTick'),'uni',0));   % 機率→%
xlabel(ax,'|B_{model} - B_{FEM}| (mT)','FontWeight','bold');  ylabel(ax,'Percentage (%)','FontWeight','bold');
lg = legend(ax, {sprintf('18-param bias  (mean %.3f mT, CV %.2f)', mb, cvb)}, ...
            'Location','northeast','Interpreter','tex');  set(lg,'FontSize',14,'FontWeight','bold','Box','on');

outpng = fullfile(figdir,'Dmodel_field_err_hist.png');
exportgraphics(fig, outpng, 'Resolution',300);
fprintf('已輸出 %s\n', outpng);
