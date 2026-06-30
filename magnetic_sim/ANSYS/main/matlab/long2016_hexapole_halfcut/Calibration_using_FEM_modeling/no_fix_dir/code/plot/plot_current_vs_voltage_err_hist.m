%% plot_current_vs_voltage_err_hist.m -- 電荷(電流)模型 vs 電壓(sensor)模型 場誤差疊圖
% =========================================================================
% Current_base（電荷模型）：B = S·(gB·K̄·F)  → 電荷 Gc = gB·K̄·F（no_fix_dir gauge_KI）。
% Voltage_base（電壓模型）：B = S·(g_V·D̄·V) → 電荷 Gv = g_V·D̄·Vmat（Hall_sensor_base_no_fix_dir）。
% 兩者數學上恆等於 Dv=(AᵀA)⁻¹AᵀBstack（gB·K̄·F=Dv、g_V·D̄=D 且 D·Vmat=Dv）→ 場/誤差/變異數相同
% （電壓重建 ≡ 電流重建）。各自獨立算、疊圖呈現重合。
% 誤差 = 逐點逐激發 |B_model−B_FEM|（向量差大小，mT）；bin 0.005；選項①粗體框圖。
% =========================================================================
clear; clc;
OUTDIR = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
          'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\no_fix_dir\figures'];   % <== 本副本輸出位置

VARIANT='gap200um_mueq'; R_select=150e-6; ell0=0.5e-3;
CAL = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
       'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
addpath(fullfile(CAL,'no_fix_dir','code','function'));     % load_coils_actuator/select_ball/fit_bias/make_Pc/build_A/gauge_KI
model='long2016_hexapole_halfcut'; cnst=mt_constants(); apdl_to_paper_idx=[1,3,6,5,2,4];
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

% ---- no_fix bias 擬合（一次）→ A, Bstack, Dv ----
D = load_coils_actuator(model, cnst, apdl_to_paper_idx, 'all', VARIANT);
[P, Bstack, npts] = select_ball(D, R_select);
s_sink=ones(1,6); for j=1:6, if ismember(apdl_to_paper_idx(j),[1 3 6]), s_sink(j)=-1; end; end
Bstack = (-Bstack).*s_sink;     % flip-sink all-source（與 Hall no_fix main_Dmatrix 同慣例，否則上極號不一致）
[ell, e_hat] = fit_bias(P, Bstack, D.Pc_base, ell0);
Pc = make_Pc(e_hat, D.Pc_base);
A  = build_A(ell, Pc, P);
Dv = (A.'*A)\(A.'*Bstack);

% ---- Current_base：Gc = gB·K̄·F ----
[Kbar, gB] = gauge_KI(ell, Pc, P, Bstack, D.F);
Gc = gB * Kbar * D.F;
err_c = vecnorm(reshape(A*Gc - Bstack, 3, npts*6), 2, 1).';

% ---- Voltage_base：Gv = g_V·D̄·Vmat ----
H = load(fullfile(CAL,'Hall_sensor_base_no_fix_dir','data',sprintf('calib_D_%s.mat',VARIANT)));  % Dmat,D_bar,Vmat_p
g_V  = (6/5)*H.Dmat(1,1);
Vmat = H.Vmat_p(:, apdl_to_paper_idx);                  % sensor paper × 激發 APDL
Gv = g_V * H.D_bar * Vmat;
err_v = vecnorm(reshape(A*Gv - Bstack, 3, npts*6), 2, 1).';

fprintf('自檢：max|Gc−Gv|=%.2e ；max|Gc−Dv|=%.2e ；max|Gv−Dv|=%.2e\n', ...
        max(abs(Gc(:)-Gv(:))), max(abs(Gc(:)-Dv(:))), max(abs(Gv(:)-Dv(:))));
ec=err_c*1e3; ev=err_v*1e3;                              % mT
fprintf('Current_base: var=%.4g mT^2, median=%.4f mT ；Voltage_base: var=%.4g mT^2, median=%.4f mT\n', ...
        var(ec), median(ec), var(ev), median(ev));

%% ---- 疊圖（選項①粗體框圖、半透明、bin 0.005）----
XMAX=0.8; edges=0:0.005:XMAX;
fig=figure('Color','w','Position',[100 100 820 580]); ax=axes(fig); hold(ax,'on');
histogram(ax, ec, edges, 'FaceColor',[0.85 0.33 0.10], 'FaceAlpha',0.55, 'EdgeColor','w');
histogram(ax, ev, edges, 'FaceColor',[0.20 0.40 0.70], 'FaceAlpha',0.55, 'EdgeColor','w');
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]); box(ax,'on'); grid(ax,'off');
xlim(ax,[0 XMAX]); xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end)); yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
xlabel(ax,'|B_{model} - B_{FEM}| (mT)','FontWeight','bold'); ylabel(ax,'Count','FontWeight','bold');
lc=sprintf('Current\\_base  (\\sigma^2 = %.4g mT^2)', var(ec));
lv=sprintf('Voltage\\_base  (\\sigma^2 = %.4g mT^2)', var(ev));
lg=legend(ax,{lc,lv},'Location','northeast','Interpreter','tex'); set(lg,'FontSize',14,'FontWeight','bold','Box','on');

out=fullfile(OUTDIR,'current_vs_voltage_err_hist.png');
exportgraphics(fig,out,'Resolution',600); fprintf('已輸出 %s\n', out);
