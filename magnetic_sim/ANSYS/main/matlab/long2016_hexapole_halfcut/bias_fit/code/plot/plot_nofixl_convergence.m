%% plot_nofixl_convergence.m
%  ============================================================================
%  no_fix_l(18-param bias)模型的「fit 收斂」圖:cost(殘差平方和)vs lsqnonlin 疊代數,
%  每個取樣半徑一條曲線。比照既有 plot_KI_convergence_gB50.m / KI_cost_convergence_gB50.png。
%
%  收斂擷取:lsqnonlin + OutputFcn,每代記 optimValues.resnorm(= sum r^2 = cost);
%  容差全設 0 + MaxIterations=NIT → 硬跑滿 NIT 代(收斂後 cost 持平,尾段變水平線)。
%  等節點數比較:各殼抽相同 Neq 點(= 最小殼節點數),隔離模型 misfit、非節點數差異。
%  註:bias 模型 g_j 已 profile 掉、無 gB0(故無 reference 的 gB50 multistart),收斂只解 [ell; e]。
%  ============================================================================
clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
model = 'long2016_hexapole_halfcut';

cnst      = mt_constants();
R_list    = 50:50:500;                       % 10 個取樣半徑殼 [um]
NIT       = 25;                              % 固定疊代數(硬跑滿,看收斂)
apdl_to_paper_idx = [1,3,6,5,2,4];
dataset   = 'all';
fig_dir   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
data_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\bias_fit\data';

% ---- actuator 旋轉 + 理想格點 ----------------------------------------------
tip   = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat  = tip ./ vecnorm(tip);
R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';                    % û=P1, v̂=P3, ŵ=P5
Pc_base = [ 1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1];

% ---- 載 6 顆 coil 一次,旋進 actuator 框 ------------------------------------
fprintf('loading 6 coils ...\n');
N_I = 6;
d1  = import_ansys_data(ansys_path(model,'data','coil1','standard'), dataset, 'coil1');
air1= filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
P_meas = [d1.x(air1), d1.y(air1), d1.z(air1)-cnst.SPH_OFST];
Pa  = (R_act * P_meas.').';
r2  = sum(P_meas.^2, 2);
Ba  = zeros(size(P_meas,1), 3, N_I);
for k = 1:N_I
    if k==1, dk=d1; airk=air1;
    else
        cn=sprintf('coil%d',k);
        dk=import_ansys_data(ansys_path(model,'data',cn,'standard'), dataset, cn);
        airk=filter_iron_nodes(dk.x,dk.y,dk.z,cnst,struct('visualize',false));
    end
    Bk = -[dk.bx(airk), dk.by(airk), dk.bz(airk)];             % all-source
    Ba(:,:,k) = (R_act * Bk.').';
end

% ---- 等節點數:Neq = 各殼 in-ball 數的最小值 -------------------------------
nAvail = arrayfun(@(R) nnz(r2 < (R*1e-6)^2), R_list);
Neq    = min(nAvail);
fprintf('Neq (equal nodes/shell) = %d\n', Neq);

% ---- 逐半徑跑收斂 ----------------------------------------------------------
hist = cell(1, numel(R_list));
for ri = 1:numel(R_list)
    R   = R_list(ri)*1e-6;
    idx = find(r2 < R^2);
    rng(0); idx = idx(randperm(numel(idx), Neq));             % 抽 Neq 個(可重現)
    P   = Pa(idx,:);
    Bstack = zeros(3*Neq, N_I);
    for j = 1:N_I, Bstack(:,j) = reshape(Ba(idx,:,j).', [], 1); end
    hist{ri} = fit_with_history(P, Bstack, Pc_base, NIT);
    fprintf('  R=%3d um: cost %.3e -> %.3e (%d iters)\n', R_list(ri), hist{ri}(1), hist{ri}(end), numel(hist{ri}));
end

% ---- 畫圖 ------------------------------------------------------------------
f = figure('Color','w','Position',[100 100 760 560]); hold on;
cmap = turbo(numel(R_list));
h = gobjects(1, numel(R_list));
for ri = 1:numel(R_list)
    h(ri) = semilogy(0:numel(hist{ri})-1, hist{ri}, '-o', 'Color', cmap(ri,:), ...
                     'LineWidth',1.4, 'MarkerSize',3, 'MarkerFaceColor', cmap(ri,:));
end
set(gca,'YScale','log'); grid on; box on; set(gca,'FontSize',11);
xlabel('iteration'); ylabel('cost');
title('no\_fix\_l charge-model fit convergence');
legend(h, compose('R = %d \\mum', R_list), 'Location','northeast', 'NumColumns',2);
exportgraphics(f, fullfile(fig_dir,'nofixl_cost_convergence.png'),'Resolution',150);
fprintf('wrote %s\n', fullfile(fig_dir,'nofixl_cost_convergence.png'));

save(fullfile(data_dir,'nofixl_convergence.mat'),'hist','R_list','Neq','NIT');

%% ===========================================================================
%  本地函式
%  ===========================================================================
function hist = fit_with_history(P, Bstack, Pc_base, NIT)
% FIT_WITH_HISTORY  固定疊代數的 bias fit,用 OutputFcn 逐代記 cost(殘差平方和)。
    hist = [];
    x0 = [0.5e-3*1e3; zeros(17,1)];                            % 初值:ell=0.5mm、e=0
    opts = optimoptions('lsqnonlin','Display','off', ...
        'MaxFunctionEvaluations',1e7,'MaxIterations',NIT, ...
        'FunctionTolerance',0,'StepTolerance',0,'OptimalityTolerance',0, ...  % 不提前停
        'OutputFcn',@outfun);
    lsqnonlin(@(x) bias_resid(x, P, Bstack, Pc_base), x0, [], [], opts);
    function stop = outfun(~, ov, state)
        stop = false;
        if strcmp(state,'iter'), hist(end+1) = ov.resnorm; end %#ok<AGROW>
    end
end

function r = bias_resid(x, P, Bstack, Pc_base)
% BIAS_RESID  profiled 殘差(g_j 用最小二乘解掉)。
    A = build_A(x(1)*1e-3, make_Pc(x(2:18), Pc_base), P);
    M = A.' * A;
    r = [];
    for j = 1:size(Bstack,2)
        gj = M \ (A.' * Bstack(:,j));
        r  = [r; A*gj - Bstack(:,j)]; %#ok<AGROW>
    end
end

function A = build_A(ell, Pc, P)
    Np = size(P,1); pbar = P/ell; A = zeros(3*Np,6);
    for k = 1:6
        D = pbar - Pc(:,k).'; r3 = sum(D.^2,2).^1.5; Sk = D ./ r3;
        A(:,k) = reshape(Sk.', 3*Np, 1);
    end
end

function Pc = make_Pc(e17, Pc_base)
    E = zeros(3,6);
    E(:,1)=e17(1:3); E(:,2)=e17(4:6); E(:,3)=e17(7:9);
    E(:,4)=e17(10:12); E(:,5)=e17(13:15); E(1,6)=e17(16); E(2,6)=e17(17);
    E(3,6)=e17(1)-e17(4)+e17(8)-e17(11)+e17(15);
    Pc = Pc_base + E;
end
