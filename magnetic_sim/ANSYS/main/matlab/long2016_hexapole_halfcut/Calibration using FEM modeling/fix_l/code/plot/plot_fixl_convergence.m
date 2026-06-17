%% plot_fixl_convergence.m
%  ============================================================================
%  fix_l(剛性格點,無 bias)模型的「fit 收斂」圖:cost vs lsqnonlin 疊代數,每半徑一條。
%  與 plot_nofixl_convergence.m 同方法(profiled g_j、Neq 等節點、turbo、NIT=25),
%  唯一差別:**電荷鎖在理想格點 Pc_base(ê ≡ 0)、只 fit 特徵長度 ell**。
%  → 兩張圖直接對比,差異 = 加 bias 的效益(no_fix_l 的 cost floor 更低)。
%
%  說明:fix_l 把 6 顆電荷固定在 ±ell·軸,自由度只剩 ell + 各 coil 的電荷量 g_j(profile 掉)。
%  這等同於既有「自由 6x6 K̂」模型的 cost(g_j = gB·K̂ 的欄),只是用 profiled 公式表達。
%  ============================================================================
clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');
model = 'long2016_hexapole_halfcut';

cnst    = mt_constants();
R_list  = 50:50:500;                         % 同 no_fix_l
NIT     = 25;
dataset = 'all';
fig_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
data_dir = matlab_path(model, 'charge_fit', 'fitting_trend');

% ---- actuator 旋轉 + 理想格點 ----------------------------------------------
tip   = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat  = tip ./ vecnorm(tip);
R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';
Pc_base = [ 1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1];

% ---- 載 6 顆 coil 一次,旋進 actuator 框 ------------------------------------
fprintf('loading 6 coils ...\n');
N_I = 6;
d1  = import_ansys_data(ansys_path(model,'coil1','standard'), dataset, 'coil1');
air1= filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
P_meas = [d1.x(air1), d1.y(air1), d1.z(air1)-cnst.SPH_OFST];
Pa  = (R_act * P_meas.').';
r2  = sum(P_meas.^2, 2);
Ba  = zeros(size(P_meas,1), 3, N_I);
for k = 1:N_I
    if k==1, dk=d1; airk=air1;
    else
        cn=sprintf('coil%d',k);
        dk=import_ansys_data(ansys_path(model,cn,'standard'), dataset, cn);
        airk=filter_iron_nodes(dk.x,dk.y,dk.z,cnst,struct('visualize',false));
    end
    Bk = -[dk.bx(airk), dk.by(airk), dk.bz(airk)];
    Ba(:,:,k) = (R_act * Bk.').';
end

nAvail = arrayfun(@(R) nnz(r2 < (R*1e-6)^2), R_list);
Neq    = min(nAvail);
fprintf('Neq (equal nodes/shell) = %d\n', Neq);

% ---- 逐半徑跑收斂(只 fit ell,ê ≡ 0)--------------------------------------
hist = cell(1, numel(R_list));
for ri = 1:numel(R_list)
    R   = R_list(ri)*1e-6;
    idx = find(r2 < R^2);
    rng(0); idx = idx(randperm(numel(idx), Neq));
    P   = Pa(idx,:);
    Bstack = zeros(3*Neq, N_I);
    for j = 1:N_I, Bstack(:,j) = reshape(Ba(idx,:,j).', [], 1); end
    hist{ri} = fit_with_history(P, Bstack, Pc_base, NIT);
    fprintf('  R=%3d um: cost %.3e -> %.3e (%d iters)\n', R_list(ri), hist{ri}(1), hist{ri}(end), numel(hist{ri}));
end

% ---- 畫圖(與 no_fix_l 同樣式)---------------------------------------------
f = figure('Color','w','Position',[100 100 760 560]); hold on;
cmap = turbo(numel(R_list));
h = gobjects(1, numel(R_list));
for ri = 1:numel(R_list)
    h(ri) = semilogy(0:numel(hist{ri})-1, hist{ri}, '-o', 'Color', cmap(ri,:), ...
                     'LineWidth',1.4, 'MarkerSize',3, 'MarkerFaceColor', cmap(ri,:));
end
set(gca,'YScale','log'); grid on; box on; set(gca,'FontSize',11);
xlabel('iteration'); ylabel('cost');
title('fix\_l charge-model fit convergence');
legend(h, compose('R = %d \\mum', R_list), 'Location','northeast', 'NumColumns',2);
exportgraphics(f, fullfile(fig_dir,'fixl_cost_convergence.png'),'Resolution',150);
fprintf('wrote %s\n', fullfile(fig_dir,'fixl_cost_convergence.png'));

save(fullfile(data_dir,'fixl_convergence.mat'),'hist','R_list','Neq','NIT');

%% ===========================================================================
%  本地函式
%  ===========================================================================
function hist = fit_with_history(P, Bstack, Pc_base, NIT)
% FIT_WITH_HISTORY  固定疊代數的 fix_l fit(只 ell),OutputFcn 逐代記 cost。
    hist = [];
    x0 = 0.5e-3*1e3;                                           % 只有 ell(mm);ê ≡ 0
    opts = optimoptions('lsqnonlin','Display','off', ...
        'MaxFunctionEvaluations',1e7,'MaxIterations',NIT, ...
        'FunctionTolerance',0,'StepTolerance',0,'OptimalityTolerance',0, ...
        'OutputFcn',@outfun);
    lsqnonlin(@(x) resid_fixl(x, P, Bstack, Pc_base), x0, [], [], opts);
    function stop = outfun(~, ov, state)
        stop = false;
        if strcmp(state,'iter'), hist(end+1) = ov.resnorm; end %#ok<AGROW>
    end
end

function r = resid_fixl(x, P, Bstack, Pc_base)
% RESID_FIXL  profiled 殘差;電荷固定在 Pc_base(無 bias),只 ell 變。
    A = build_A(x(1)*1e-3, Pc_base, P);                        % Pc = Pc_base(ê ≡ 0)
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
