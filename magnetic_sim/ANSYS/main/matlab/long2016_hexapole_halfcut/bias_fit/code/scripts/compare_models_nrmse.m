%% compare_models_nrmse.m
%  用同一把 canonical NRMSE 尺，比較兩個場模型：
%    (A) K̂_I 自由擬合（36 DOF，R=150 那組 ell/gB/K̂）
%    (B) Hall-sensor d-model（6 DOF，d/Vmat，no-gain：b=S·V·d，ell_hat=0.856）
%  NRMSE 定義照搬 sweep_alln_vs_R.m:60-66：
%    每 coil  NRMSE_kc = sqrt(mean_i ‖B_mod−B_FEM‖²)/max_i‖B_FEM‖×100
%    整體     = worst-coil = max over 6
%  區域：50µm 操作球（canonical）+ R≤150µm（fit 區，對照）。
%  自洽檢查：K̂_I 的 50µm worst-coil NRMSE 須重現 sweep 的 NRMSE_R(R=150)≈1.23%。
clear; clc;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
ddir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\bias_fit\data';   % 規則#2：sweep_alln 已移至 bias_fit/data

cnst = mt_constants();
a2p  = [1,3,6,5,2,4];
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

% ---- K̂_I 模型參數（R=150）----
SW = load(fullfile(ddir,'sweep_alln_vs_R.mat'));   % R_um, ell_R, gB_R, Ksave, NRMSE_R
gi = find(SW.R_um==150,1);
ellK = SW.ell_R(gi); gBK = SW.gB_R(gi); Khat = SW.Ksave(:,:,gi);
fprintf('K_I (R=150): ell=%.4f mm, gB=%.4e, sweep NRMSE_R=%.3f%% (自洽目標)\n', ellK*1e3, gBK, SW.NRMSE_R(gi));

% ---- sensor 模型參數（no-gain：模型 b = S·V·d，無 g_H）----
SD = load(['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\' ...
           'Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir\data\calib_sensor_d.mat']);  % 規則#2：Hall_sensor_base_fix_dir/data/
d = SD.d; Vmat = SD.Vmat; ellS = SD.ell_hat;
exc_sign = ones(1,6); for j=1:6, if ismember(a2p(j),[1 3 6]), exc_sign(j)=-1; end, end
fprintf('sensor   : ell_hat=%.4f mm (no-gain d)\n\n', ellS*1e3);

% ---- 載 6 coils 'wp' ----
fprintf('loading 6 coils ...\n');
C = struct('P',{},'Bn',{});
for k=1:6
    da = import_ansys_data(fullfile(results_root, sprintf('coil%d',k), 'standard'),'wp',sprintf('coil%d',k));
    air = filter_iron_nodes(da.x,da.y,da.z,cnst,struct('visualize',false));
    C(k).P  = [da.x(air), da.y(air), da.z(air)-cnst.SPH_OFST];
    C(k).Bn = -[da.bx(air), da.by(air), da.bz(air)];     % negate 慣例（= K̂_I target）
end

regions = struct('name',{'50um ball','R<=150um'},'R',{50e-6,150e-6});
res = struct();
for rr = 1:2
    Rcut = regions(rr).R;
    idx  = find(sum(C(1).P.^2,2) <= Rcut^2);
    Ptest = C(1).P(idx,:); Nt = numel(idx);

    evK = zeros(1,6); evS = zeros(1,6);
    PNk = Ptest/ellK; PNs = Ptest/ellS;
    for kc = 1:6
        Bf = C(kc).Bn(idx,:);                 % -B_FEM at test nodes
        denom = max(vecnorm(Bf,2,2));

        % (A) K̂_I model
        Iv = zeros(6,1); Iv(a2p(kc))=1; w = gBK*Khat*Iv; Bm = zeros(Nt,3);
        for i=1:6, Dm=PNk-dhat(:,i).'; r3=(sum(Dm.^2,2)).^1.5; Bm=Bm+w(i)*(Dm./r3); end
        evK(kc) = sqrt(mean(sum((Bm-Bf).^2,2)))/denom*100;

        % (B) sensor d-model:  bm = S*diag(Vmat(:,kc))*d ;  bf = -exc_sign*Bn （no-gain，無 g_H）
        wS = Vmat(:,kc).*d;                   % 6x1 等效權重 = V_kc ⊙ d
        Bm2 = zeros(Nt,3);
        for i=1:6, Dm=PNs-dhat(:,i).'; r3=(sum(Dm.^2,2)).^1.5; Bm2=Bm2+wS(i)*(Dm./r3); end
        bf2 = -exc_sign(kc)*Bf;               % calib_fem 的 all-source target（norm 與 Bf 同）
        evS(kc) = sqrt(mean(sum((Bm2-bf2).^2,2)))/denom*100;
    end
    res(rr).Nt=Nt; res(rr).evK=evK; res(rr).evS=evS;

    fprintf('\n=============== NRMSE @ %s  (N=%d nodes/coil) ===============\n', regions(rr).name, Nt);
    fprintf('  coil(P-excite):    P1     P3     P6     P5     P2     P4   | worst\n');
    fprintf('  K_I  (36 DOF) : %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f | %6.3f %%\n', evK, max(evK));
    fprintf('  sensor (6 DOF): %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f | %6.3f %%\n', evS, max(evS));
end

fprintf('\n=============== 摘要（worst-coil NRMSE）===============\n');
fprintf('  region        K_I(36DOF)   sensor(6DOF)\n');
fprintf('  50um ball     %7.3f %%    %7.3f %%\n', max(res(1).evK), max(res(1).evS));
fprintf('  R<=150um      %7.3f %%    %7.3f %%\n', max(res(2).evK), max(res(2).evS));
