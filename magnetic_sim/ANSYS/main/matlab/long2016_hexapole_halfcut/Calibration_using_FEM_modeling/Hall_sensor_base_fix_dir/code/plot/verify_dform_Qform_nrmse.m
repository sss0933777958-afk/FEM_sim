%% verify_dform_Qform_nrmse.m
%  驗證「同一個 sensor 模型的兩種寫法」在 R<=50µm 的 NRMSE 完全相同：
%    form-1 (d):  b = g_H · S · diag(d) · V_j
%    form-2 (Q):  b = (k_m/ℓ̂²) · S · Q̂_j ,  Q̂ = (1/μ0)·diag(d)·V
%  因 g_H = k_m/(ℓ̂²μ0)  ⇒  兩 form 逐點相同 ⇒ NRMSE 必相等（機器精度）。
%  NRMSE 定義同 sweep_alln_vs_R.m：每 coil sqrt(mean‖Δ‖²)/max‖B_FEM‖×100，worst-coil。
%  讀 calib_sensor_d.mat（由 ../main/main.m 產出）。
clear; clc;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
ddir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit';

cnst = mt_constants();
a2p  = [1,3,6,5,2,4];
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

SD = load(fullfile(ddir,'calib_sensor_d.mat'));        % d, gH, Vmat, ell_hat
d = SD.d; gH = SD.gH; Vmat = SD.Vmat; ell = SD.ell_hat;
km_l2 = cnst.k_m/ell^2;                                % k_m/ℓ̂²
Q = (1/cnst.mu_0)*diag(d)*Vmat;                        % 磁荷 Q̂（6x6）
exc_sign = ones(1,6); for j=1:6, if ismember(a2p(j),[1 3 6]), exc_sign(j)=-1; end, end

fprintf('loading 6 coils ...\n');
C = struct('P',{},'Bn',{});
for k=1:6
    da = import_ansys_data(fullfile(results_root, sprintf('coil%d',k), 'standard'),'wp',sprintf('coil%d',k));
    air = filter_iron_nodes(da.x,da.y,da.z,cnst,struct('visualize',false));
    C(k).P  = [da.x(air), da.y(air), da.z(air)-cnst.SPH_OFST];
    C(k).Bn = -[da.bx(air), da.by(air), da.bz(air)];
end

% R<=50µm 操作球
idx = find(sum(C(1).P.^2,2) <= (50e-6)^2);
Ptest = C(1).P(idx,:); Nt = numel(idx); PN = Ptest/ell;
fprintf('R<=50um ball: N=%d nodes/coil\n\n', Nt);

evD = zeros(1,6); evQ = zeros(1,6);
for kc = 1:6
    Bf = C(kc).Bn(idx,:); denom = max(vecnorm(Bf,2,2));
    bf = -exc_sign(kc)*Bf;                             % all-source target

    wD = gH*(Vmat(:,kc).*d);                           % form-1 等效權重
    wQ = km_l2*Q(:,kc);                                % form-2 等效權重（應 == wD）
    BmD = zeros(Nt,3); BmQ = zeros(Nt,3);
    for i=1:6
        Dm=PN-dhat(:,i).'; r3=(sum(Dm.^2,2)).^1.5;
        BmD = BmD + wD(i)*(Dm./r3);
        BmQ = BmQ + wQ(i)*(Dm./r3);
    end
    evD(kc) = sqrt(mean(sum((BmD-bf).^2,2)))/denom*100;
    evQ(kc) = sqrt(mean(sum((BmQ-bf).^2,2)))/denom*100;
end

fprintf('=============== R<=50µm NRMSE：兩 form 對比 ===============\n');
fprintf('  coil(激發極):    P1     P3     P6     P5     P2     P4   | worst\n');
fprintf('  form-1 (d)   : %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f | %6.3f %%\n', evD, max(evD));
fprintf('  form-2 (Q)   : %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f | %6.3f %%\n', evQ, max(evQ));
fprintf('  逐 coil max|NRMSE_d - NRMSE_Q| = %.3e %%%%\n', max(abs(evD-evQ)));
fprintf('  權重 max|wD - wQ|（全部 coil）  = %.3e\n', max(abs(km_l2*Q(:) - gH*reshape(Vmat.*d,[],1))));
fprintf('==========================================================\n');
