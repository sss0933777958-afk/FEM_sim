%% verify_dform_Qform_nrmse.m  (18-param bias 版)
%  驗證「同一個 sensor 模型的兩種寫法」在 R<=50µm 的 NRMSE 完全相同：
%    form-1 (d):  b = g_H · S · diag(d) · V_j
%    form-2 (Q):  b = (k_m/ℓ̂²) · S · Q̂_j ,  Q̂ = (1/μ0)·diag(d)·V
%  因 g_H = k_m/(ℓ̂²μ0)  ⇒  兩 form 逐點相同 ⇒ NRMSE 必相等（機器精度）。
%  bias 版：空間核用 18-param 電荷格點 Pc_18、節點旋進 actuator 框（讀 calib_sensor_d_no_fix_dir.mat）。
clear; clc;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
matf = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\calib_sensor_d_no_fix_dir.mat';

cnst = mt_constants();
a2p  = [1,3,6,5,2,4];
SD = load(matf);                                       % d, gH, Vmat, ell_hat, Pc(=Pc_18), Rrot
d = SD.d; gH = SD.gH; Vmat = SD.Vmat; ell = SD.ell_hat; Pc = SD.Pc; Rrot = SD.Rrot;
km_l2 = cnst.k_m/ell^2;                                % k_m/ℓ̂²
Q = (1/cnst.mu_0)*diag(d)*Vmat;                        % 磁荷 Q̂（6x6）
exc_sign = ones(1,6); for j=1:6, if ismember(a2p(j),[1 3 6]), exc_sign(j)=-1; end, end

fprintf('loading 6 coils ...\n');
C = struct('P',{},'Bn',{});
for k=1:6
    da = import_ansys_data(fullfile(results_root, sprintf('coil%d',k), 'standard'),'wp',sprintf('coil%d',k));
    air = filter_iron_nodes(da.x,da.y,da.z,cnst,struct('visualize',false));
    C(k).P  = [da.x(air), da.y(air), da.z(air)-cnst.SPH_OFST];   % measure 框
    C(k).Bn = -[da.bx(air), da.by(air), da.bz(air)];
end

% R<=50µm 操作球(|p| 對旋轉不變,measure 框選點)
idx = find(sum(C(1).P.^2,2) <= (50e-6)^2);
PN = (Rrot * C(1).P(idx,:).').' / ell;                 % 旋進 actuator 框後正規化
Nt = numel(idx);
fprintf('R<=50um ball: N=%d nodes/coil\n\n', Nt);

evD = zeros(1,6); evQ = zeros(1,6);
for kc = 1:6
    Bf  = (Rrot * C(kc).Bn(idx,:).').';                % all-source 場(actuator 框,=-FEM)
    bf  = -exc_sign(kc)*Bf;                            % all-source target(物理)
    denom = max(vecnorm(Bf,2,2));

    wD = gH*(Vmat(:,kc).*d);                           % form-1 等效權重
    wQ = km_l2*Q(:,kc);                                % form-2 等效權重（應 == wD）
    BmD = zeros(Nt,3); BmQ = zeros(Nt,3);
    for i=1:6
        Dm=PN-Pc(:,i).'; r3=(sum(Dm.^2,2)).^1.5;      % Pc = Pc_18(actuator)
        BmD = BmD + wD(i)*(Dm./r3);
        BmQ = BmQ + wQ(i)*(Dm./r3);
    end
    evD(kc) = sqrt(mean(sum((BmD-bf).^2,2)))/denom*100;
    evQ(kc) = sqrt(mean(sum((BmQ-bf).^2,2)))/denom*100;
end

fprintf('=============== R<=50µm NRMSE：兩 form 對比 (18-param bias) ===============\n');
fprintf('  coil(激發極):    P1     P3     P6     P5     P2     P4   | worst\n');
fprintf('  form-1 (d)   : %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f | %6.3f %%\n', evD, max(evD));
fprintf('  form-2 (Q)   : %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f | %6.3f %%\n', evQ, max(evQ));
fprintf('  逐 coil max|NRMSE_d - NRMSE_Q| = %.3e %%%%\n', max(abs(evD-evQ)));
fprintf('  權重 max|wD - wQ|（全部 coil）  = %.3e\n', max(abs(km_l2*Q(:) - gH*reshape(Vmat.*d,[],1))));
fprintf('==========================================================\n');
