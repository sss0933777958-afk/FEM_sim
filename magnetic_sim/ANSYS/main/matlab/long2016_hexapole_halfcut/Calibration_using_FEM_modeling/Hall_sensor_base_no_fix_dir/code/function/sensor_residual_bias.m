function nrmse_sensor = sensor_residual_bias(A, Bstack, Vmat, exc_sign, d, gH, N_I)
% SENSOR_RESIDUAL_BIAS  sensor 模型殘差(actuator 框,build_A 堆疊版)。
%   模型:bm_j = g_H · A · (V_j .* d)  (= 各點 g_H S_i V_j d 的堆疊,含增益,all-source)
%   目標:bf_j = exc_sign(j)·(−Bstack_j)  (all-source 物理場;Bstack = −FEM, actuator 框)
%   回傳全域相對 RMSE [%](對 R_select 球內所有點 × 6 sim)。
%   與 solve_d 的法方程一致:d 最小化 Σ‖g_H A diag(V_j) d − bf_j‖²。
    num = 0; den = 0;                                        % 相對 RMSE 的分子/分母
    for j = 1:N_I                                            % 逐模擬
        bm = gH * (A * (Vmat(:,j) .* d));                   % 3Np x 1 sensor 模型場
        bf = exc_sign(j) * (-Bstack(:,j));                  % 3Np x 1 all-source FEM 場
        num = num + sum((bm - bf).^2);
        den = den + sum(bf.^2);
    end
    nrmse_sensor = sqrt(num/den) * 100;                     % 全域相對 RMSE [%]
end
