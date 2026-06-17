function nrmse_sensor = sensor_residual(P, B, Vmat, exc_sign, ell_hat, Pc, d, gH, N_I)
% SENSOR_RESIDUAL  sensor 模型殘差(R<=fit 球的點上:b_ij vs g_H S_i V_j d)。
%   回傳全域相對 RMSE [%]。 (ported verbatim from calib_fem.m PAGE 2; calls build_S)
    Np = size(P,1);                                          % 點數
    num = 0; den = 0;                                        % 相對 RMSE 的分子/分母
    for i = 1:Np                                             % 逐點
        Si = build_S(P(i,:), ell_hat, Pc);                  % 3x6(在 ell_hat)
        for j = 1:N_I                                        % 逐模擬
            bm = gH * Si * diag(Vmat(:,j)) * d;             % sensor 模型場 g_H·S_i V_j d(含增益,all-source)
            bf = exc_sign(j) * (-squeeze(B(i,:,j)).');      % all-source FEM 場 b_ij(−B=物理,×exc_sign=all-source)
            num = num + sum((bm-bf).^2); den = den + sum(bf.^2);
        end
    end
    nrmse_sensor = sqrt(num/den)*100;                       % 全域相對 RMSE [%]
end
