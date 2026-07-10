function J = sensor_residual_bias(P, B, Vmat, exc_sign, ell_hat, Pc, d, N_I)
% SENSOR_RESIDUAL_BIAS  sensor 模型成本 J（actuator 框,build_S 逐點版,no-gain）。
%   模型:bm_ij = S_i V_j d            (V_j = diag(Vmat(:,j))；無 g_H)
%   目標:bf_ij = exc_sign(j)·(−B(i,:,j)) (all-source 物理場;B = 頁1 負號版 −FEM、actuator 框)
%   回傳 cost J = Σ_ij‖bm_ij − bf_ij‖²  [T²]（與 Hall_sensor_base_fix_dir 對齊:只回 cost,無 RMSE/g_H）。
%   唯一與 fix_dir 不同 = 點/電荷在 actuator 框、Pc = Pc_18(18-param bias)。 (calls build_S)
    Np = size(P,1);                                         % 工作點數
    J  = 0;                                                 % 成本累加器
    for j = 1:N_I                                           % 逐模擬
        Vj = diag(Vmat(:,j));                              % V_j = diag(該 sim 的 6 個 all-source sensor 電壓)
        for i = 1:Np                                       % 逐工作點
            Si = build_S(P(i,:), ell_hat, Pc);            % 3×6 空間核（actuator 框、Pc_18）
            bm = Si * Vj * d;                             % sensor 模型場 b_model = S_i V_j d（no-gain）
            bf = exc_sign(j) * (-squeeze(B(i,:,j)).');    % all-source 物理場 b_FEM
            J  = J + sum((bm - bf).^2);                   % 累加殘差平方
        end
    end
end
