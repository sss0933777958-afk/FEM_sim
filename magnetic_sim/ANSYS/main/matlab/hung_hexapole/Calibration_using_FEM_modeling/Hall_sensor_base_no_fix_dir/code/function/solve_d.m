function d = solve_d(Vmat, exc_sign, M, c, N_I)
% SOLVE_D  解每極常數 d（無增益 no-gain,與 Hall_sensor_base_fix_dir 對齊;無 g_H、無 K_H）。
%   模型:b_ij = S_i V_j d;閉式解 d = (Σ_j V_j M V_j)^{-1}(Σ_j V_j cc_j)。
%   M、c 重用 page-1 在 ell_hat 下算好的(c = Σ S^T B,B 為頁1 負號版場);
%   cc_j = exc_sign(j)·(−c_j) = all-source 物理場。唯一與 fix_dir 不同處 = M,c 來自 18-param bias。
    A = zeros(6,6); rhs = zeros(6,1);                         % 法矩陣與右端累加器
    for j = 1:N_I                                             % 逐模擬
        Vj  = diag(Vmat(:,j));                                % V_j = diag(該 sim 的 6 個 all-source sensor 電壓)
        ccj = exc_sign(j) * (-c(:,j));                        % all-source cc_j(b 同步翻:−c=物理,×exc_sign=all-source)
        A   = A   + Vj * M * Vj;                              % Σ V_j M V_j
        rhs = rhs + Vj * ccj;                                 % Σ V_j cc_j
    end
    d = A \ rhs;                                             % no-gain d（不再除 g_H）
end
