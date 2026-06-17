function [d, gH] = solve_d(Vmat, exc_sign, M, c, ell_hat, cnst, N_I)
% SOLVE_D  解 d(依 d_final.pdf,含增益 g_H)。
%   模型:b_ij = g_H · S_i V_j d;閉式解 d = (1/g_H)(Σ_j V_j M V_j)^{-1}(Σ_j V_j cc_j)。
%   g_H = k_m/(ℓ̂²μ_0) = 1/(4πℓ̂²)。M、c 重用 page-1 在 ell_hat 下算好的(c = Σ S^T B,
%   B 為頁1 負號版場);cc_j = exc_sign(j)·(−c_j) = all-source 物理場。
%   回傳 d(=d_final,含增益)與 gH。 (ported verbatim from calib_fem.m PAGE 2)
    gH = cnst.k_m / (ell_hat^2 * cnst.mu_0);                  % 增益 g_H = 1/(4πℓ̂²)  [1/m^2]
    A = zeros(6,6); rhs = zeros(6,1);                         % 法矩陣與右端累加器
    for j = 1:N_I                                             % 逐模擬
        Vj  = diag(Vmat(:,j));                                % V_j = diag(該 sim 的 6 個 all-source sensor 電壓)
        ccj = exc_sign(j) * (-c(:,j));                        % all-source cc_j(b 同步翻:−c=物理,×exc_sign=all-source)
        A   = A   + Vj * M * Vj;                              % Σ V_j M V_j
        rhs = rhs + Vj * ccj;                                 % Σ V_j cc_j
    end
    d = (A \ rhs) / gH;                                       % d_final = (1/g_H)(ΣVMV)^{-1}(ΣVcc) = 4πℓ̂²·(無增益 d)
end
