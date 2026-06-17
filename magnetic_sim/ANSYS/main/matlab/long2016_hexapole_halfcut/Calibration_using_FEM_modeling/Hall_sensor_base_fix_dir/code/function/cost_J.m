function J = cost_J(ell, P, B, Pc)
% COST_J  步驟7 的剖面化最小二乘成本 J(ell),對所有點與所有模擬加總。
%   ell : 純量 ell_hat [m]
%   P   : Np x 3 工作點
%   B   : Np x 3 x N_I 記錄的場 b_ij
%   Pc  : 3 x 6 電荷位置矩陣
%   J   : 純量成本(各模擬剖面化後殘差平方和之總和)
%   (ported verbatim from calib_fem.m local function costJ; calls build_S)
    Np  = size(P,1);                                          % 點數
    N_I = size(B,3);                                          % 模擬數
    M   = zeros(6,6);                                         % 累加器:sum_i S_i^T S_i(6x6,與 j 無關)
    c   = zeros(6,N_I);                                       % 累加器:sum_i S_i^T b_ij(每模擬 6x1)→ 6xN_I
    bb  = zeros(1,N_I);                                       % 累加器:sum_i b_ij^T b_ij(每模擬純量)
    for i = 1:Np                                              % 逐點 p_i
        Si = build_S(P(i,:), ell, Pc);                       % 該點的 3x6 空間矩陣
        M  = M + Si.' * Si;                                   % 把 S_i^T S_i(6x6)加進共用法矩陣
        for j = 1:N_I                                         % 逐模擬
            bij    = squeeze(B(i,:,j)).';                     % 3x1 該點記錄的場 b_ij
            c(:,j) = c(:,j) + Si.' * bij;                     % 把 S_i^T b_ij(6x1)加進模擬 j 的累加器
            bb(j)  = bb(j) + bij.' * bij;                     % 把 b_ij^T b_ij(純量)加進模擬 j 的累加器
        end
    end
    J = 0;                                                    % 初始化總成本
    for j = 1:N_I                                            % 加總各模擬剖面化後的殘差
        J = J + ( bb(j) - c(:,j).' * (M \ c(:,j)) );         % bb - c^T M^{-1} c(M\c 即解 M g = c)
    end
end
