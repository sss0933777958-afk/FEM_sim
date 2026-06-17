function S = build_S(p_row, ell, Pc)
% BUILD_S  步驟6 的單點空間函數矩陣 S(pbar; ell)。
%   p_row : 1x3 工作點 p_i [m]
%   ell   : 純量特徵長度 ell_hat [m]
%   Pc    : 3x6 電荷位置矩陣(單位距離的極方向)
%   S     : 3x6,第 k 欄 = (pbar - Pc(:,k)) / ||pbar - Pc(:,k)||^3
%   (ported verbatim from calib_fem.m local function)
    pbar = p_row(:) / ell;                                    % 3x1 正規化點 pbar_i = p_i/ell
    D    = pbar - Pc;                                         % 3x6 各電荷的差 (pbar - pbar_ck)(廣播相減)
    nrm  = vecnorm(D);                                        % 1x6 各欄歐氏範數 ||pbar - pbar_ck||
    S    = D ./ (nrm.^3);                                     % 3x6 每欄除以自己範數的立方 → s_ik
end
