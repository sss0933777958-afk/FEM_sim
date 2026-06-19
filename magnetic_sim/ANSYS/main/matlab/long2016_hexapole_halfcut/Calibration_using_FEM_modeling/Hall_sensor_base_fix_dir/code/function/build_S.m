function S = build_S(p_row, ell, Pc)
% BUILD_S  建構「單一工作點」的點電荷空間函數矩陣 S(p̄; ℓ̂)。
% -------------------------------------------------------------------------
% 物理背景（Long 2016 等效磁荷模型；本包為 Hall-sensor 版）：
%   把每顆磁極等效成一顆放在「極尖端方向、距 WP 中心 ℓ̂」處的磁單極(magnetic charge)。
%   某點 p 的磁場 = 6 顆電荷的庫倫式場線性疊加：
%       b(p) = S(p) · w           （w 為 6×1 每極電荷「強度」）
%   本函式只負責算「與強度 w 無關」的純幾何核 S：
%       S(:,k) = (p̄ − p̄_ck) / ‖p̄ − p̄_ck‖³，  p̄ = p/ℓ̂，p̄_ck = Pc(:,k)（第 k 極方向）
%   S 是線性模型的「設計矩陣」。本包(Hall-sensor 模型)的強度為
%       w = g_H · V_j · d         （g_H = 1/(4πℓ̂²)、V_j = sensor 電壓、d = 每極常數）
%   故第 j 次模擬的場 b_ij = g_H · S_i · V_j · d。
%   注意：本包「沒有 gB」——gB 是 fix_dir charge-fit 的增益，不屬於此 Hall-sensor 模型。
%
% 在 pipeline 的位置：
%   main.m PAGE1 步驟6 逐點呼叫，用來組法矩陣 M=ΣSᵀS、右端 c=ΣSᵀb；
%   sensor 模型(solve_d / sensor_residual)也重用同一個 S。
%
% 輸入：
%   p_row : 1×3（或 3×1）工作點座標 p_i [m]（WP 框）
%   ell   : 純量，特徵長度 ℓ̂ [m]（電荷離 WP 中心的距離尺度）
%   Pc    : 3×6 電荷位置矩陣，各欄為「單位距離」的極方向 d̂_k（= 極尖單位向量）
% 輸出：
%   S     : 3×6，第 k 欄 = (p̄ − Pc(:,k)) / ‖p̄ − Pc(:,k)‖³ = 單位電荷在 p_i 的場核
%
%   (ported verbatim from calib_fem.m local function)
% -------------------------------------------------------------------------
    pbar = p_row(:) / ell;                                    % 正規化工作點 p̄_i = p_i/ℓ̂（3×1；強制成欄向量）
    D    = pbar - Pc;                                         % 3×6：各電荷的位移 (p̄ − p̄_ck)（pbar 對 Pc 廣播相減）
    nrm  = vecnorm(D);                                        % 1×6：各欄歐氏範數 ‖p̄ − p̄_ck‖（= 電荷到工作點距離）
    S    = D ./ (nrm.^3);                                     % 3×6：每欄除以自己範數立方 → 庫倫式核 s_ik（1/r² 衰減的方向場）
end
