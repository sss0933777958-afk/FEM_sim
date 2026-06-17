function A = build_A(ell, Pc, P)
% BUILD_A  堆疊空間函數矩陣 A(3Np x 6):第 k 欄為各點 (pbar-Pc_k)/||pbar-Pc_k||^3 的
%   x,y,z 沿點堆疊(點-major、xyz 交錯),與 Bstack 的堆疊一致。pbar = p/ell。
%   P、Pc 皆在 actuator 框(18-param bias 模型)。
%   (ported verbatim from bias_fit/calib_fem_bias.m local function build_A)
    Np   = size(P, 1);                                       % 點數
    pbar = P / ell;                                          % Np x 3 正規化點
    A    = zeros(3*Np, 6);                                   % 預配置
    for k = 1:6                                              % 逐電荷(欄)
        D  = pbar - Pc(:,k).';                               % Np x 3:pbar - Pc_k(廣播)
        r3 = sum(D.^2, 2).^1.5;                              % Np x 1:||.||^3
        Sk = D ./ r3;                                        % Np x 3:該電荷的核
        A(:,k) = reshape(Sk.', 3*Np, 1);                    % 沿點堆疊成 3Np x 1
    end
end
