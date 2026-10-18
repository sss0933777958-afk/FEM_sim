function [e, l_hat, J] = fitting(P, Bstack, Pc_base, l0, USE_BIAS)
%FITTING  變數投影優化器：解有效長度 l_hat（+ 偏移 e）。
%   [e, l_hat, J] = FITTING(P, Bstack, Pc_base, l0, USE_BIAS)
%     P       : Np×3 取樣點 [m]
%     Bstack  : 3Np×6 all-source 場 [mT]（欄 = 各激發）
%     Pc_base : 3×6 理想電荷格
%     l0      : l_hat 初值 [m]
%     USE_BIAS: e 開關。false = 只解 l_hat（e≡0）；true = 解 l_hat + e(17)。
%
%   殘差模型（先建 Pc_bar 與 S、再丟優化器）：給 (l,e) → make_Pc(e) → build_S(l,Pc)=S
%   → 每激發欄 LS profile 電荷 g_j=(SᵀS)\(Sᵀb_j) → 殘差 S·g_j − b_j 疊 6 欄（變數投影，
%   電荷被消去、優化器只疊代幾何 l(+e)）。
%
%   回：e(17×1)、l_hat[m]、J（殘差平方和 [mT²]）。l 內部打包成 mm（尺度佳）；不回 Pc。
    if nargin < 5, USE_BIAS = true; end
    opts = optimoptions('lsqnonlin', 'Display','off', ...
        'MaxFunctionEvaluations',1e5, 'MaxIterations',4e3, ...
        'FunctionTolerance',1e-20, 'StepTolerance',1e-12);

    if USE_BIAS
        x0 = [l0*1e3; zeros(17,1)];                                  % l(mm) + e(17)
        xf = lsqnonlin(@(x) resid(x, P, Bstack, Pc_base), x0, [], [], opts);
        l_hat = xf(1)*1e-3;   e = xf(2:18);
    else
        x0 = l0*1e3;                                                 % 只 l(mm)，e 凍結 0
        xf = lsqnonlin(@(x) resid([x; zeros(17,1)], P, Bstack, Pc_base), x0, [], [], opts);
        l_hat = xf(1)*1e-3;   e = zeros(17,1);
    end
    J = sum(resid([l_hat*1e3; e], P, Bstack, Pc_base).^2);
end

% ---- 殘差：g_j 逐激發欄 LS 消去（變數投影）----
function r = resid(x, P, Bstack, Pc_base)
    S = build_S(x(1)*1e-3, make_Pc(x(2:18), Pc_base), P);   % 3Np×6
    M = S.' * S;
    n = size(S, 1);
    r = zeros(numel(Bstack), 1);
    for j = 1:size(Bstack, 2)
        gj = M \ (S.' * Bstack(:,j));
        r((j-1)*n+1 : j*n) = S*gj - Bstack(:,j);
    end
end

% ---- 電荷格 Pc_bar = Pc_base + E(e)（含 e6z 約束）----
function Pc = make_Pc(e17, Pc_base)
    E = zeros(3, 6);
    E(:,1) = e17(1:3);     E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);     E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);      E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);   % e6z 約束（paper step 2c）
    Pc = Pc_base + E;
end

% ---- 每點無因次 Coulomb kernel S（3Np×6）----
function S = build_S(l_hat, Pc, P)
    Np   = size(P, 1);
    pbar = P / l_hat;
    S    = zeros(3*Np, 6);
    for k = 1:6
        d  = pbar - Pc(:,k).';
        r3 = sum(d.^2, 2).^1.5;
        S(:,k) = reshape((d ./ r3).', 3*Np, 1);
    end
end
