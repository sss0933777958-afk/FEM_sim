function [gF, KH, Ra, H] = compute_KH(dd, Vmat, F, mu0, Nc)
% COMPUTE_KH  依 derivation.pdf 流程,從 sensor 的 d 求 ĝ_F、K̂_H、R_a。
%   磁荷 Q̂_j = (1/μ0) diag(dd) V_j ;  Q̂ = ĝ_F K̂_H F
%   ⇒ ĝ_F K̂_H = (1/μ0) D V F^T (F F^T)^{-1} = H
%   gauge (K̂_H)_11 = 5/6 ⇒ ĝ_F = (6/5) h11, K̂_H = H/ĝ_F, R_a = N_c/(μ0 ĝ_F)。
%   (ported from gen_KH_latex.m loop body)
    G  = (1/mu0) * diag(dd) * Vmat;                   % 6x6 磁荷矩陣 Q̂(列=極 P1..P6,欄=sim/coil)
    H  = G * F.' / (F * F.');                          % H = Q̂ F^T (F F^T)^{-1}(欄重排成激發極 P1..P6)
    gF = (6/5) * H(1,1);                               % ĝ_F = (6/5) h11
    KH = H / gF;                                       % K̂_H = H/ĝ_F((1,1)=5/6)
    Ra = Nc / (mu0 * gF);                              % R_a = N_c/(μ0 ĝ_F)
end
