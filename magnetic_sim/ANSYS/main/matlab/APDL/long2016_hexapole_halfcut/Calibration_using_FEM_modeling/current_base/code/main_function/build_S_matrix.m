function A = build_S_matrix(ell, Pc, P)
%BUILD_S_MATRIX  Stacked per-point S-matrix (nondimensional kernel), 3Np x 6.
%   A = BUILD_S_MATRIX(ell, Pc, P)
%   Column k holds the kernel S_k(p) = (pbar - Pc_k)/||pbar - Pc_k||^3 stacked over
%   the Np sample points (pbar = P/ell, actuator frame).  (paper step 3b.)
%   (Was build_A; renamed 2026-07-16 for clarity — this is the per-point S_matrix.)
    Np   = size(P, 1);
    pbar = P / ell;
    A    = zeros(3*Np, 6);
    for k = 1:6
        D  = pbar - Pc(:,k).';
        r3 = sum(D.^2, 2).^1.5;
        Sk = D ./ r3;
        A(:,k) = reshape(Sk.', 3*Np, 1);
    end
end
