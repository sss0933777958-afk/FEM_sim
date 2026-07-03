function A = build_A(ell, Pc, P)
%BUILD_A  Stacked nondimensional spatial-function matrix A (3Np x 6).
%   A = BUILD_A(ell, Pc, P)
%   Column k holds the kernel (pbar - Pc_k)/||pbar - Pc_k||^3 stacked over the Np
%   sample points (pbar = P/ell, actuator frame).  (no_fix_l.pdf step 3b.)
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
