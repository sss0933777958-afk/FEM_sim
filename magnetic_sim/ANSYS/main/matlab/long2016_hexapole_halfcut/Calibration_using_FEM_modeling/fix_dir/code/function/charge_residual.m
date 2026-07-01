function r = charge_residual(x, coil, dhat, I, freemask)
%CHARGE_RESIDUAL  Stacked residual (B_model - B_FEM) of the point-charge model.
%   r = CHARGE_RESIDUAL(x, coil, dhat, I, freemask)
%   Model (document fix-l, charges fixed on axis at ell*dhat, no bias):
%     B(p) = ^Bg_I * sum_i ( sum_j K_bar_ij I_j ) * (p/ell - dhat_i)/||p/ell - dhat_i||^3
%   Inputs:
%     x        : packed [ell_mm; ghat_I_B; K_bar(free)]  (see unpack_params)
%     coil     : struct array, coil(k).p (Nx3 [m]), .bfem (3N x1 [mT]), .pj (paper idx)
%     dhat     : 3x6 unit pole-tip directions
%     I        : scalar drive current [A]
%     freemask : 6x6 logical of free Khat entries
%   Output r : stacked [Bx;By;Bz] residual over all coils.
    [ell, ghat_I_B, K_bar] = unpack_params(x, freemask);
    r = [];
    for k = 1:numel(coil)
        pn = coil(k).p / ell;  N = size(pn,1);
        B  = zeros(3*N,1);
        w  = ghat_I_B * K_bar(:, coil(k).pj) * I;   % 6x1 charge weights (^Bg_I * K_bar * I)
        for i = 1:6
            dx = pn(:,1)-dhat(1,i); dy = pn(:,2)-dhat(2,i); dz = pn(:,3)-dhat(3,i);
            r3 = (dx.^2 + dy.^2 + dz.^2).^1.5;
            B  = B + w(i) * [dx./r3; dy./r3; dz./r3];
        end
        r = [r; B - coil(k).bfem]; %#ok<AGROW>
    end
end
