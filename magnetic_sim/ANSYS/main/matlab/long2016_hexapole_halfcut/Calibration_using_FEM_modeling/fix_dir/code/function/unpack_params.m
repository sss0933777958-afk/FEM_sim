function [ell, ghat_I_B, K_bar] = unpack_params(x, freemask)
%UNPACK_PARAMS  Unpack optimizer vector into model parameters.
%   [ell, ghat_I_B, K_bar] = UNPACK_PARAMS(x, freemask)
%     x        : packed vector [ell_mm; ghat_I_B; K_bar(free entries)]
%     freemask : 6x6 logical, true where K_bar is free (K_bar(1,1) fixed = 5/6)
%   Returns ell [m], ghat_I_B (= ^B g_I [mT/A]), K_bar (= K-bar, 6x6, with K_bar(1,1)=5/6).
    ell      = x(1)*1e-3;          % mm -> m（fit 在 SI；呼叫端再 ×1e6 成 µm）
    ghat_I_B = x(2);               % ^B g_I (current-side gain, mT/A)
    K_bar    = zeros(6);
    K_bar(1,1)      = 5/6;         % gauge: k11 fixed to break the degeneracy
    K_bar(freemask) = x(3:end);
end
