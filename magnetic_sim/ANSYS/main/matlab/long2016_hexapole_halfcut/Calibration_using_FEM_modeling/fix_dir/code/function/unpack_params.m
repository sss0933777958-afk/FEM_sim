function [ell, gB, Khat] = unpack_params(x, freemask)
%UNPACK_PARAMS  Unpack optimizer vector into model parameters.
%   [ell, gB, Khat] = UNPACK_PARAMS(x, freemask)
%     x        : packed vector [ell_mm; gB; Khat(free entries)]
%     freemask : 6x6 logical, true where Khat is free (Khat(1,1) fixed = 5/6)
%   Returns ell [m], gB (scalar), Khat (6x6, with Khat(1,1)=5/6).
    ell  = x(1)*1e-3;          % mm -> m
    gB   = x(2);
    Khat = zeros(6);
    Khat(1,1)       = 5/6;     % gauge: k11 fixed to break the degeneracy
    Khat(freemask)  = x(3:end);
end
