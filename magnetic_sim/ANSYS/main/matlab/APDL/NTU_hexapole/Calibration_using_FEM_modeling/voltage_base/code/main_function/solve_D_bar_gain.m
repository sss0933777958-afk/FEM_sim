function [D_bar, ghat_V_B, H_V] = solve_D_bar_gain(Dv, Vmat)
%SOLVE_D_BAR_GAIN  Back out D-bar and the voltage-base gain ^Bg_V from G and V.
%   [D_bar, ghat_V_B, H_V] = SOLVE_D_BAR_GAIN(Dv, Vmat)
%   Recover Ĥ_V = G Vᵀ(VVᵀ)⁻¹ (code Dv = G = D^v), then gauge with D_bar(1,1)=5/6:
%     ^Bg_V = (6/5) H_V(1,1),  D_bar = (5/(6 H_V(1,1))) H_V.
%   (paper step 9.  Voltage-side twin of current-side solve_KI_bar_gain; new 2026-07-16.
%   H_V (= Ĥ_V, the un-gauged recovery) is returned for internal use only, not printed.)
    H_V      = (Dv * Vmat.') / (Vmat * Vmat.');   % Ĥ_V = G Vᵀ(VVᵀ)⁻¹
    D_bar    = H_V * (5 / (6 * H_V(1,1)));         % gauge D_bar(1,1)=5/6
    ghat_V_B = (6/5) * H_V(1,1);                   % ^Bg_V (voltage-side gain, mT/mV)
end
