%% phi_max.m -- flux of every pole when one coil carries I_max (single-parameter calibration, long2016).
%  flux_force_unit.pdf Sec. 1.3-1.4, 5.3:  Phi = (B_ghat_I / B_ghat_Phi) Kbar_I I,   B_ghat_Phi = 1 / (4 pi l_hat^2).
%  Units: B_ghat_I [mT/A], l_hat [um] -> B_ghat_Phi [1/um^2] = [mT/(mT um^2)];  1 mT um^2 = 1e-15 Wb = 1e-6 nWb.
%  Coil j alone at I_max (I = I_max e_j):  Phi(:, j) = (B_ghat_I / B_ghat_Phi) I_max Kbar_I(:, j)  [nWb], pole order P1..P6.
%  Output utils/data/long2016_hexapole_halfcut/phi_max.mat: I_max [A], B_ghat_I [mT/A], l_hat [um], B_ghat_Phi [mT/nWb],
%    Kbar_I (6 x 6), Phi (6 x 6 [nWb], column j = coil j alone at I_max), Phi_self (6 x 1 [nWb], diag(Phi))
if ~exist('I_max', 'var'), I_max = 3.5; end                  % [A]
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
cal = load_flux_calib('long2016_hexapole_halfcut', 'tip40um', 'calib_current_maxwell_axshN43_R150_single', 'current');
B_ghat_I = cal.gB_m;  l_hat = cal.l_m;  Kbar_I = cal.Mbar_m;
B_ghat_Phi = 1 / (4*pi*l_hat^2) * 1e6;                         % [mT/nWb]  (1 nWb = 1e6 mT um^2)
Phi = (B_ghat_I / B_ghat_Phi) * I_max * Kbar_I;                % [nWb]
Phi_self = diag(Phi);
save(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', 'phi_max.mat'), ...
     'I_max', 'B_ghat_I', 'l_hat', 'B_ghat_Phi', 'Kbar_I', 'Phi', 'Phi_self');
fprintf('B_ghat_I %.5g mT/A, l_hat %.5g um, B_ghat_Phi %.5g mT/nWb, ratio %.5g nWb/A\n', B_ghat_I, l_hat, B_ghat_Phi, B_ghat_I / B_ghat_Phi);
disp('Phi [nWb] (row = pole P1..P6, column = coil P1..P6 alone at I_max):');  disp(round(Phi, 1));
