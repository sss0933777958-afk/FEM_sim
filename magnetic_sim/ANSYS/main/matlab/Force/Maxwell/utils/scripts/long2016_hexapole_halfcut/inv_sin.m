%% inv_sin.m -- six fluxes of the 2D Fourier inverse model for a sinusoidal force at one position.
%  Force command: f(t) = [F0 sin(2 pi FR t); 0; 0] pN at p = PX (actuator frame, um).
%  Inverse model (ls_fit.m DIR f2d, PH6, fourier2d_new.pdf Sec. 4.4-4.5): unit-force flux
%    Phi_hat_h(p, phi, theta) = P' D_h(phi, theta) P,  P = [pbar; 1],  pbar = p / l_hat,
%    D_h = sum_m sum_n [A_h,mn cos(m phi) + B_h,mn sin(m phi)] g_mn(theta), g_mn = cos(n theta) (m even) / sin(n theta) (m odd),
%    phi = azimuth from +x, theta = angle from +z.  Fits without the DEF field use the old definition (fourier2d_old.pdf:
%    (X, X', g_n) = (A, Ap, cos n phi_lat) for m + n even, (B, Bp, sin n phi_lat) for m + n odd, theta = azimuth).
%  The force is quadratic in Phi, so a force of size |f| along r_hat needs Phi = sqrt(|f| / 1 pN) Phi_hat(p, r_hat),
%  r_hat = f / |f| (+x_a for f_x > 0, -x_a for f_x < 0; Phi = 0 at f = 0).
%  Check: f_hat_k = Phi' L_k(pbar) Phi against the command.
%  Input  utils/data/long2016_hexapole_halfcut/<FIT>.mat (A, B, DEF, M, N, l_hat; old fits A, Ap, B, Bp)
%  Output utils/data/long2016_hexapole_halfcut/<OUTN>.mat (default inv_sin.mat); t = 0:DT:TEND
%    t (1 x N_t [s]), f_cmd (3 x N_t [pN]), Phi (6 x N_t, pole order P1..P6; sqrt(|f| / 1 pN) Phi_hat, gain-free so
%    that f = Phi' L_k Phi), Phi_nWb (6 x N_t [nWb], physical flux sqrt(|f| / F_ghat_Phi) Phi_hat), f_hat (3 x N_t [pN]),
%    e_f (1 x N_t [%])
if ~exist('FIT', 'var'), FIT = 'fit_f2d6_single_lam0_M1N1'; end   % M = N = 1, LAM = 0, single
if ~exist('PX', 'var'), PX = [0 0 0]; end                          % position [um], actuator frame
if ~exist('F0', 'var'), F0 = 5; end                                % amplitude [pN]
if ~exist('FR', 'var'), FR = 100; end                              % frequency [Hz]
if ~exist('TEND', 'var'), TEND = 0.02; end                         % duration [s] (two periods)
if ~exist('DT', 'var'), DT = 1e-5; end                            % sample time [s] (1 sample = 10 us)
if ~exist('OUTN', 'var'), OUTN = 'inv_sin'; end                    % output name (no .mat)
% VCAL: voltage calibration of Flux/Maxwell (no .mat, e.g. 'calib_voltage_maxwell_axshN43_soff3mm_R150_single') -> also the
%   coil voltages V_d = (B_ghat_Phi / B_ghat_V) Mbar^-1 Phi_nWb [mV] (Mbar = D_bar, B_ghat_V = gV_hat [mT/mV]); '' = none
if ~exist('VCAL', 'var'), VCAL = ''; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
DAT = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut');
vn_ = who('-file', fullfile(DAT, [FIT '.mat']));
F = load(fullfile(DAT, [FIT '.mat']), 'A', 'B', 'M', 'N', 'l_hat', vn_{ismember(vn_, {'Ap', 'Bp', 'DEF', 'CAL'})});
cal = load_flux_calib('long2016_hexapole_halfcut', 'tip40um', 'calib_current_maxwell_axshN43_R150_single', 'current');
assert(abs(cal.l_m - F.l_hat) < 1e-9, 'inv_sin:cal', 'fit and calibration l_hat differ');
P = [PX(:) / F.l_hat; 1];                                          % P = [pbar; 1]

t = (0:round(TEND/DT)) * DT;  NT = numel(t);
f_cmd = [F0 * sin(2*pi*FR*t); zeros(2, NT)];
Phi = zeros(6, NT);
for i = 1:NT
    fn = norm(f_cmd(:,i));
    if fn < 1e-12, continue, end
    r_hat = f_cmd(:,i) / fn;
    Phi(:,i) = sqrt(fn) * phihat_(F, P, r_hat);
end
L = build_L(PX(:).', cal.l_m, cal.e_m, cal.Pc_base);                % 1 x 6 x 6 x 3
f_hat = zeros(3, NT);
for k = 1:3
    L_k = squeeze(L(1,:,:,k));  L_k = (L_k + L_k.') / 2;
    f_hat(k,:) = sum(Phi .* (L_k * Phi), 1);
end
e_f = vecnorm(f_hat - f_cmd, 2, 1) ./ max(vecnorm(f_cmd, 2, 1), eps) * 100;
e_f(vecnorm(f_cmd, 2, 1) < 1e-12) = 0;
M = F.M;  N = F.N;
% physical flux (flux_force_unit.pdf Sec. 1.3, 2, 3): Phi_phys = sqrt(||f_d|| / F_ghat_Phi) Phi_hat,
% F_ghat_Phi = M_ghat_B / (2 l_hat) (B_ghat_Phi)^2,  B_ghat_Phi = 1 / (4 pi l_hat^2)
% units: l_hat [um], M_ghat_B [A um^2/mT], Phi [mT um^2] = 1e-15 Wb; UF converts A mT um to pN (as nn_train.m)
M_ghat_B = 0.0451;  UF = 1e3;
B_ghat_Phi = 1 / (4*pi*cal.l_m^2);                                 % [1/um^2]
F_ghat_Phi = UF * M_ghat_B / (2*cal.l_m) * B_ghat_Phi^2;           % [pN / (mT um^2)^2]
Phi_nWb = Phi / sqrt(F_ghat_Phi) * 1e-15 * 1e9;                    % [nWb]
save(fullfile(DAT, [OUTN '.mat']), 't', 'f_cmd', 'Phi', 'Phi_nWb', 'F_ghat_Phi', 'B_ghat_Phi', 'f_hat', 'e_f', 'FIT', 'PX', ...
     'F0', 'FR', 'M', 'N', 'DT');
if ~isempty(VCAL)                                                  % coil voltages from the flux (voltage calibration)
    VC = load(fullfile(fileparts(fileparts(CALROOT)), 'Flux', 'Maxwell', 'data', 'long2016_hexapole_halfcut', '.mat', [VCAL '.mat']), ...
              'D_bar', 'gV_hat');
    Mbar = VC.D_bar;  B_ghat_V = VC.gV_hat;                         % [-], [mT/mV]
    B_ghat_Phi_mTnWb = B_ghat_Phi * 1e6;                            % [1/um^2] -> [mT/nWb] (1 mT um^2 = 1e-6 nWb)
    V_d  = (B_ghat_Phi_mTnWb / B_ghat_V) * (Mbar \ Phi_nWb);        % [mV]
    save(fullfile(DAT, [OUTN '.mat']), 'V_d', 'Mbar', 'B_ghat_V', 'B_ghat_Phi_mTnWb', 'VCAL', '-append');
    fprintf('V_d [mV] at the last sample: %s (B_ghat_Phi/B_ghat_V = %.4f mV/nWb)\n', mat2str(V_d(:,end).', 5), B_ghat_Phi_mTnWb/B_ghat_V);
end
fprintf('F_ghat_Phi = %.4g pN/(mT um^2)^2, Phi_nWb range [%.2f, %.2f] nWb\n', F_ghat_Phi, min(Phi_nWb(:)), max(Phi_nWb(:)));
fprintf('%s at p = %s um: max |e_f| %.4f %%, Phi range [%.3f, %.3f] -> %s\n', FIT, mat2str(PX), max(e_f), ...
        min(Phi(:)), max(Phi(:)), fullfile(DAT, [OUTN '.mat']));
fprintf('Phi_hat(+x_a) = %s\nPhi_hat(-x_a) = %s\n', mat2str(phihat_(F, P, [1;0;0]).', 4), mat2str(phihat_(F, P, [-1;0;0]).', 4));

function z = phihat_(F, P, r_hat)
% unit-force flux P' D_h(phi, theta) P of the fitted 2D Fourier model (6 x 1)
Nk = size(F.A, 5);  z = zeros(Nk, 1);
if isfield(F, 'DEF')                                               % fourier2d_new.pdf
    ph = atan2(r_hat(2), r_hat(1));  th = acos(max(-1, min(1, r_hat(3))));
    for k = 1:Nk
        D = zeros(4);
        for m = 0:F.M
            for n = 0:F.N
                if mod(m, 2) == 0, g = cos(n*th);  else, g = sin(n*th); end
                D = D + (cos(m*ph) * F.A(:,:,m+1,n+1,k) + sin(m*ph) * F.B(:,:,m+1,n+1,k)) * g;
            end
        end
        z(k) = P.' * D * P;
    end
    return
end
th = atan2(r_hat(2), r_hat(1));  ph = asin(max(-1, min(1, r_hat(3))));   % old definition (fourier2d_old.pdf)
for k = 1:Nk
    D = zeros(4);
    for m = 0:F.M
        for n = 0:F.N
            if mod(m + n, 2) == 0, X = F.A;  Xp = F.Ap;  g = cos(n*ph);  else, X = F.B;  Xp = F.Bp;  g = sin(n*ph); end
            D = D + (cos(m*th) * X(:,:,m+1,n+1,k) + sin(m*th) * Xp(:,:,m+1,n+1,k)) * g;
        end
    end
    z(k) = P.' * D * P;
end
end
