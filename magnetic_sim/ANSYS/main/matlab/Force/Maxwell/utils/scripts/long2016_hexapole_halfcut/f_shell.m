% f_shell.m -- force shell at a point P_UM [um] in the actuator frame (default: the origin = common-sphere centre of
% the six pole tips).  TAG names the output ('' for the origin).
% Current-based force model  f_k = I' * Q_k * I,  Q_k = F H_I' * L_k(p/l_hat, e) * F H_I  (6x6 symmetric),
% F H_I = sqrt(F g_I) * K_I bar,  F g_I = U_F * M g_B / (2 l_hat) * (B g_I)^2.
% For every direction r_hat = [cos(phi)cos(theta), cos(phi)sin(theta), sin(phi)] (phi elevation, theta azimuth):
%   Q_r = sum_k r_hat_k Q_k,  eigen-decomposition -> lambda_max and its unit eigenvector I* (||I*||_2 = 1),
%   then the full force f = [I*' Q_1 I*, I*' Q_2 I*, I*' Q_3 I*]  (r_hat' f = lambda_max).
% Flux parameters: Flux/Maxwell calib_current_maxwell_axshN55_R150_eighteen (long2016, eighteen parameters).
% Output: utils/data/long2016_hexapole_halfcut/f_shell[_<TAG>].mat  (also P_UM)
%   PHI, TH [deg] grids (NPHI x NTH); F [pN] NPHI x NTH x 3; LAM [pN] lambda_max; IS 6 x NPHI x NTH eigenvectors;
%   Q 6x6x3 [pN/A^2]; l_hat [um], e, gB, Kbar, gF, mgB, UF, src;
%   R_act 3x3 (global -> actuator, from model_config; a force in the global (Maxwell) frame is R_act' * f)
if ~exist('P_UM','var'), P_UM = [0 0 0]; end
if ~exist('TAG','var'),  TAG  = ''; end
clearvars -except P_UM TAG;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
MATNAME = 'calib_current_maxwell_axshN55_R150_eighteen';
mgB = 0.0451;  UF = 1e3;                    % M g_B [A um^2/mT], unit factor to pN (package defaults)
DSTEP = 0.25;                               % direction grid step [deg]

cal = load_flux_calib(MODEL, GEOM, MATNAME, 'current');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;  % F g_I [pN/A^2]
FH  = sqrt(gF) * cal.Mbar_m;                % F H_I hat [sqrt(pN)/A]
L   = build_L(P_UM(:).', cal.l_m, cal.e_m, cal.Pc_base);   % 1 x 6 x 6 x 3 at P_UM (actuator frame, um)
Q   = zeros(6, 6, 3);
for k = 1:3
    Lk = squeeze(L(1,:,:,k));  Q(:,:,k) = FH.' * Lk * FH;  Q(:,:,k) = (Q(:,:,k) + Q(:,:,k).')/2;
end

phi = -90:DSTEP:90;  th = 0:DSTEP:360;      % theta = 360 repeats 0 so the shell closes
[TH, PHI] = meshgrid(th, phi);
NP = numel(phi);  NT = numel(th);
F = zeros(NP, NT, 3);  LAM = zeros(NP, NT);  IS = zeros(6, NP, NT);
for a = 1:NP
    for b = 1:NT
        r  = [cosd(phi(a))*cosd(th(b)); cosd(phi(a))*sind(th(b)); sind(phi(a))];
        Qr = r(1)*Q(:,:,1) + r(2)*Q(:,:,2) + r(3)*Q(:,:,3);
        [V, D] = eig((Qr + Qr.')/2);
        [lam, im] = max(diag(D));
        v = V(:, im);  if v(find(abs(v) == max(abs(v)), 1)) < 0, v = -v; end   % sign fixed (f is even in I)
        LAM(a,b) = lam;  IS(:,a,b) = v;
        for k = 1:3, F(a,b,k) = v.' * Q(:,:,k) * v; end
    end
end
chk = max(abs(sum(F .* cat(3, cosd(PHI).*cosd(TH), cosd(PHI).*sind(TH), sind(PHI)), 3) - LAM), [], 'all');
fprintf('P = %s um (actuator) | l_hat %.2f um | B g_I %.4f mT/A | F g_I %.5g pN/A^2 | grid %d x %d\n', mat2str(P_UM), cal.l_m, cal.gB_m, gF, NP, NT);
fprintf('lambda_max: min %.5g, max %.5g pN | |f|: min %.5g, max %.5g pN | max|r''f - lambda| %.2g\n', ...
        min(LAM(:)), max(LAM(:)), min(vecnorm(F,2,3),[],'all'), max(vecnorm(F,2,3),[],'all'), chk);
fprintf('eig(Q_k) k=1..3:\n');  for k = 1:3, fprintf('  %s\n', mat2str(eig(Q(:,:,k)).', 4)); end
R_act = cal.R_act;  l_hat = cal.l_m;  e = cal.e_m;  gB = cal.gB_m;  Kbar = cal.Mbar_m;  src = cal.src;
SFX = '';  if ~isempty(TAG), SFX = ['_' TAG]; end
save(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', ['f_shell' SFX '.mat']), 'P_UM', 'PHI', 'TH', 'F', 'LAM', 'IS', 'Q', 'DSTEP', 'R_act', ...
     'l_hat', 'e', 'gB', 'Kbar', 'gF', 'mgB', 'UF', 'src', 'MODEL', 'MATNAME');
