% gamma_map.m -- force-generation consistency Gamma(p_i) = f_min / f_max inside the design sphere (actuator frame).
% Current-based force model, ||I||_2 <= 1:  Q_k(p) = F H_I' * L_k(p/l_hat, e) * F H_I,  Q_r = sum_k r_hat_k Q_k,
%   f_min(p) = min over r_hat of lambda_max(Q_r)   (the force that can be reached in every direction)
%   f_max(p) = max over r_hat of lambda_max(Q_r),   Gamma = f_min / f_max.
% Positions: spherical grid r = 0:DR:RMAX, phi_p = -90:DA:90, theta_p = 0:DA:360 (actuator frame, um).
% Directions: Fibonacci set of ND unit vectors (coarse search), then fminsearch over (phi, theta) from the best
% coarse direction for both the min and the max (the min of lambda_max often sits on a kink where the two largest
% eigenvalues cross, so a grid alone is first-order accurate there).
% Flux parameters: Flux/Maxwell calib_current_maxwell_axshN55_R150_eighteen (long2016, eighteen parameters).
% TAG = 'eighteen' (default, ...axshN55_R150_eighteen) | 'single' (...axshN43_R150_single, e = 0).
% Output: utils/data/long2016_hexapole_halfcut/gamma_map.mat (eighteen) | gamma_map_single.mat
%   r [um] (NR), ph, th [deg] (NPH, NTH); FMIN, FMAX, GAM  NR x NPH x NTH; R_act; l_hat, gF, mgB, UF, src
if ~exist('DR','var'), DR = 20; end
if ~exist('DA','var'), DA = 5; end
if ~exist('NDIR','var'), NDIR = 2000; end
if ~exist('TAG','var'), TAG = 'eighteen'; end
clearvars -except DR DA NDIR TAG;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
MN = struct('eighteen','calib_current_maxwell_axshN55_R150_eighteen', 'single','calib_current_maxwell_axshN43_R150_single');
MATNAME = MN.(TAG);  SFX = '';  if ~strcmp(TAG, 'eighteen'), SFX = ['_' TAG]; end
mgB = 0.0451;  UF = 1e3;  RMAX = 500;

cal = load_flux_calib(MODEL, GEOM, MATNAME, 'current');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;
FH  = sqrt(gF) * cal.Mbar_m;

r  = 0:DR:RMAX;  ph = -90:DA:90;  th = 0:DA:360;
[RR, PP, TT] = ndgrid(r, ph, th);
P  = [RR(:).*cosd(PP(:)).*cosd(TT(:)), RR(:).*cosd(PP(:)).*sind(TT(:)), RR(:).*sind(PP(:))];
% unique positions only (r = 0, the poles and theta = 360 repeat)
[Pu, ~, iu] = unique(round(P, 9), 'rows');
L  = build_L(Pu, cal.l_m, cal.e_m, cal.Pc_base);            % Nu x 6 x 6 x 3

k  = (0:NDIR-1).' + 0.5;  zz = 1 - 2*k/NDIR;  az = pi*(1 + sqrt(5))*k;
D  = [sqrt(1 - zz.^2).*cos(az), sqrt(1 - zz.^2).*sin(az), zz].';   % 3 x ND, Fibonacci directions
Nu = size(Pu, 1);  fmin = zeros(Nu,1);  fmax = zeros(Nu,1);
opt = optimset('TolX',1e-6, 'TolFun',1e-9, 'MaxFunEvals',400, 'Display','off');
t0 = tic;
for i = 1:Nu
    Q = zeros(6,6,3);
    for kk = 1:3, Lk = squeeze(L(i,:,:,kk));  Q(:,:,kk) = FH.' * Lk * FH;  Q(:,:,kk) = (Q(:,:,kk) + Q(:,:,kk).')/2; end
    Qr  = reshape(reshape(Q, 36, 3) * D, 6, 6, []);
    lam = max(real(pageeig(Qr)), [], 1);  lam = lam(:);   % real(): degenerate pages can come back with round-off imaginary parts, and max() of complex compares |.|
    [~, a] = min(lam);  [~, b] = max(lam);
    g  = @(x) lmax_(Q, x);
    fmin(i) =  fminsearch_val_(g, dir2ang_(D(:,a)), opt);
    fmax(i) = -fminsearch_val_(@(x) -g(x), dir2ang_(D(:,b)), opt);
    fmin(i) = min(fmin(i), lam(a));  fmax(i) = max(fmax(i), lam(b));
    if mod(i, 2000) == 0, fprintf('%d / %d  (%.0f s)\n', i, Nu, toc(t0)); end
end
FMIN = reshape(fmin(iu), size(RR));  FMAX = reshape(fmax(iu), size(RR));  GAM = FMIN ./ FMAX;
fprintf('grid %d x %d x %d (%d unique points) | %.0f s\n', numel(r), numel(ph), numel(th), Nu, toc(t0));
fprintf('origin: f_min %.4f, f_max %.4f, Gamma %.4f | Gamma range %.4f .. %.4f\n', ...
        FMIN(1,1,1), FMAX(1,1,1), GAM(1,1,1), min(GAM(:)), max(GAM(:)));
R_act = cal.R_act;  l_hat = cal.l_m;  src = cal.src;
save(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', ['gamma_map' SFX '.mat']), 'TAG', 'r', 'ph', 'th', 'FMIN', 'FMAX', 'GAM', 'R_act', ...
     'DR', 'DA', 'NDIR', 'l_hat', 'gF', 'mgB', 'UF', 'src', 'MODEL', 'MATNAME');

% ---------------------------------------------------------------------------------------
function v = lmax_(Q, x)
% lambda_max of Q_r for the direction (phi, theta) = x [rad]
    rr = [cos(x(1))*cos(x(2)); cos(x(1))*sin(x(2)); sin(x(1))];
    M  = rr(1)*Q(:,:,1) + rr(2)*Q(:,:,2) + rr(3)*Q(:,:,3);
    v  = max(eig((M + M.')/2));
end

function x = dir2ang_(d)
    x = [asin(max(-1, min(1, d(3)))); atan2(d(2), d(1))];
end

function v = fminsearch_val_(f, x0, opt)
    [~, v] = fminsearch(f, x0, opt);
end
