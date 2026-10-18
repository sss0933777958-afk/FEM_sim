% fmin_ray.m -- radius of the f_min = FLEV contour along every direction of the gamma_map grid (actuator frame).
% For each direction (phi, theta) the first sample interval of utils/data/long2016_hexapole_halfcut/gamma_map.mat where f_min drops below FLEV
% brackets the crossing; it is then refined by bisection to TOL [um], evaluating f_min(p) exactly as gamma_map.m does
% (Fibonacci direction search + fminsearch on the min of lambda_max(Q_r)).  Directions that never drop below FLEV
% inside RMAX keep r = RMAX.
% TAG = 'eighteen' (default) | 'single': which gamma_map[_<TAG>].mat (and so which calibration) is used.
% Output: utils/data/long2016_hexapole_halfcut/fmin_ray_<FLEV*10>[_<TAG>].mat   ph, th [deg]; RL [um] NPH x NTH; FLEV, TOL, R_act
if ~exist('FLEV','var'), FLEV = 6.5; end
if ~exist('TAG','var'), TAG = 'eighteen'; end
clearvars -except FLEV TAG;
SFX = '';  if ~strcmp(TAG, 'eighteen'), SFX = ['_' TAG]; end
close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
G   = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', ['gamma_map' SFX '.mat']));
cal = load_flux_calib(G.MODEL, 'tip40um', G.MATNAME, 'current');
FH  = sqrt(G.gF) * cal.Mbar_m;
TOL = 1;  NDIR = G.NDIR;
k  = (0:NDIR-1).' + 0.5;  zz = 1 - 2*k/NDIR;  az = pi*(1 + sqrt(5))*k;
D  = [sqrt(1 - zz.^2).*cos(az), sqrt(1 - zz.^2).*sin(az), zz].';
opt = optimset('TolX',1e-6, 'TolFun',1e-9, 'MaxFunEvals',400, 'Display','off');
fmin_at = @(p) fmin_point_(p, cal, FH, D, opt);

NPH = numel(G.ph);  NTH = numel(G.th);  RL = zeros(NPH, NTH);  nb = 0;  t0 = tic;
for j = 1:NPH
    for kk = 1:NTH
        if kk == NTH, RL(j,kk) = RL(j,1); continue; end                  % theta = 360 repeats 0
        if abs(G.ph(j)) == 90 && kk > 1, RL(j,kk) = RL(j,1); continue; end   % poles: one direction
        g = G.FMIN(:, j, kk);  i = find(g < FLEV, 1);
        if isempty(i), RL(j,kk) = G.r(end); continue; end
        if i == 1, RL(j,kk) = 0; continue; end
        u  = [cosd(G.ph(j))*cosd(G.th(kk)), cosd(G.ph(j))*sind(G.th(kk)), sind(G.ph(j))];
        lo = G.r(i-1);  hi = G.r(i);                                   % f_min(lo) >= FLEV > f_min(hi)
        while hi - lo > TOL
            m = (lo + hi)/2;
            if fmin_at(m*u) >= FLEV, lo = m; else, hi = m; end
            nb = nb + 1;
        end
        RL(j,kk) = (lo + hi)/2;
    end
end
fprintf('f_min = %g pN: r %.1f .. %.1f um (median %.1f) | %d bisection evaluations, %.0f s\n', ...
        FLEV, min(RL(:)), max(RL(:)), median(RL(:)), nb, toc(t0));
ph = G.ph;  th = G.th;  R_act = G.R_act;
save(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', sprintf('fmin_ray_%d%s.mat', round(10*FLEV), SFX)), 'ph', 'th', 'RL', 'FLEV', 'TOL', 'R_act');

% ---------------------------------------------------------------------------------------
function f = fmin_point_(p, cal, FH, D, opt)
% f_min at one point p [um, actuator frame]: same procedure as gamma_map.m
    L = build_L(p(:).', cal.l_m, cal.e_m, cal.Pc_base);
    Q = zeros(6,6,3);
    for k = 1:3, Lk = squeeze(L(1,:,:,k));  Q(:,:,k) = FH.' * Lk * FH;  Q(:,:,k) = (Q(:,:,k) + Q(:,:,k).')/2; end
    lam = max(real(pageeig(reshape(reshape(Q, 36, 3) * D, 6, 6, []))), [], 1);   % real(): see gamma_map.m
    [l0, a] = min(lam(:));
    x0 = [asin(max(-1, min(1, D(3,a)))); atan2(D(2,a), D(1,a))];
    [~, f] = fminsearch(@(x) lmax_(Q, x), x0, opt);
    f = min(f, l0);
end

function v = lmax_(Q, x)
    rr = [cos(x(1))*cos(x(2)); cos(x(1))*sin(x(2)); sin(x(1))];
    M  = rr(1)*Q(:,:,1) + rr(2)*Q(:,:,2) + rr(3)*Q(:,:,3);
    v  = max(eig((M + M.')/2));
end
