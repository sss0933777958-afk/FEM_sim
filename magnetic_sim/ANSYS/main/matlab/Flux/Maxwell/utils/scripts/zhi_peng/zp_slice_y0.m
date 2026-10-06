%% zp_slice_y0.m -- Zhi-peng pair pole, P1 excitation: magnetic-circuit arrows on the plane y = 0 around the centre.
%  Same construction as the 10-01 / 10-02 circuit figures (pair_slice.m, ntu_circuit.m): a jittered grid of
%  ~NARR arrow positions over the window, field interpolated from the y = 0 .fld nodes (scatteredInterpolant),
%  arrow length ~ (|b|/max)^0.25 in the x-z direction.
%  Frame: origin = centre between the tips, Maxwell global (0, 0, 0.2886751) mm; x, z in um.
%  Output: utils/data/zhi_peng/zp_slice_y0.mat
%    Xs, Zs, Uq, Wq, Bm  arrows [um] and |b| [mT] (field as exported, 1 A excitation);  C0, WIN, NARR, JIT, f
clearvars;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT  = fullfile(CALROOT, 'utils', 'data', 'zhi_peng');
cdir = fullfile(CALROOT, 'config', 'zhi_peng', 'Pair_pole');
addpath(cdir, '-begin');  clear mt_constants;  cfg = mt_constants();  rmpath(cdir);

C0   = [0, 0, 288.6751];                              % centre [um] (Maxwell global)
WIN  = 600;  MRG = 250;                               % window half width, source margin [um]
NARR = 2000;  JIT = 0.9;  LEN = 20;                   % arrows, jitter fraction, max arrow length [um]

f = fullfile(cfg.fld_dir, cfg.fld_files{1});          % P1 excitation
fid = fopen(f, 'r');  fgetl(fid); fgetl(fid);
C = textscan(fid, '%f %f %f %f %f %f', 'CollectOutput', true);  fclose(fid);
D = C{1};  clear C
D(:,1:3) = D(:,1:3) * 1e3;                            % mm -> um
X = D(:,1) - C0(1);  Z = D(:,3) - C0(3);
src = abs(D(:,2)) < 1e-6 & abs(X) < WIN + MRG & abs(Z) < WIN + MRG;   % y = 0 nodes around the window
B = cfg.s_source(1) * 1e3 * D(src, 4:6);             % mT (field as exported, all-source)
X = X(src);  Z = Z(src);  clear D
Fx = scatteredInterpolant(X, Z, B(:,1), 'linear', 'none');
Fz = Fx;  Fz.Values = B(:,3);
Fm = Fx;  Fm.Values = vecnorm(B, 2, 2);

cs  = sqrt((2*WIN)^2 / NARR);  ng = round(2*WIN / cs);  h = 2*WIN / ng;
[gx, gz] = meshgrid(0:ng-1, 0:ng-1);  rng(0);
Xs = -WIN + (gx(:) + 0.5 + JIT*(rand(numel(gx),1) - 0.5)) * h;
Zs = -WIN + (gz(:) + 0.5 + JIT*(rand(numel(gz),1) - 0.5)) * h;
Bx = Fx(Xs, Zs);  Bz = Fz(Xs, Zs);  Bm = Fm(Xs, Zs);
ok = isfinite(Bx) & isfinite(Bz) & isfinite(Bm) & Bm > 1e-4;
Xs = Xs(ok);  Zs = Zs(ok);  Bx = Bx(ok);  Bz = Bz(ok);  Bm = Bm(ok);
bxz = hypot(Bx, Bz);  bxz(bxz == 0) = 1e-12;
scl = LEN * (Bm ./ max(Bm)).^0.25 ./ bxz;
Uq  = Bx .* scl;  Wq = Bz .* scl;
save(fullfile(DAT, 'zp_slice_y0.mat'), 'Xs', 'Zs', 'Uq', 'Wq', 'Bm', 'C0', 'WIN', 'NARR', 'JIT', 'LEN', 'f');
fprintf('y = 0 circuit arrows: %d (from %d nodes), |b| %.2f .. %.1f mT\n', numel(Xs), numel(X), min(Bm), max(Bm));
