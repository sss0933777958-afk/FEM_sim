%% zp_axis_fit.m -- Zhi-peng pair pole, P1 excitation: 1-D charge fit on 17 axis samples + profile cost.
%  Same sampling and model as the NTU scripts (utils/scripts/NTU_hexapole/sp_axis_calc, sp_axis_fit, sp_cost_b):
%    origin = P1 tip apex at mid-tongue height (config/zhi_peng/Pair_pole), pole axis +x (into P1)
%    samples x = -20 .. -340 um every 20 um on the line y = 0, z = 0 (tip frame)
%    data y = |b_x| [mT] as exported (I = 1 A, 70 turns), model H(x) = a*I/(x - l_hat)^2, a closed form for each l_hat
%  Output: utils/data/zhi_peng/zp_axis.mat    (x, d, b, xc, g, J, r, dc, Hc -- same fields as sp_axis.mat)
%          utils/data/zhi_peng/zp_cost_b.mat  (bs, C, A, bmin, Cmin, amin -- same fields as sp_cost_b.mat)
if ~exist('SAMP','var'), SAMP = 'tip'; end
% SAMP = 'tip' : horizontal line through the P1 tip at mid-tongue height, x = -20 .. -340 um (as NTU)
%        'xa'  : the ORIGINAL 10-01 line (plot/zhi_peng/pair_sampling.m): actuator axis xa through the centre
%                (0,0,0.2887) mm toward the P1 tongue's lower corner, 17 points on +-150 um around the centre.
%                1-D coordinate x = position along xa measured from that corner (+ into P1), so the samples
%                sit at x = -650 .. -350 um; data = |b . u| (field component along the line).
clearvars -except SAMP;  close all;
SFX = '';  if ~strcmp(SAMP, 'tip'), SFX = ['_' SAMP]; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT  = fullfile(CALROOT, 'utils', 'data', 'zhi_peng');
cdir = fullfile(CALROOT, 'config', 'zhi_peng', 'Pair_pole');
addpath(cdir, '-begin');  clear mt_constants;  cfg = mt_constants();  rmpath(cdir);

NP = 17;  I = cfg.I_exc;
switch SAMP
    case 'tip'                                        % horizontal, tip frame (origin = cfg.WP)
        DX = 20e-6;  x = -(1:NP).' * DX;
        u  = [1 0 0];  P = [x, zeros(NP,2)];           % sample positions relative to cfg.WP [m]
    case 'xa'
        C0 = [0, 0, 288.6751e-6];                      % centre (Maxwell global) [m]
        u  = [sqrt(2/3), 0, -1/sqrt(3)];               % xa, centre -> P1 tongue lower corner
        K  = [cfg.P1_TIP_X, 0, cfg.P1_TIP_Z(1)];       % that corner (0.408, 0, 0) mm
        s0 = dot(K - C0, u);                           % corner position along xa (~ 500 um)
        s  = linspace(-150e-6, 150e-6, NP).';  DX = s(2) - s(1);
        x  = s - s0;                                   % 1-D coordinate from the corner [m]
        P  = C0 + s*u - cfg.WP.';                      % positions relative to cfg.WP [m]
    otherwise
        error('unknown SAMP ''%s''', SAMP);
end
d  = -x;                                              % distance from the origin along the line [m]

% ---- field on the samples (P1 excitation) ----
f = fullfile(cfg.fld_dir, cfg.fld_files{1});
fid = fopen(f, 'r');  fgetl(fid); fgetl(fid);
C = textscan(fid, '%f %f %f %f %f %f', 'CollectOutput', true);  fclose(fid);
D = C{1};  D(:,1:3) = D(:,1:3) * 1e-3;                % mm -> m
xs = unique(D(:,1)) - cfg.WP(1);  ys = unique(D(:,2)) - cfg.WP(2);  zs = unique(D(:,3)) - cfg.WP(3);
nx = numel(xs);  ny = numel(ys);  nz = numel(zs);
assert(nx*ny*nz == size(D,1) && D(2,3) ~= D(1,3), 'zp_axis_fit: %s is not a z-fastest full grid', f);
if strcmp(SAMP, 'tip')   % the stencils must stay in front of the tip apex (x < 0 = air)
    ix = arrayfun(@(v) find(xs <= v, 1, 'last'), x);
    assert(all(xs(ix+1) < 0), 'zp_axis_fit: a sample stencil reaches the tip apex');
else                      % xa samples are within 150 um of the centre, inside the iron-free radius
    assert(max(vecnorm(P + cfg.WP.' - C0, 2, 2)) < 423e-6, 'zp_axis_fit: xa samples leave the iron-free ball');
end
b = zeros(NP, 3);
for c = 1:3
    G = permute(reshape(cfg.s_source(1) * 1e3 * D(:,3+c), nz, ny, nx), [3 2 1]);   % T -> mT, all-source
    F = griddedInterpolant({xs, ys, zs}, G, 'linear', 'none');
    b(:,c) = F(P);
end
clear D C

% ---- 1-D fit and profile cost (um units, as in the NTU scripts) ----
xi = x * 1e6;  y = abs(b * u.');                    % |b . u| [mT]
af = @(bb) sum(y ./ (xi-bb).^2) / (I * sum(1 ./ (xi-bb).^4));   % closed-form a for a given l_hat
Jf = @(bb) sum((y - af(bb) * I ./ (xi-bb).^2).^2);
scan = linspace(max(xi)+0.01, 3000, 300001);
[~, im] = min(arrayfun(Jf, scan));
xc = fminbnd(Jf, scan(max(im-1,1)), scan(min(im+1,end)), optimset('TolX',1e-9));
g  = af(xc);  J = Jf(xc);  r = y - g * I ./ (xi-xc).^2;
dc = linspace(0, 400, 801).';  Hc = g * I ./ (-dc - xc).^2;
save(fullfile(DAT, ['zp_axis' SFX '.mat']), 'SAMP', 'u', 'P', 'x', 'd', 'b', 'I', 'DX', 'NP', 'f', 'xc', 'g', 'J', 'r', 'dc', 'Hc');

BR = 500;  if strcmp(SAMP, 'xa'), BR = 2000; end    % sweep 0..BR um (same as NTU 'c500')
bs = linspace(0, BR, 50001).';
C  = arrayfun(Jf, bs);  A = arrayfun(af, bs);
[~, im] = min(C);
bmin = fminbnd(Jf, bs(max(im-1,1)), bs(min(im+1,end)), optimset('TolX',1e-9));
Cmin = Jf(bmin);  amin = af(bmin);
save(fullfile(DAT, ['zp_cost_b' SFX '.mat']), 'bs', 'C', 'A', 'bmin', 'Cmin', 'amin');

btr = vecnorm(b - (b*u.')*u, 2, 2);                   % field component across the sampling line
fprintf('Zhi-peng pair, P1 excitation, SAMP = %s, %d samples at x = %.1f .. %.1f um\n', SAMP, NP, xi(1), xi(end));
fprintf('  |b.u| : %.3f .. %.3f mT (min %.3f), max transverse %.3f mT\n', y(1), y(end), min(y), max(btr));
fprintf('  l_hat = %.4f um  a = %.6g mT*um^2/A  cost = %.6g  RMS = %.4g mT\n', xc, g, J, sqrt(J/NP));
fprintf('  profile: l_hat_min = %.4f um, J/J_min at l_hat = 0: %.4g, at %g um: %.4g\n', bmin, C(1)/Cmin, BR, C(end)/Cmin);
