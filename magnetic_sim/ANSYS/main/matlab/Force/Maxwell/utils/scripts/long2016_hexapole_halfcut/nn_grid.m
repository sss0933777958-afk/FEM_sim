%% nn_grid.m -- training positions of the NN inverse model: cubic grid in the ACTUATOR frame.
%  Grid step H [um], axes along x_a, y_a, z_a (the frame the network input p/R and rho_solve use),
%  centred on the tip-sphere centre, keeping the nodes with |p| <= R.  The same points are also
%  stored in the measurement frame (P_meas = R_act' * P_act).
%  Output: utils/data/long2016_hexapole_halfcut/nn_grid.mat  (H, R other than 30, 150: nn_grid_h<H>_r<R>.mat,
%          e.g. H = 15, R = 45 -> nn_grid_h15_r45.mat, the 15 um grid inside the 45 um sphere)
%    P_act   N x 3 [um]   actuator frame, origin = sphere centre (the grid)
%    P_meas  N x 3 [um]   measurement frame
%    R_act   3 x 3        rows = P1, P3, P5 tip unit vectors (actuator-frame.md)
%    H, R               grid step and radius [um]
%  NR set: spherical-coordinate layout instead of the grid, uniform density in the ball |p| <= R: the ball is cut into
%  NR bands of equal thickness D = R / NR; band n ((n-1) D <= |p| <= n D) gets a share of the N points proportional to
%  its volume (n^3 - (n-1)^3), placed on one shell at the band's volume midpoint r_n = D ((n^3 + (n-1)^3) / 2)^(1/3);
%  the directions on a shell from a Fibonacci sphere lattice (equal-area, rotated by the golden angle from shell to
%  shell so the shells do not line up)  -> nn_grid_sph_r<R>.mat
%  MID = true (with NR): test positions between the training shells: shells at the band boundaries r = n D
%  (n = 1..NR-1), N points shared in proportion to r^2 (equal areal density), Fibonacci directions with the azimuth
%  shifted by half a golden angle  -> nn_grid_sph_r<R>_mid.mat
%  OCT = A: first-octant cube grid 0:H:A on x_a, y_a, z_a instead (Long et al. 2021 TMECH, Fig. 5: H = 15, A = 45,
%  4^3 = 64 points incl. the origin)  -> nn_grid_oct_h<H>_a<A>.mat, also PAR (N x 1): the position-relay parent
%  (as chain.m): node (a,b,c) = P_act / H steps back along z if c > 0, else along y if b > 0, else along x;
%  PAR = 0 for the origin (the root)
%  AXS = Nr: six-axis shells (Flux conv_design_ws axes_shells): the centre + Nr equally spaced points (step R / Nr) on
%  each of +-x_a, +-y_a, +-z_a, N = 6 Nr + 1  -> nn_grid_axs<Nr>_r<R>.mat
if ~exist('AXS','var'), AXS = []; end
if ~exist('H','var'), H = 30; end
if ~exist('OCT','var'), OCT = []; end                         % first-octant cube side A [um] ([] = off)
if ~exist('R','var'), R = 150; end
if ~exist('NR','var'), NR = []; end                           % number of shells ([] = cubic grid)
if ~exist('N','var'), N = 343; end                            % points in the ball (NR set)
if ~exist('MID','var'), MID = false; end                      % NR set: shells between the training shells
clearvars -except H R NR N MID OCT AXS;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
MATNAME = 'calib_current_maxwell_axshN55_R150_eighteen';
FN  = 'nn_grid.mat';  if H ~= 30 || R ~= 150, FN = sprintf('nn_grid_h%g_r%g.mat', H, R); end
if ~isempty(NR), FN = sprintf('nn_grid_sph_r%g.mat', R); end
if ~isempty(NR) && MID, FN = sprintf('nn_grid_sph_r%g_mid.mat', R); end
if ~isempty(OCT), FN = sprintf('nn_grid_oct_h%g_a%g.mat', H, OCT); end
if ~isempty(AXS), FN = sprintf('nn_grid_axs%d_r%g.mat', AXS, R); end
OUT = fullfile(CALROOT, 'utils', 'data', MODEL, FN);

cal   = load_flux_calib(MODEL, GEOM, MATNAME, 'current');
R_act = cal.R_act;
assert(abs(det(R_act) - 1) < 1e-6, 'nn_grid:Ract', 'R_act is not a proper rotation (det %.6g)', det(R_act));

if ~isempty(AXS)                                                % six-axis shells, N = 6 Nr + 1
    E = [eye(3); -eye(3)];  P_act = [0 0 0; kron(E, (1:AXS).' * R / AXS)];  Ns = [];  r_n = (1:AXS) * R / AXS;
elseif ~isempty(OCT)                                            % first-octant cube (Fig. 5)
    [X, Y, Z] = ndgrid((0:floor(OCT / H)) * H);
    P_act = [X(:) Y(:) Z(:)];  Ns = [];  r_n = [];
    K = round(P_act / H);  PAR = zeros(size(K, 1), 1);         % relay parent: back along z, else y, else x
    for i = 1:size(K, 1)
        k = K(i,:);  ax = find(k > 0, 1, 'last');
        if isempty(ax), continue; end                           % origin
        k(ax) = k(ax) - 1;  PAR(i) = find(all(K == k, 2));
    end
elseif isempty(NR)
    n = floor(R / H);
    [X, Y, Z] = ndgrid((-n:n) * H);
    P_act  = [X(:) Y(:) Z(:)];
    P_act  = P_act(vecnorm(P_act, 2, 2) <= R + 1e-9, :);
    Ns = [];  r_n = [];
else
    n   = 1:NR;  D = R / NR;
    r_n = D * ((n.^3 + (n-1).^3) / 2).^(1/3);                   % shell radii [um]: volume midpoints of the bands
    w   = (n.^3 - (n-1).^3) / NR^3;                             % band volume shares
    a0  = 0;                                                    % extra azimuth shift
    if MID                                                      % between the training shells
        r_n = (1:NR-1) * D;  w = r_n.^2 / sum(r_n.^2);  a0 = pi * (3 - sqrt(5)) / 2;
    end
    Ns  = round(w * N);  Ns(end) = N - sum(Ns(1:end-1));        % points per shell
    ga  = pi * (3 - sqrt(5));                                   % golden angle
    P_act = zeros(0, 3);
    for s = 1:numel(r_n)
        q  = (0:Ns(s)-1).';
        zq = 1 - (2*q + 1) / Ns(s);                             % Fibonacci lattice: equal-area bands
        aq = ga * q + s * ga + a0;                              % azimuth, shifted per shell
        u  = [sqrt(1 - zq.^2) .* cos(aq), sqrt(1 - zq.^2) .* sin(aq), zq];
        P_act = [P_act; r_n(s) * u];                            %#ok<AGROW>
    end
end
P_meas = (R_act.' * P_act.').';

save(OUT, 'P_meas', 'P_act', 'R_act', 'H', 'R', 'NR', 'N', 'Ns', 'r_n', 'MID', 'OCT', 'AXS', 'MODEL', 'MATNAME');
if ~isempty(OCT), save(OUT, 'PAR', '-append'); end
if isempty(NR), fprintf('%d points (step %g um, |p| <= %g um) -> %s\n', size(P_meas, 1), H, R, OUT);
else, fprintf('%d points (%d shells, r = %s um, per shell %s) -> %s\n', size(P_meas, 1), numel(r_n), mat2str(r_n, 4), mat2str(Ns), OUT); end
