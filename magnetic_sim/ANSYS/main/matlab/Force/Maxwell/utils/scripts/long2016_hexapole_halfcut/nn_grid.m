%% nn_grid.m -- training positions of the NN inverse model: cubic grid in the ACTUATOR frame.
%  Grid step H [um], axes along x_a, y_a, z_a (the frame the network input p/R and rho_solve use),
%  centred on the tip-sphere centre, keeping the nodes with |p| <= R.  The same points are also
%  stored in the measurement frame (P_meas = R_act' * P_act).
%  Output: utils/data/long2016_hexapole_halfcut/nn_grid.mat
%    P_act   N x 3 [um]   actuator frame, origin = sphere centre (the grid)
%    P_meas  N x 3 [um]   measurement frame
%    R_act   3 x 3        rows = P1, P3, P5 tip unit vectors (actuator-frame.md)
%    H, R               grid step and radius [um]
if ~exist('H','var'), H = 30; end
if ~exist('R','var'), R = 150; end
clearvars -except H R;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'));
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
MATNAME = 'calib_current_maxwell_axshN55_R150_eighteen';
OUT = fullfile(CALROOT, 'utils', 'data', MODEL, 'nn_grid.mat');

cal   = load_flux_calib(MODEL, GEOM, MATNAME, 'current');
R_act = cal.R_act;
assert(abs(det(R_act) - 1) < 1e-6, 'nn_grid:Ract', 'R_act is not a proper rotation (det %.6g)', det(R_act));

n = floor(R / H);
[X, Y, Z] = ndgrid((-n:n) * H);
P_act  = [X(:) Y(:) Z(:)];
P_act  = P_act(vecnorm(P_act, 2, 2) <= R + 1e-9, :);
P_meas = (R_act.' * P_act.').';

save(OUT, 'P_meas', 'P_act', 'R_act', 'H', 'R', 'MODEL', 'MATNAME');
fprintf('%d points (step %g um, |p| <= %g um) -> %s\n', size(P_meas, 1), H, R, OUT);
