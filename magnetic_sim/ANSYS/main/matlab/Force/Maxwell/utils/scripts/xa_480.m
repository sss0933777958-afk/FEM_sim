%% xa_480.m -- RBF at R = 480, gradient at x_a = +150 vs the eighteen-parameter model
%  R = 480 was not in the range scan (its ladder goes ...450, 500), so the kernel is
%  INTERPOLATED between the two neighbours, which share lam = 1e-2:
%      R=450 -> rho 610      R=500 -> rho 740      => R=480 -> rho 690
%  That is not this range's own optimum; it is the best available stand-in.
%  [ADDED 2026-09-17, diagnostic]
clear; clc;
here=fileparts(mfilename('fullpath')); FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX))); FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
Pq = [linspace(-150,150,601).' zeros(601,2)];
I  = zeros(6,1); I(1) = 1;
t0 = tic;
[~,gq,iq] = rbf_field(480, struct('rho',690,'lam',1e-2));
G = gq(Pq, I, 'du');  g = G(:,1);
fprintf(['RBF R=480 (rho 690, lam 1e-2, %d nodes, %.0f s)' newline], size(iq.P,1), toc(t0));
fprintf(['  x_a=+150: %.4f   x_a=-150: %.4f   mT^2/um' newline], g(end), g(1));
ref = 1.0040;  refm = 0.1874;                      % eighteen-parameter model
fprintf(['  rel. err vs eighteen parameters:  +150: %+.2f %%   -150: %+.2f %%' newline], ...
        (g(end)-ref)/ref*100, (g(1)-refm)/refm*100);
