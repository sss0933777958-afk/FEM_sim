%% xa_smooth_check.m -- do the four all-node RBF fits actually have a smooth
%%                       gradient on the three actuator axes?
%
%  This is a CHECK, not a search.  Each range is fitted on ALL its nodes with the
%  kernel the earlier scan picked for it, and then the resulting gradient is tested
%  for smoothness -- d3(b.b)/ds^3 keeping one sign along each actuator axis:
%
%      x_a <- P1 excited      y_a <- P3      z_a <- P6
%
%  Two test spans are reported for every fit, because they answer different things:
%      |s| <= 150   the working region, which is all the figure shows
%      |s| <= R     the range's own extent, which is what the original scan used
%  A fit can pass the first and fail the second; that difference is the whole reason
%  the large-R kernels were forced to be so heavily damped.
%
%  Kernels (from range_scan.mat, selected there on an 80/20 split):
%      R=150 rho 460 lam 1e-7     R=350 rho 500 lam 1e-3
%      R=250 rho 280 lam 1e-5     R=450 rho 610 lam 1e-2
%
%  Prints, per range and per axis, the number of sign changes in d3 over each span.
%  0 = smooth.  [ADDED 2026-09-17]

clear;  clc;

RS   = [150 250 350 450];
RHO  = [460 280 500 610];
LAM  = [1e-7 1e-5 1e-3 1e-2];
POLE = [1 3 6];   SG = [1 1 -1];   NQ = 301;   AXN = {'x_a','y_a','z_a'};

here=fileparts(mfilename('fullpath'));  FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX)));  FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));

fprintf(['%5s %6s %8s %8s | %-22s | %-22s' newline], ...
        'R','rho','lam','nodes','sign changes |s|<=150','sign changes |s|<=R');
fprintf([repmat('-',1,84) newline]);
for q = 1:numel(RS)
    R0 = RS(q);
    [~, gq, iq] = rbf_field(R0, struct('rho',RHO(q),'lam',LAM(q)));
    nch = zeros(2,3);
    for w = 1:2
        if w == 1, SPAN = 150; else, SPAN = R0; end
        sA = linspace(-SPAN, SPAN, NQ).';
        for a = 1:3
            Pv = zeros(NQ,3);  Pv(:,a) = sA;
            I  = zeros(6,1);   I(POLE(a)) = 1;
            d3 = gq(Pv, I, 'd3');
            v  = SG(a)*d3(:,a);   g = sign(v);   nz = find(g~=0);
            nch(w,a) = sum(diff(g(nz)) ~= 0);
        end
    end
    v1 = 'PASS'; if any(nch(1,:)), v1 = 'FAIL'; end
    v2 = 'PASS'; if any(nch(2,:)), v2 = 'FAIL'; end
    fprintf(['%5d %6g %8.0e %8d | [%d %d %d] %-4s            | [%d %d %d] %-4s' newline], ...
            R0, RHO(q), LAM(q), size(iq.P,1), nch(1,:), v1, nch(2,:), v2);
    clear gq iq
end
fprintf([newline 'axes: %s <- P1, %s <- P3, %s <- P6 (z_a flipped, as everywhere)' newline], AXN{:});
