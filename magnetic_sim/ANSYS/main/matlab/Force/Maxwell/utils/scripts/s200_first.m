%% s200_first.m -- the FIRST converged degree at R = 200, instead of the highest one.
%
%  The uncapped "highest still-smooth degree, patience 10" rule ran to the
%  mathematical ceiling at R = 200 (L = 63, K = 4095 against 4164 nodes -- 69 degrees
%  of freedom left, effectively an interpolant).  This script takes the other reading
%  of the gate instead: sweep L upward and stop at the FIRST failure, keeping the
%  degree just below it.  At R = 150 that rule returns L = 9, the same answer the
%  patience rule gave there, so the two agree wherever the sweep converges at all.
%
%  Then the same comparison as s150_200.m: the gradient on the three actuator axes at
%  1001 equidistant points, sph(R=200, L*) against sph(R=150, L=9).
%
%  Output: data/long2016_hexapole_halfcut/.mat/s200_first.mat   [ADDED 2026-09-17]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');

SPAN = 150;   NQ = 1001;   LCAP = 40;
POLE = [1 3 6];   SG = [1 1 -1];   DSGN = [1 1 -1];   ANM = {'x_a','y_a','z_a'};

s   = linspace(-SPAN, SPAN, NQ).';
AXQ = cell(1,3);
for a = 1:3, Pv = zeros(NQ,3);  Pv(:,a) = s;  AXQ{a} = Pv;  end

%% ---- sweep R = 200 upward, stop at the first gate failure --------------------
[~,~,i200] = rbf_field(200, struct('rho',400,'lam',1e-8));
fprintf(['R = 200: %d nodes' newline], i200.Np);
Lstar = NaN;   PASS = false(1,LCAP);
for Lq = 1:LCAP
    [~, gq, iq] = sph_field(200, struct('P',i200.P,'B',i200.B6), ...
                            struct('L',Lq,'Rn',200,'quiet',true,'condmax',1500));
    ok = true;
    for a = 1:3
        I = zeros(6,1);  I(POLE(a)) = 1;
        d3 = gq(AXQ{a}, I, 'd3');
        if nsc_(SG(a)*d3(:,a)) > 0, ok = false;  break, end
    end
    PASS(Lq) = ok;
    vs = {'FAIL','pass'};
    fprintf(['  L = %2d (K = %4d): %s' newline], Lq, iq.K, vs{ok+1});
    if ~ok, Lstar = Lq - 1;  break, end
end
assert(~isnan(Lstar), 's200_first:noFail', 'no degree failed up to L = %d', LCAP);
fprintf([newline 'first failure at L = %d  ->  first converged degree L* = %d' newline], ...
        Lstar+1, Lstar);

%% ---- the two models and their gradients --------------------------------------
[~,~,i150] = rbf_field(150, struct('rho',400,'lam',1e-8));
S150 = load(fullfile(DAT,'sph_v150.mat'));  r150 = S150.RES([S150.RES.R]==150);
[~, g150] = sph_field(150, struct('P',i150.P,'B',i150.B6), ...
                      struct('L',r150.L,'Rn',150,'quiet',true,'condmax',1500));
[~, g200] = sph_field(200, struct('P',i200.P,'B',i200.B6), ...
                      struct('L',Lstar,'Rn',200,'quiet',true,'condmax',1500));
G150 = zeros(NQ,3);   G200 = zeros(NQ,3);
for a = 1:3
    I = zeros(6,1);  I(POLE(a)) = 1;
    Q = g150(AXQ{a}, I, 'du');   G150(:,a) = DSGN(a)*Q(:,a);
    Q = g200(AXQ{a}, I, 'du');   G200(:,a) = DSGN(a)*Q(:,a);
end
D = G200 - G150;

fprintf([newline 'sph R=150 L=%d   vs   sph R=200 L=%d' newline], r150.L, Lstar);
fprintf(['%6s %16s %16s %16s %14s' newline], ...
        'axis','mean(200-150)','mean|200-150|','max|200-150|','mean|.|/mean|grad|');
mabs = zeros(1,3);  msig = zeros(1,3);  mrel = zeros(1,3);
for a = 1:3
    msig(a)=mean(D(:,a));  mabs(a)=mean(abs(D(:,a)));
    mrel(a)=mabs(a)/mean(abs(G150(:,a)))*100;
    fprintf('%6s %16.3e %16.3e %16.3e %13.3f %%\n', ANM{a}, msig(a), mabs(a), ...
            max(abs(D(:,a))), mrel(a));
end
mean_error = mean(abs(D(:)));
fprintf([newline 'mean_error = %.4e mT^2/um  (%.3f %% of the mean |gradient|)' newline], ...
        mean_error, mean_error/mean(abs(G150(:)))*100);

save(fullfile(DAT,'s200_first.mat'), 's','G150','G200','D','mean_error','mabs','msig', ...
     'mrel','Lstar','PASS','ANM','POLE','DSGN','NQ','SPAN');

function n = nsc_(v)
    g = sign(v);  nz = find(g~=0);  n = sum(diff(g(nz)) ~= 0);
end
