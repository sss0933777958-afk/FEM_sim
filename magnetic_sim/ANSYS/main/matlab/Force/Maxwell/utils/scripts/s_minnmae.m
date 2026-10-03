%% s_minnmae.m -- among the degrees that PASS the smoothness gate, the one with the
%%                 lowest working-region NMAE.
%
%  Sweeps L upward at each data range, records the gate verdict and the NMAE over the
%  1771 nodes inside R <= 150 um (six excitations, vector-norm convention, the same
%  measure used everywhere else), and picks the lowest-NMAE passer.
%
%  RANGE OF THE SWEEP.  L runs to LTOP, which is deliberately modest: at R = 200 the
%  gate passes L = 1..35 solidly and then only sporadically (39, 41, 43, 50 pass; 36,
%  37, 38, 40, 42, 44..49 fail), and each degree past 50 costs minutes.  So this
%  covers the solid block and the answer it gives is the best WITHIN IT -- a sporadic
%  high degree could in principle score lower and is not ruled out here.
%
%  Output: data/long2016_hexapole_halfcut/.mat/s_minnmae.mat   [ADDED 2026-09-18]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));  FMX = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');

RL = [150 200];   LTOP = [19 35];
SPAN = 150;   NQ = 1001;
POLE = [1 3 6];   SG = [1 1 -1];

s = linspace(-SPAN,SPAN,NQ).';
AXQ = cell(1,3);
for a = 1:3, Pv = zeros(NQ,3);  Pv(:,a) = s;  AXQ{a} = Pv;  end

[~,~,iW] = rbf_field(150, struct('rho',400,'lam',1e-8));
Pw = iW.P;   Bw = iW.B6;   Nw = iW.Np;
nBw = sum(vecnorm(reshape(permute(Bw,[1 3 2]),[],3),2,2));

TAB = cell(1,2);   BEST = zeros(1,2);
for k = 1:2
    R0 = RL(k);
    [~,~,iN] = rbf_field(R0, struct('rho',400,'lam',1e-8));
    T = nan(LTOP(k), 4);                       % L, K, pass, NMAE
    for Lq = 1:LTOP(k)
        [bq, gq, iq] = sph_field(R0, struct('P',iN.P,'B',iN.B6), ...
                                 struct('L',Lq,'Rn',R0,'quiet',true,'condmax',1500));
        ok = true;
        for a = 1:3
            I = zeros(6,1);  I(POLE(a)) = 1;
            d3 = gq(AXQ{a}, I, 'd3');
            if nsc_(SG(a)*d3(:,a)) > 0, ok = false;  break, end
        end
        Bp = zeros(Nw,3,6);
        for j = 1:6, I = zeros(6,1); I(j) = 1;  Bp(:,:,j) = bq(Pw, I);  end
        nm = sum(vecnorm(reshape(permute(Bp-Bw,[1 3 2]),[],3),2,2)) / nBw * 100;
        T(Lq,:) = [Lq, iq.K, ok, nm];
    end
    TAB{k} = T;
    p = T(T(:,3)==1, :);
    [~,ib] = min(p(:,4));   BEST(k) = p(ib,1);
    fprintf([newline 'R = %d : %d of %d degrees pass;  lowest-NMAE passer = L %d ' ...
             '(K = %d), NMAE = %.4f %%' newline], R0, nnz(T(:,3)), LTOP(k), ...
            p(ib,1), p(ib,2), p(ib,4));
    fprintf(['%4s %6s %6s %10s' newline], 'L','K','gate','NMAE %');
    vs = {'FAIL','pass'};
    for i = 1:size(T,1)
        fprintf('%4d %6d %6s %10.4f\n', T(i,1), T(i,2), vs{T(i,3)+1}, T(i,4));
    end
end
save(fullfile(DAT,'s_minnmae.mat'), 'TAB','BEST','RL','LTOP','NQ','SPAN','POLE','SG');

function n = nsc_(v)
    g = sign(v);  nz = find(g~=0);  n = sum(diff(g(nz)) ~= 0);
end
