%% s_L2.m -- the degree chosen as "the first that passes, excluding L = 1".
%
%  L = 1 is excluded because it cannot fail: a linear potential gives a constant b,
%  so the gradient is identically zero and its third derivative never changes sign.
%  The next degree up is therefore the first that carries any information, and at
%  both ranges the gate passes from there on (R = 150: L = 1..9; R = 200: L = 1..35),
%  so the rule returns L = 2 for both.  [user's call, 2026-09-18]
%
%  Reports, for R = 150 and R = 200: the working-region NMAE of each fit, the
%  gradient on the three actuator axes at the same 1001 equidistant points, and the
%  mean |R200 - R150| in the same convention as every other mean_error here.
%
%  Output: data/long2016_hexapole_halfcut/.mat/s_L2.mat   [ADDED 2026-09-18]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));  FMX = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');

RL = [150 200];   LSEL = 2;
SPAN = 150;   NQ = 1001;
POLE = [1 3 6];   SG = [1 1 -1];   DSGN = [1 1 -1];   ANM = {'x_a','y_a','z_a'};

s = linspace(-SPAN,SPAN,NQ).';
AXQ = cell(1,3);
for a = 1:3, Pv = zeros(NQ,3);  Pv(:,a) = s;  AXQ{a} = Pv;  end

% the R <= 150 nodes are the NMAE set for both fits, as everywhere else
[~,~,iW] = rbf_field(150, struct('rho',400,'lam',1e-8));
Pw = iW.P;   Bw = iW.B6;   Nw = iW.Np;
nBw = sum(vecnorm(reshape(permute(Bw,[1 3 2]),[],3),2,2));

G = cell(1,2);   NM = zeros(1,2);   NP = zeros(1,2);   RIP = zeros(2,3);
for k = 1:2
    R0 = RL(k);
    [~,~,iN] = rbf_field(R0, struct('rho',400,'lam',1e-8));
    NP(k) = iN.Np;
    [bq, gq, iq] = sph_field(R0, struct('P',iN.P,'B',iN.B6), ...
                             struct('L',LSEL,'Rn',R0,'quiet',true,'condmax',1500));
    Bp = zeros(Nw,3,6);
    for j = 1:6, I = zeros(6,1); I(j) = 1;  Bp(:,:,j) = bq(Pw, I);  end
    NM(k) = sum(vecnorm(reshape(permute(Bp-Bw,[1 3 2]),[],3),2,2)) / nBw * 100;
    Gk = zeros(NQ,3);
    for a = 1:3
        I = zeros(6,1);  I(POLE(a)) = 1;
        Q  = gq(AXQ{a}, I, 'du');   Gk(:,a) = DSGN(a)*Q(:,a);
        d3 = gq(AXQ{a}, I, 'd3');   RIP(k,a) = nsc_(SG(a)*d3(:,a));
    end
    G{k} = Gk;
    fprintf(['sph @ R = %3d, L = %d (K = %d), %d nodes : NMAE(R<=150) = %.4f %%, ' ...
             'ripples %s' newline], R0, LSEL, iq.K, NP(k), NM(k), mat2str(RIP(k,:)));
end

D = G{2} - G{1};
fprintf([newline '%6s %16s %16s %16s %14s' newline], ...
        'axis','mean(200-150)','mean|200-150|','max|200-150|','mean|.|/mean|grad|');
mabs = zeros(1,3);  msig = zeros(1,3);  mrel = zeros(1,3);
for a = 1:3
    msig(a)=mean(D(:,a));  mabs(a)=mean(abs(D(:,a)));
    mrel(a)=mabs(a)/mean(abs(G{1}(:,a)))*100;
    fprintf('%6s %16.3e %16.3e %16.3e %13.3f %%\n', ANM{a}, msig(a), mabs(a), ...
            max(abs(D(:,a))), mrel(a));
end
mean_error = mean(abs(D(:)));
fprintf([newline 'mean_error (sph R200 L2 vs sph R150 L2) = %.4e mT^2/um  (%.3f %%)' newline], ...
        mean_error, mean_error/mean(abs(G{1}(:)))*100);

save(fullfile(DAT,'s_L2.mat'), 's','G','D','mean_error','mabs','msig','mrel', ...
     'NM','NP','RIP','LSEL','RL','SPAN','NQ','POLE','DSGN','ANM');

function n = nsc_(v)
    g = sign(v);  nz = find(g~=0);  n = sum(diff(g(nz)) ~= 0);
end
