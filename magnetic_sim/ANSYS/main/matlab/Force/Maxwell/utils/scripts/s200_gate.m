%% s200_gate.m -- the gate verdict for every admissible degree at R = 200 um.
%
%  L = 1..35 and the failure at 36 are already known from s200_first.m; this fills in
%  36..63 (63 is the ceiling: K = 4095 against 4164 nodes, L = 64 would exceed it) and
%  merges the two into one table.  The per-degree log the original sweep kept was lost
%  to a save() that omitted LTAB, since fixed in sph_strat.m.
%
%  Output: data/long2016_hexapole_halfcut/.mat/s200_gate.mat   [ADDED 2026-09-18]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));  FMX = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');

SPAN = 150;   NQ = 1001;   LMAXA = 63;
POLE = [1 3 6];   SG = [1 1 -1];
s = linspace(-SPAN,SPAN,NQ).';
AXQ = cell(1,3);
for a = 1:3, Pv = zeros(NQ,3);  Pv(:,a) = s;  AXQ{a} = Pv;  end

[~,~,i200] = rbf_field(200, struct('rho',400,'lam',1e-8));
S = load(fullfile(DAT,'s200_first.mat'));      % L = 1..36 already measured
PASS = nan(1,LMAXA);   PASS(1:35) = 1;   PASS(36) = 0;
KK   = arrayfun(@(L)(L+1)^2-1, 1:LMAXA);
for Lq = 37:LMAXA
    t0 = tic;
    [~, gq] = sph_field(200, struct('P',i200.P,'B',i200.B6), ...
                        struct('L',Lq,'Rn',200,'quiet',true,'condmax',1500));
    ok = true;
    for a = 1:3
        I = zeros(6,1);  I(POLE(a)) = 1;
        d3 = gq(AXQ{a}, I, 'd3');
        if nsc_(SG(a)*d3(:,a)) > 0, ok = false;  break, end
    end
    PASS(Lq) = ok;
    vs = {'FAIL','pass'};
    fprintf(['  L = %2d (K = %4d): %s   [%.0f s]' newline], Lq, KK(Lq), vs{ok+1}, toc(t0));
end
fprintf([newline 'passes: %d of %d' newline], nnz(PASS==1), LMAXA);
fprintf(['failing degrees: %s' newline], mat2str(find(PASS==0)));
save(fullfile(DAT,'s200_gate.mat'), 'PASS','KK','LMAXA','NQ','SPAN','POLE','SG');

function n = nsc_(v)
    g = sign(v);  nz = find(g~=0);  n = sum(diff(g(nz)) ~= 0);
end
