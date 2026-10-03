%% xa_nmae.m -- NMAE of the two R <= 150 um models against the trilinear b on x_a.
%
%  Reference: xa_interp.mat -- b at 1001 equidistant actuator-frame points on the
%  x_a axis, six excitations, from the project's regular-grid trilinear.
%  Models: RBF and spherical harmonics, refitted on the same 1771 nodes with the
%  settings that passed the smoothness gate (read from rbf_v150.mat / sph_v150.mat).
%
%      NMAE = sum_j sum_i ||b_model - b_interp|| / sum_j sum_i ||b_interp|| * 100
%
%  Output: data/long2016_hexapole_halfcut/.mat/xa_nmae.mat   [ADDED 2026-09-17]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');

R0 = 150;
X  = load(fullfile(DAT,'xa_interp.mat'));          % s, P_act, B (1001 x 3 x 6)
SR = load(fullfile(DAT,'rbf_v150.mat'));
[~,iw] = min(SR.TAB.nmae);   RHO = SR.TAB.rho(iw);   LAM = SR.TAB.lam(iw);
SS = load(fullfile(DAT,'sph_v150.mat'));
r  = SS.RES([SS.RES.R] == R0);   LDEG = r.L;

[bR,~,iR] = rbf_field(R0, struct('rho',RHO,'lam',LAM));
assert(iR.Np == 1771, 'xa_nmae:nodes', 'expected 1771 nodes, got %d', iR.Np);
bS = sph_field(R0, struct('P',iR.P,'B',iR.B6), struct('L',LDEG,'Rn',R0,'quiet',true));

NQ = size(X.P_act,1);   NC = size(X.B,3);
BR = zeros(NQ,3,NC);   BS = zeros(NQ,3,NC);
for j = 1:NC
    I = zeros(6,1);  I(j) = 1;
    BR(:,:,j) = bR(X.P_act, I);
    BS(:,:,j) = bS(X.P_act, I);
end

nm  = @(D) sum(vecnorm(reshape(permute(D,[1 3 2]),[],3),2,2));
den = nm(X.B);
NMAE_R = nm(BR - X.B) / den * 100;
NMAE_S = nm(BS - X.B) / den * 100;
nmj_R = zeros(1,NC);   nmj_S = zeros(1,NC);
for j = 1:NC
    d = sum(vecnorm(X.B(:,:,j),2,2));
    nmj_R(j) = sum(vecnorm(BR(:,:,j)-X.B(:,:,j),2,2)) / d * 100;
    nmj_S(j) = sum(vecnorm(BS(:,:,j)-X.B(:,:,j),2,2)) / d * 100;
end

fprintf([newline 'RBF  (rho %g, lam %.0e) : NMAE = %.4f %%' newline], RHO, LAM, NMAE_R);
fprintf(['sph  (L = %d)              : NMAE = %.4f %%' newline newline], LDEG, NMAE_S);
fprintf(['%6s %10s %10s' newline], 'exc', 'RBF %', 'sph %');
for j = 1:NC, fprintf('%6s %10.4f %10.4f\n', sprintf('P%d',j), nmj_R(j), nmj_S(j)); end

save(fullfile(DAT,'xa_nmae.mat'), 'NMAE_R','NMAE_S','nmj_R','nmj_S','BR','BS', ...
     'RHO','LAM','LDEG','R0');
