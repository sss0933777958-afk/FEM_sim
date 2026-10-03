%% s150_200.m -- the SAME question as x200.m, asked of the spherical harmonics:
%%                what does widening the data range 150 -> 200 um do to the gradient?
%
%  Both models are solid harmonics.  Each is fitted on every node inside its own
%  radius, at the degree that passes the smoothness gate there (read from
%  sph_v150.mat / sph_v200.mat), and both are evaluated on the SAME 1001 equidistant
%  points per actuator axis over |s| <= 150 um.
%
%      mean_error = mean |grad_200 - grad_150|,  the same definition used for
%                   |RBF - sph| in x150.m and x200.m, so the three are comparable.
%
%  SIGNS: z_a is the P5 direction and the excited pole is P6, at negative z_a, so
%  that axis is negated -- identical to the other scripts.
%
%  Output: data/long2016_hexapole_halfcut/.mat/s150_200.mat   [ADDED 2026-09-17]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');

RL   = [150 200];
SPAN = 150;   NQ = 1001;
POLE = [1 3 6];   DSGN = [1 1 -1];   ANM = {'x_a','y_a','z_a'};

s  = linspace(-SPAN, SPAN, NQ).';
G  = cell(1,2);   LD = zeros(1,2);   NP = zeros(1,2);
for k = 1:2
    R0 = RL(k);
    S  = load(fullfile(DAT, sprintf('sph_v%d.mat', R0)));
    r  = S.RES([S.RES.R] == R0);
    LD(k) = r.L;
    fprintf(['sph @ R = %3d : %d of %d degrees pass, L = %d (K = %d), NMAE(R<=150) = %.4f %%' newline], ...
            R0, r.npass, r.Ltried, r.L, r.K, r.nmae);
    % the node set: taken through rbf_field so it is the same selection the other
    % scripts used (no iron inside 500 um, so this is simply |p| <= R0)
    [~,~,iN] = rbf_field(R0, struct('rho',400,'lam',1e-8));
    NP(k) = iN.Np;
    [~, gS] = sph_field(R0, struct('P',iN.P,'B',iN.B6), ...
                        struct('L',LD(k), 'Rn',R0, 'quiet',true, 'condmax',1500));
    % condmax: cond() is an SVD of a 3N x K matrix and is only a diagnostic; at the
    % degrees this sweep reaches (K ~ 4000) it costs minutes and nothing reads it.
    Gk = zeros(NQ,3);
    for a = 1:3
        Pq = zeros(NQ,3);   Pq(:,a) = s;
        I  = zeros(6,1);    I(POLE(a)) = 1;
        Gq = gS(Pq, I, 'du');   Gk(:,a) = DSGN(a) * Gq(:,a);
    end
    G{k} = Gk;
end
fprintf(['nodes: %d at R = 150, %d at R = 200' newline], NP(1), NP(2));

D = G{2} - G{1};                                   % 200 minus 150
fprintf([newline '%6s %16s %16s %16s %14s' newline], ...
        'axis','mean(200-150)','mean|200-150|','max|200-150|','mean|.|/mean|grad|');
mabs = zeros(1,3);   msig = zeros(1,3);   mrel = zeros(1,3);
for a = 1:3
    msig(a) = mean(D(:,a));
    mabs(a) = mean(abs(D(:,a)));
    mrel(a) = mabs(a) / mean(abs(G{1}(:,a))) * 100;
    fprintf('%6s %16.3e %16.3e %16.3e %13.3f %%\n', ANM{a}, msig(a), mabs(a), ...
            max(abs(D(:,a))), mrel(a));
end
mean_error = mean(abs(D(:)));
fprintf([newline 'mean_error (sph R200 vs sph R150) = %.4e mT^2/um  (%.3f %% of the mean ' ...
         '|gradient|)' newline], mean_error, mean_error/mean(abs(G{1}(:)))*100);

save(fullfile(DAT,'s150_200.mat'), 's','G','D','mean_error','mabs','msig','mrel', ...
     'LD','NP','RL','SPAN','NQ','POLE','DSGN','ANM');
