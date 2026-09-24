%% sph_fit.m -- solid-harmonic fit of a model's FEM field:  B = J*C,  r = J*C - B
%
%  Phi(p) = sum_{k=1..K} c_k p_k(x,y,z),   B = grad(Phi),   K = (L+1)^2 - 1
%
%  The six excitations are fitted TOGETHER: they share one design matrix J, so the
%  unknown is the K x 6 coefficient matrix C, one column per excited pole.
%
%      J   3Np x K    row (node, component c), column k = dp_k/dx_c
%      C   K   x 6    <- what is solved for; left in the workspace as C
%      B   3Np x 6    the FEM field, stacked in the same row order
%      r = J*C - B    the residual; the fit minimises its 2-norm, and the number
%                     reported below is rms(r) = sqrt(mean(r.^2)) in mT
%
%  The nodes inside R are split TRAIN / TEST (default 80 / 20). Only the train
%  part enters the fit; the test part is scored afterwards with the same r, which
%  is what tells overfitting from a genuine improvement when L is raised.
%
%  Run it and inspect C (and info.lab, which names the (l, m, cos/sin) of each row).

clear;  clc;

%% ---- per-run knobs ---------------------------------------------------------
MODEL   = 'long2016_hexapole_halfcut';
GEOM    = 'tip40um';
VARIANT = 'maxwell';        % No gap baseline
R       = 150;              % fitted region, |p| <= R [um]
L       = 10;               % maximum degree -> K = (L+1)^2 - 1
FTRAIN  = 0.8;              % fraction of the nodes used for the fit
SEED    = 0;                % rng seed for the split, so the run repeats

%% ---- paths -----------------------------------------------------------------
here = fileparts(mfilename('fullpath'));            % .../Force/Maxwell/utils
FMX  = fileparts(here);                             % .../Force/Maxwell
MAIN = fileparts(fileparts(fileparts(FMX)));        % .../ANSYS/main
FLX  = fullfile(MAIN, 'matlab', 'Flux', 'Maxwell');
addpath(fullfile(FMX,'utils'), fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));

%% ---- FEM field: actuator frame, um, mT, all-source --------------------------
cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);    % .fld, Maxwell frame, T
ad  = build_actuator_data(raw, cfg);                % actuator frame, mT, iron filtered
P   = ad.Pa * 1e6;                                  % [m] -> [um]
Bf  = ad.Ba;                                        % Np x 3 x 6 [mT]

in  = vecnorm(P, 2, 2) <= R;
P   = P(in,:);   Bf = Bf(in,:,:);   Np = size(P,1);

%% ---- 80 / 20 split ---------------------------------------------------------
rng(SEED);
ix   = randperm(Np);
ntr  = round(FTRAIN * Np);
itr  = ix(1:ntr);         ite = ix(ntr+1:end);
Ptr  = P(itr,:);   Btr = Bf(itr,:,:);
Pte  = P(ite,:);   Bte = Bf(ite,:,:);
fprintf(['nodes |p| <= %g um : %d   ->  train %d (%.0f %%) / test %d' newline], ...
        R, Np, numel(itr), 100*numel(itr)/Np, numel(ite));

%% ---- fit (train only) ------------------------------------------------------
[b, grad, info] = sph_field(R, struct('P',Ptr,'B',Btr), struct('L',L,'Rn',R));
C = info.C;                                          % K x 6 -- the answer

%% ---- residual r = J*C - B on both parts ------------------------------------
NC   = size(C,2);
Bp   = zeros(size(Bte));
for k = 1:NC
    I = zeros(NC,1);  I(k) = 1;
    Bp(:,:,k) = b(Pte, I);
end
r_te = reshape(permute(Bp - Bte, [1 2 3]), [], NC);  % same elements as J*C - B
rms_te = sqrt(mean(r_te.^2, 1));

fprintf([newline '%-8s %12s %12s' newline], 'coil', 'train rms', 'test rms');
for k = 1:NC
    fprintf('%-8d %12.5g %12.5g\n', k, info.rms_exc(k), rms_te(k));
end
fprintf('%-8s %12.5g %12.5g   [mT]\n', 'all', info.rms, sqrt(mean(r_te.^2,'all')));
fprintf('K = %d, cond(J) = %.2e, NMAE(train) = %.4f %%\n', info.K, info.cnd, info.NMAE_all);

%% ---- the coefficient matrix ------------------------------------------------
lab = info.lab;                                      % [l m t], t = 1 cos, 0 sin
nm  = arrayfun(@(q) sprintf('l%d_m%d%s', lab(q,1), lab(q,2), char('s'*(lab(q,3)==0) + 'c'*(lab(q,3)==1))), ...
               (1:info.K).', 'UniformOutput', false);
Ctab = array2table(C, 'VariableNames', compose("coil%d", 1:NC), 'RowNames', nm);
disp(Ctab);                                          % C is also plain in the workspace
