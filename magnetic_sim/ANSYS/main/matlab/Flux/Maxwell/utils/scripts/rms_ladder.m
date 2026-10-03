%% rms_ladder.m -- calibration RMS vs number of calibration points (long2016, current)
%  Question: how does the calibration error RMS = sqrt(J/(3*6*N)) change along the
%  axes-shell ladder, for the single-parameter and the eighteen-parameter model?
%  Input : Maxwell baseline field (model_config long2016_hexapole_halfcut / tip40um, variant maxwell)
%  Output: utils/data/rms_ladder.mat  (Nr, N, J, RMS, l_hat for both models)
%  Plot  : plot/long2016_hexapole_halfcut/current/rms_ladder.m

MODEL   = 'long2016_hexapole_halfcut';
GEOM    = 'tip40um';
VARIANT = 'maxwell';
R       = 150e-6;          % sampling sphere radius [m]
l0      = 0.5e-3;          % l_hat initial guess [m]
NR      = 1:400;           % N = 6*Nr + 1 -> 7 ... 2401

here = fileparts(mfilename('fullpath'));            % .../utils/scripts
CAL  = fileparts(fileparts(here));                  % .../Flux/Maxwell
addpath(fullfile(CAL, 'function'), fullfile(CAL, 'common_path'));

cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
opt = struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, 'frame','actuator', 'raw',raw);

nq = numel(NR);
N  = nan(nq,1);   J = nan(nq,2);   l_hat = nan(nq,2);   % columns: single, eighteen
BRMS = nan(nq,1);                                    % field RMS at the calibration points [mT]
t0 = tic;
for q = 1:nq
    [P, Bs, gi] = conv_design_ws(NR(q), R, opt);
    N(q) = gi.npts_kept;
    BRMS(q) = sqrt(sum(Bs(:).^2) / (3*6*N(q)));      % same points and divisor as RMS
    for m = 1:2
        [~, l_hat(q,m), J(q,m)] = fitting(P, Bs, ad.Pc_base, l0, m == 2);
    end
    fprintf('Nr=%3d N=%3d | RMS single %.5f  eighteen %.5f mT | %.0f s\n', ...
            NR(q), N(q), sqrt(J(q,1)/(18*N(q))), sqrt(J(q,2)/(18*N(q))), toc(t0));
end
RMS = sqrt(J ./ (3*6*N));      % [mT] per residual: 3 components x 6 excitations x N points

out = fullfile(CAL, 'utils', 'data', 'rms_ladder.mat');
Nr = NR(:);
save(out, 'MODEL','GEOM','VARIANT','R','l0','Nr','N','J','RMS','l_hat','BRMS');
fprintf('saved %s\n', out);
