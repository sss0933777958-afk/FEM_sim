function S = axsh_sweep(force, MODEL, GEOM, VARIANT, NMAX, Rum, NRV)
%AXSH_SWEEP  Calibrate along the actuator-axis shell ladder and cache the trend.
% =========================================================================
%   Engine only -- draws nothing. The three plot scripts read its .mat:
%       axsh_ell.m    l_hat        vs point count
%       axsh_gain.m   g_I_hat      vs point count
%       axsh_kfro.m   ||K - K_full||_F / ||K_full||_F   vs point count
%
%   Ladder (sample_axes_shells): the six actuator-axis directions on Nr
%   equally spaced shells r_k = (R/Nr)*k, plus the centre -> 6*Nr+1 points.
%   Nr runs 1 ... floor((NMAX-1)/6), so NMAX = 500 gives Nr <= 83 (499 points).
%
%   Reference: every real .fld lattice point inside R is calibrated too, for
%   both models. That is the "truth" the trends are compared against --
%   plotted as a dashed line by axsh_ell / axsh_gain, and used as the
%   denominator of the Frobenius error by axsh_kfro.
%
%   Cached to  figures/paper_fig_plot/data/axsh_R<Rum>[_<model>].mat
%   Fields  N ell1 ell2 gI1 gI2 K1 K2 ki1 ki2 fro1 fro2
%           ref1 ref2 Kref1 Kref2 nref R MODEL GEOM VARIANT
%     ...1 = single parameter, ...2 = eighteen parameters.
%
%   Usage
%     axsh_sweep;                                            % long2016, R=150, N<=500
%     axsh_sweep(true, 'zhi_peng','R500','maxwell_split');   % zhi-Peng
% =========================================================================
    clc;
    if nargin < 1 || isempty(force),   force   = false;                       end
    if nargin < 2 || isempty(MODEL),   MODEL   = 'long2016_hexapole_halfcut'; end
    if nargin < 3 || isempty(GEOM),    GEOM    = 'tip40um';                   end
    if nargin < 4,                     VARIANT = '';                          end
    if nargin < 5 || isempty(NMAX),    NMAX    = 500;                         end   % point cap
    if nargin < 6 || isempty(Rum),     Rum     = 150;                         end   % sampling radius [um]
    % [ADDED 2026-08-27] NRV: explicit list of Nr values. Cost is O(sum N) = O(Nr^2),
    %   so a uniform 1:Nr_max is unaffordable once the tail reaches thousands of points.
    %   Default = dense 1..floor((NMAX-1)/6); pass a custom vector to add a SPARSE tail
    %   that pins the steady-state error, e.g.
    %       [1:83, 90:10:200, 220:20:400, 450:50:800]   -> N up to 4801
    if nargin < 7 || isempty(NRV), NRV = 1:floor((NMAX - 1)/6); end
    NRV = unique(NRV(:).');

    % [MODIFIED 2026-08-30] 本檔由 temp_code/scripts/ 搬進 paper_fig_plot/plot/（產論文圖的腳本
    %   一律住 plot/，見 figure-output.md）。路徑改成相對自身推導，不再寫死絕對路徑：
    %   plot/ -> paper_fig_plot/ -> figures/ -> main/
    MAIN = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    CAL  = fullfile(MAIN,'matlab','Flux','Maxwell');
    addpath(fullfile(CAL,'function'), fullfile(CAL,'common_path'), fullfile(CAL,'utils'));
    % [REMOVED 2026-08-30] 原本 addpath temp_code/scripts 是為了 sample_axes_shells；
    %   它已搬進 matlab/Flux/Maxwell/function/（上一行的 CAL addpath 已涵蓋）-> 該行是死碼。

    msuf = '';
    if ~strcmp(MODEL,'long2016_hexapole_halfcut')
        msuf = ['_' regexprep(MODEL,'_.*$','')];
    end
    cachef = fullfile(MAIN,'figures','paper_fig_plot','data', ...
                      sprintf('axsh_R%d%s.mat', Rum, msuf));
    if exist(cachef,'file') && ~force
        S = load(cachef);   fprintf('由快取載入 %s\n', cachef);   return
    end

    R   = Rum*1e-6;
    l0  = 0.5e-3;                                   % l_hat initial guess [m]
    cfg = model_config(MODEL, GEOM);
    if isempty(VARIANT), VARIANT = cfg.default_variant; end
    F = zeros(6, cfg.N_I);
    for j = 1:cfg.N_I, F(cfg.apdl_to_paper_idx(j), j) = 1; end
    o = struct('frame','actuator','model',MODEL,'geom',GEOM,'variant',VARIANT,'quiet',true);

    % ---- full-grid reference -------------------------------------------------
    raw = extract_maxwell_data(cfg, 'all', VARIANT);
    ad  = build_actuator_data(raw, cfg);
    [P0, B0, nref] = cfg.select_ball(ad, R);
    [er,lr] = fitting(P0,B0,cfg.Pc_base,l0,false);
    [Kr,gr] = solve_current(lr,er,cfg.Pc_base,P0,B0,F,[]);
    S.ref1 = [lr*1e6, gr];   S.Kref1 = Kr;
    [er,lr] = fitting(P0,B0,cfg.Pc_base,l0,true);
    [Kr,gr] = solve_current(lr,er,cfg.Pc_base,P0,B0,F,[]);
    S.ref2 = [lr*1e6, gr];   S.Kref2 = Kr;
    S.nref = nref;
    fprintf('全格點基準（%d 點）：1p l=%.2f um g=%.4f mT/A | 18p l=%.2f um g=%.4f mT/A\n', ...
            nref, S.ref1, S.ref2);

    % ---- ladder --------------------------------------------------------------
    NRMAX = numel(NRV);
    [N, ell1, ell2, gI1, gI2, fro1, fro2] = deal(nan(1,NRMAX));
    K1 = nan(6,6,NRMAX);   K2 = nan(6,6,NRMAX);
    [ki1, ki2] = deal(false(1,NRMAX));
    for q = 1:NRMAX
        Pq = conv_design_ws(NRV(q), R, struct('points_only',true,'R_act',cfg.R_act,'quiet',true));
        evalc('[P, Bs] = conv_design_ws([], R, setfield(o,''query'',Pq));');
        N(q) = size(P,1);
        [ell1(q), K1(:,:,q), ki1(q), gI1(q)] = fit_one(P, Bs, cfg.Pc_base, l0, false, F);
        [ell2(q), K2(:,:,q), ki2(q), gI2(q)] = fit_one(P, Bs, cfg.Pc_base, l0, true,  F);
        fro1(q) = frob_rel(K1(:,:,q), S.Kref1);
        fro2(q) = frob_rel(K2(:,:,q), S.Kref2);
        fprintf(['  Nr=%3d  N=%4d | 1p l=%7.2f g=%7.4f dK=%6.3f%% ki=%d' ...
                 ' | 18p l=%7.2f g=%7.4f dK=%6.3f%% ki=%d\n'], ...
                NRV(q), N(q), ell1(q), gI1(q), fro1(q), ki1(q), ell2(q), gI2(q), fro2(q), ki2(q));
    end

    S.N = N;   S.ell1 = ell1;   S.ell2 = ell2;   S.gI1 = gI1;   S.gI2 = gI2;
    S.K1 = K1; S.K2 = K2;       S.ki1 = ki1;     S.ki2 = ki2;
    S.fro1 = fro1;   S.fro2 = fro2;
    S.Nr = NRV;
    S.R = R;   S.Rum = Rum;   S.NMAX = NMAX;
    S.MODEL = MODEL;   S.GEOM = GEOM;   S.VARIANT = VARIANT;
    save(cachef, '-struct', 'S');
    fprintf('已存 %s\n', cachef);
end

% ============================================================================
function [ell_um, K, kiok, gI] = fit_one(P, Bstack, Pc_base, l0, USE_BIAS, F)
% One design -> l_hat [um], K_I_bar, physicality flag, g_I_hat [mT/A].
    try
        evalc(['[e, l] = fitting(P, Bstack, Pc_base, l0, USE_BIAS);' ...
               '[K, gI] = solve_current(l, e, Pc_base, P, Bstack, F, []);']);
        ell_um = l*1e6;
        [~, am] = max(abs(K), [], 2);   od = K(~eye(6));
        kiok = all(diag(K) > 0) && isequal(am(:).', 1:6) && all(od < 0);
    catch
        ell_um = NaN;   K = nan(6);   kiok = false;   gI = NaN;
    end
end

% ============================================================================
function v = frob_rel(K, Kref)
% ||K - Kref||_F / ||Kref||_F * 100  [%]
    if any(~isfinite(K(:))) || any(~isfinite(Kref(:))), v = NaN; return; end
    v = 100 * norm(K - Kref, 'fro') / norm(Kref, 'fro');
end
