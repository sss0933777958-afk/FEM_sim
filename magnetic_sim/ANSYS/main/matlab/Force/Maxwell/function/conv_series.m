function out = conv_series(model, geom, variant, Rum, opt)
%CONV_SERIES  Walk the axes-shell ladder and record BOTH calibrations per level.
%
%   out = CONV_SERIES(model, geom, variant, Rum)
%   out = CONV_SERIES(model, geom, variant, Rum, opt)
%
%   For every ladder level Nr = 1..Nrmax (N = 6*Nr+1 sample points) this runs
%     flux   fitting -> solve_current      on the trilinearly sampled field
%     force  fitting_force -> solve_force  on f = 0.5*mgB*grad(b.b) from a
%                                          harmonic model of the same export
%   and stores the four trends against the point count, so a convergence figure
%   can be drawn later without re-running anything (plot-scripts-pure).
%
%   The two sides are given the SAME points: the flux call returns them in
%   metres, the force call in um, and the two are cross-checked before use.
%
%   opt fields (all optional)
%     .LDEG     harmonic degree for the measurement model (default 9)
%     .Nrmax    highest ladder level                      (default 40)
%     .USE_BIAS eighteen-parameter fit                    (default true)
%     .mgB      particle magnetisation constant [A.um^2/mT] (default 0.0451)
%     .UF       force unit factor                         (default 1e3)
%     .save     write the .mat                            (default true)
%
%   The record also carries the FULL-NODE calibration of both sides, which is
%   what a convergence figure draws as its reference line.
%
%   Model-agnostic on purpose: it takes the model/geom/variant, so it belongs in
%   utils/ root rather than under a per-model folder.
%
%   See also CONV_DESIGN_WS, FITTING, SOLVE_CURRENT, FITTING_FORCE, SOLVE_FORCE.

    if nargin < 5, opt = struct(); end
    g = @(f,d) getfield_(opt, f, d);
    LDEG = g('LDEG', 9);   Nrmax = g('Nrmax', 40);   USE_BIAS = g('USE_BIAS', true);
    mgB  = g('mgB', 0.0451);   UF = g('UF', 1e3);    dosave = g('save', true);

    here = fileparts(mfilename('fullpath'));          % .../Force/Maxwell/utils
    FMX  = fileparts(here);
    ROOT = fileparts(fileparts(FMX));                 % .../main/matlab
    CALR = fullfile(ROOT, 'Flux', 'Maxwell');
    addpath(fullfile(CALR,'function'), fullfile(CALR,'utils'), fullfile(CALR,'common_path'), ...
            fullfile(FMX,'function'), fullfile(FMX,'function'));

    % ---- data, one harmonic fit, one excitation set -------------------------
    cfg = model_config(model, geom);
    assert(all(cfg.s_source == 1), 'conv_series:sSource', ...
           'raw is not all-source; b = sum_k I_k b_k would superpose flipped fields');
    raw = extract_maxwell_data(cfg, 'all', variant);
    ad  = build_actuator_data(raw, cfg);
    assert(isequal(ad.F, eye(6)), 'conv_series:Fmap', 'coil->pole map is not identity');
    inb = find(ad.r2 < (Rum*1e-6)^2);
    Pn  = ad.Pa(inb,:)*1e6;   Bn = ad.Ba(inb,:,:);
    [~, gq, si] = sph_field(Rum, struct('P',Pn,'B',Bn), ...
                            struct('L',LDEG,'quiet',true,'condmax',1500));
    pr = nchoosek(1:6,2);   u = [eye(6), zeros(6,size(pr,1))];
    for q = 1:size(pr,1), u(pr(q,:), 6+q) = 1; end
    Pc = ad.Pc_base;   F = ad.F;
    fprintf('[conv_series] %s / %s : %d nodes, sph L=%d NMAE %.4f %%%s', ...
            model, variant, numel(inb), si.L, si.NMAE_all, newline);

    % ---- full-node reference (what the figure draws as the dashed line) -----
    [Pm, Bst, nfull] = cfg.select_ball(ad, Rum*1e-6);
    [eR, lR] = fitting(Pm, Bst, Pc, 0.5e-3, USE_BIAS);
    [~, gR]  = solve_current(lR, eR, Pc, Pm, Bst, F, [], Pm, Bst);
    fR = fmeas_(gq, u, Pm*1e6, mgB, UF);
    [lFR, eFR, HFR] = fitting_force(fR, Pm*1e6, Pc, 500, eye(6), u, USE_BIAS);
    [~, gFR] = solve_force(HFR, mgB, lFR, UF, 'current');
    [~, ~, gBR] = solve_force(HFR, mgB, lFR, UF, 'current');
    ref = struct('npts',nfull, 'l_flux',lR*1e6, 'g_flux',gR, ...
                 'fg_flux',UF*(mgB/(2*lR*1e6))*gR^2, ...
                 'l_force',lFR, 'fg_force',gFR, 'g_force',gBR);
    fprintf('[conv_series] full node (N=%d): flux l=%.3f fg=%.5f | force l=%.3f fg=%.5f%s', ...
            nfull, ref.l_flux, ref.fg_flux, ref.l_force, ref.fg_force, newline);

    % ---- the ladder ---------------------------------------------------------
    wsopt = struct('model',model, 'geom',geom, 'variant',variant, 'frame','actuator');
    N = nan(Nrmax,1);  LF = N;  GF = N;  LC = N;  GC = N;  BF = N;  BC = N;
    for Nr = 1:Nrmax
        [P, Bs] = conv_design_ws(Nr, Rum*1e-6, wsopt);          % metres, with field
        [~,~,sp] = conv_design_ws(Nr, Rum, ...
                        struct('points_only',true,'R_act',ad.R_act,'quiet',true));
        % same design, two unit systems -- verify rather than assume. (The
        % field-engine call returns the points as its first output; only the
        % points_only path fills info.P_act.)
        assert(max(abs(P*1e6 - sp.P_act), [], 'all') < 1e-6, ...
               'conv_series:points', 'the flux and force point sets differ at Nr=%d', Nr);
        N(Nr) = size(P,1);

        [ef, lf]     = fitting(P, Bs, Pc, 0.5e-3, USE_BIAS);
        [~, gf]      = solve_current(lf, ef, Pc, P, Bs, F, [], P, Bs);
        LF(Nr) = lf*1e6;   BF(Nr) = gf;                        % B g_I [mT/A]
        GF(Nr) = UF*(mgB/(2*LF(Nr)))*gf^2;                     % flux's implied F g

        fm           = fmeas_(gq, u, sp.P_act, mgB, UF);
        [lc, ec, Hc] = fitting_force(fm, sp.P_act, Pc, 500, eye(6), u, USE_BIAS);
        [~, gc, bc]  = solve_force(Hc, mgB, lc, UF, 'current');
        LC(Nr) = lc;   GC(Nr) = gc;   BC(Nr) = bc;
        fprintf('  Nr=%2d N=%4d | flux l=%8.3f g=%8.4f | force l=%8.3f g=%8.4f%s', ...
                Nr, N(Nr), LF(Nr), BF(Nr), LC(Nr), BC(Nr), newline);
    end

    out = struct('model',model, 'geom',geom, 'variant',variant, 'R',Rum, ...
                 'LDEG',LDEG, 'K_sph',si.K, 'NMAE_sph',si.NMAE_all, ...
                 'USE_BIAS',USE_BIAS, 'mgB',mgB, 'UF',UF, 'nexc',size(u,2), ...
                 'npts',N, 'l_flux',LF, 'fg_flux',GF, 'g_flux',BF, ...
                 'l_force',LC, 'fg_force',GC, 'g_force',BC, ...
                 'ref',ref, 'made',datestr(now,'yyyy-mm-dd HH:MM'));   %#ok<TNOW1,DATST>
    if dosave
        tag = 'single';  if USE_BIAS, tag = 'eighteen'; end
        d = fullfile(FMX,'utils','data');
        if ~exist(d,'dir'), mkdir(d); end
        f = fullfile(d, sprintf('conv_R%d_%s.mat', round(Rum), tag));
        save(f, '-struct', 'out');
        fprintf('saved %s%s', f, newline);
    end
end

function f = fmeas_(gq, u, P, mgB, UF)
    Np = size(P,1);  NI = size(u,2);  f = zeros(3,Np,NI);
    for j = 1:NI, f(:,:,j) = (0.5*mgB*UF*gq(P,u(:,j),'du')).'; end
end

function v = getfield_(s, f, d)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
