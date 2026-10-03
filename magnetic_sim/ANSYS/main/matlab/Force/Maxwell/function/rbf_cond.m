function out = rbf_cond(model, geom, variant, Rum, opt)
%RBF_COND  Condition number of the Gaussian RBF Gram matrix against rho.
%
%   out = RBF_COND(model, geom, variant, Rum)
%   out = RBF_COND(model, geom, variant, Rum, opt)
%
%   For each rho it forms Phi_ij = exp(-|p_i - p_j|^2 / rho^2) on the nodes with
%   |p| <= Rum and takes ONE symmetric eigendecomposition. Every ridge then
%   follows from the same spectrum, because
%
%       cond(Phi + lam I) = (lam_max + lam) / (lam_min + lam),
%
%   so a sweep over several lam costs nothing extra. Doing it the obvious way --
%   cond() per (rho, lam) -- would be one SVD per pair.
%
%   opt fields
%     .RHO   rho values [um]   (default 10 .. 600, denser where the wall is)
%     .LAM   ridges            (default [0 1e-8 1e-4])
%     .save  write the .mat    (default true)
%
%   The lam = 0 column is reported as lam_max / max(lam_min, eps*lam_max): once
%   Phi is numerically singular its smallest eigenvalue is round-off (it can even
%   go negative), so anything the column shows above ~1e16 means "singular in
%   double precision", not a measured value.
%
%   See also CONV_SERIES, RBF_FIELD.

    if nargin < 5, opt = struct(); end
    g = @(f,d) getfield_(opt, f, d);
    RHO = g('RHO', [10 12 15 18 20 25 30 40 50 60 80 100 150 200 300 400 490 600]);
    LAM = g('LAM', [0 1e-8 1e-4]);
    dosave = g('save', true);

    here = fileparts(mfilename('fullpath'));
    FMX  = fileparts(here);
    ROOT = fileparts(fileparts(FMX));
    CALR = fullfile(ROOT, 'Flux', 'Maxwell');
    addpath(fullfile(CALR,'function'), fullfile(CALR,'common_path'));

    cfg = model_config(model, geom);
    raw = extract_maxwell_data(cfg, 'all', variant);
    ad  = build_actuator_data(raw, cfg);
    inb = find(ad.r2 < (Rum*1e-6)^2);
    P   = ad.Pa(inb,:) * 1e6;                       % um, actuator frame
    Np  = size(P,1);
    s2  = sum(P.^2, 2);
    D2  = max(s2 + s2.' - 2*(P*P.'), 0);   D2(1:Np+1:end) = 0;
    dmin = sqrt(min(D2(D2 > 1e-9)));
    fprintf('[rbf_cond] %s / %s : %d nodes, nearest spacing %.2f um\n', ...
            model, variant, Np, dmin);

    C = nan(numel(RHO), numel(LAM));   EV = nan(numel(RHO), 2);
    fprintf('%7s %12s %12s   %s\n','rho','lam_min','lam_max','cond per lam');
    for k = 1:numel(RHO)
        ev = eig((exp(-D2 / RHO(k)^2) + exp(-D2 / RHO(k)^2).')/2);   % symmetrise
        lmx = max(ev);   lmn = min(ev);
        EV(k,:) = [lmn lmx];
        for j = 1:numel(LAM)
            if LAM(j) == 0
                C(k,j) = lmx / max(lmn, eps*lmx);
            else
                C(k,j) = (lmx + LAM(j)) / (max(lmn,0) + LAM(j));
            end
        end
        fprintf('%7g %12.3e %12.3e   %s\n', RHO(k), lmn, lmx, sprintf('%10.3e',C(k,:)));
    end

    out = struct('model',model, 'geom',geom, 'variant',variant, 'R',Rum, ...
                 'Np',Np, 'dmin',dmin, 'RHO',RHO(:), 'LAM',LAM(:), ...
                 'cond',C, 'ev',EV, 'made',datestr(now,'yyyy-mm-dd HH:MM'));  %#ok<TNOW1,DATST>
    if dosave
        d = fullfile(FMX,'utils','data');
        if ~exist(d,'dir'), mkdir(d); end
        f = fullfile(d, sprintf('rbf_cond_R%d.mat', round(Rum)));
        save(f, '-struct', 'out');
        fprintf('saved %s\n', f);
    end
end

function v = getfield_(s, f, d)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
