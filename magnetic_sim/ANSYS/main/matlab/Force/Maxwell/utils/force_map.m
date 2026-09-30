function out = force_map(recfile, pole, opt)
%FORCE_MAP  Force-error magnitude of a calibrated model on the three planes
%           through the origin, for one excited pole.
%
%   out = FORCE_MAP(recfile, pole)
%   out = FORCE_MAP(recfile, pole, opt)
%
%   recfile : a force calibration record written by main/main.m (FSRC='sph')
%   pole    : which pole is driven at 1 A (1..6), or a VECTOR of poles. With a
%             vector every pole is driven on its own, the error magnitude is
%             taken per pole, and the maps are AVERAGED point by point:
%             mean_k |f_k - f^_k|. The magnitude is taken BEFORE the average --
%             averaging the error vectors first would let opposite errors
%             cancel. This is the same quantity as the histogram's mean.
%   opt.h    grid spacing [um]                      (default 1)
%   opt.save write data/<model>/.mat/emap_*.mat     (default true)
%
%   Both forces are continuous functions of position, so the error is
%   evaluated on a regular grid rather than on the FEM nodes:
%
%       measurement   f  = 0.5 * UF * mgB * grad(b.b)     harmonic model of the
%                                                         FEM field, degree and
%                                                         range from the record
%       model         f^ = (H I)' L(p; l, e) (H I)        calibrated parameters
%       plotted       |f - f^|                            [pN], the norm of the
%                                                         three components
%
%   The three planes are x_a-y_a (z_a = 0), x_a-z_a (y_a = 0) and y_a-z_a
%   (x_a = 0), each restricted to the disc r <= R of the calibration; points
%   outside the disc are NaN. Coordinates are actuator frame, um.
%
%   Nothing is refitted. The harmonic model is rebuilt from the same nodes and
%   the same degree the calibration used, which is why the FEM field is read.
%
%   See also SPH_FIELD, BUILD_L, BASE_MODEL.

    if nargin < 3, opt = struct(); end
    h      = getdef_(opt, 'h', 1);
    dosave = getdef_(opt, 'save', true);

    here = fileparts(mfilename('fullpath'));          % .../Force/Maxwell/utils
    FMX  = fileparts(here);
    ROOT = fileparts(fileparts(FMX));                 % .../main/matlab
    CALR = fullfile(ROOT, 'Flux', 'Maxwell');
    addpath(fullfile(CALR,'function'), fullfile(CALR,'utils'), fullfile(CALR,'common_path'), ...
            fullfile(FMX,'function'), fullfile(FMX,'utils'));

    r = load(recfile);
    assert(isfield(r,'FSRC') && strcmpi(r.FSRC,'sph'), 'force_map:record', ...
           'the record is not an FSRC=''sph'' calibration');
    assert(strcmp(r.base,'current'), 'force_map:base', ...
           'current-base records only: the excitation is a coil current');
    validateattributes(pole, {'numeric'}, {'vector','integer','>=',1,'<=',6}, mfilename, 'pole');
    pole = unique(pole(:).');
    R = r.R_select;

    % ---- the measurement model, rebuilt exactly as the calibration built it ---
    cfg = model_config(r.model, r.GEOM);
    raw = extract_maxwell_data(cfg, 'all', r.VARIANT);
    ad  = build_actuator_data(raw, cfg);
    inb = find(ad.r2 < (R*1e-6)^2);
    [~, gq, si] = sph_field(R, struct('P',ad.Pa(inb,:)*1e6,'B',ad.Ba(inb,:,:)), ...
                            struct('L',r.LDEG,'quiet',true,'condmax',1500));
    assert(abs(si.NMAE_all - r.NMAE_sph) < 1e-9, 'force_map:sph', ...
           'harmonic model differs from the one the record was calibrated on (NMAE %.6f vs %.6f)', ...
           si.NMAE_all, r.NMAE_sph);

    g = (-R:h:R).';   [A, B] = ndgrid(g, g);   msk = hypot(A, B) <= R + 1e-9;
    a = A(msk);   b = B(msk);   z = zeros(size(a));
    PL  = {[a b z], [a z b], [z a b]};                % x_a-y_a | x_a-z_a | y_a-z_a
    NAM = {'xy','xz','yz'};   AXN = {'x_a','y_a'; 'x_a','z_a'; 'y_a','z_a'};

    out = struct('recfile',recfile, 'model',r.model, 'R',R, 'LDEG',r.LDEG, 'pole',pole, ...
                 'h',h, 'g',g, 'npts_cal',r.npts, 'l_hat',r.l_hat, 'USE_BIAS',r.USE_BIAS, ...
                 'unit','pN', 'made',datestr(now,'yyyy-mm-dd HH:MM'));   %#ok<TNOW1,DATST>
    emax = 0;   emin = inf;
    for k = 1:3
        P  = PL{k};
        Lk = build_L(P, r.l_hat, r.e_hat, r.Pc_base);        % shared by every pole
        ek = zeros(size(P,1), numel(pole));   fk = ek;
        for q = 1:numel(pole)
            I  = zeros(6,1);   I(pole(q)) = 1;
            f  = (0.5 * r.mgB * r.UF * gq(P, I, 'du')).';                  % 3 x N
            fh = base_model(Lk, r.H_hat, I, 1);
            ek(:,q) = vecnorm(fh(:,:,1) - f, 2, 1).';
            fk(:,q) = vecnorm(f, 2, 1).';
        end
        e  = mean(ek, 2);                                    % magnitude first, then mean
        E  = nan(size(A));   E(msk) = e;
        F  = nan(size(A));   F(msk) = mean(fk, 2);
        out.(NAM{k}) = struct('E',E, 'F',F, 'ax',{AXN(k,:)}, 'mean',mean(e), 'max',max(e), ...
                              'min',min(e), 'mean_pole',mean(ek,1), 'max_pole',max(ek,[],1));
        emax = max(emax, max(e));   emin = min(emin, min(e));
        fprintf('[force_map] P%s, %s-%s plane: %d points, |error| mean %.4f pN, min %.4f, max %.4f pN%s', ...
                sprintf('%d',pole), AXN{k,1}, AXN{k,2}, numel(e), mean(e), min(e), max(e), newline);
    end
    out.emax = emax;   out.emin = emin;
    if dosave
        tg = 'single';  if r.USE_BIAS, tg = 'eighteen'; end
        fo = fullfile(FMX, 'data', r.model, '.mat', ...
                      sprintf('emap_R%d_P%s_%s.mat', round(R), ptag_(pole), tg));
        save(fo, '-struct', 'out');
        fprintf('saved %s%s', fo, newline);
    end
end

function t = ptag_(pole)
%PTAG_  file tag: the pole number, or 'all' when all six are averaged
    if isequal(pole, 1:6), t = 'all'; else, t = sprintf('%d', pole); end
end

function v = getdef_(s, f, d)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
