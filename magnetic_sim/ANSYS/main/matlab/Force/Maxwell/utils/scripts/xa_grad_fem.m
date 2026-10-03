%% xa_grad_fem.m -- b.b and d(b.b)/dx_a on the x_a axis, straight from the Maxwell
%%                   export for P1, by our own trilinear interpolation.
%
%  MESH picks which solve to read.  Both were exported on the SAME Cartesian grid
%  (x,y in +-0.6 mm, z in [-13.31,-12.11] mm, pitch 0.02 mm, 61^3 = 226,981 nodes,
%  P1 excited), so their curves are directly comparable; only the solver's
%  tetrahedral mesh differs.
%
%      MESH = '0.1'    grad(0.1)\B_p1_grad.fld   + B_p1_test.fld     tag h0p1
%      MESH = '0.02'   B_p1_grad_0.02.fld        (no matching B)     tag h0p02
%
%  OUTPUTS.  Drawn by plot/long2016_hexapole_halfcut/xa_grad_fem.m (its QTY switch):
%
%      g_xa  = d(b.b)/dx_a    [mT^2/um]   <- interpolate AEDT's Grad(Dot(B,B)),
%                                            THEN dot with the x_a unit vector
%      bb    = b.b            [mT^2]      <- interpolate B, THEN square and sum
%      res   = p - bb         [mT^2]      <- degree-DEG least-squares polynomial in
%                                            xt = s/SPAN, minus the curve it was
%                                            fitted to (i.e. the LS residual)
%
%  bb / res / the fit / the two cross-checks need the B export of the same solve.
%  Only the 0.1 mm run has one, so with MESH = '0.02' the script produces g_xa
%  alone and says so.
%
%  ORDER.  Interpolate the raw exported field first, form the derived quantity
%  second.  Trilinear interpolation does not commute with squaring, so
%  interpolating b and then squaring is NOT the same as interpolating a nodal b.b.
%
%  WHICH SOLVE.  'maxwell_test' (the _test exports of 2026-09-04) for MESH='0.1',
%  NOT the 'maxwell' baseline of 2026-07-31.  Two independent signs, established
%  2026-09-20:
%    - every B_p<k>_grad.fld was written about a minute after its B_p<k>_test.fld,
%      in one AEDT session on 2026-09-04;
%    - central differences of b.b taken from B_p1_test.fld reproduce the exported
%      gradient about twice as well as the same differences taken from
%      baseline\B_p1.fld (median 5.41 % vs 12.46 % over the interior; 2.65 % vs
%      4.49 % inside r <= 150 um).
%  The two solves are not the same field: their b.b differs by a median of -0.97 %
%  (p95 8.5 %), |B|max 2.6464 vs 2.7223 T.  So these curves are NOT strictly
%  comparable with models fitted to the baseline export.
%  MESH='0.02' is the run solved 2026-09-20 22:15 and exported 22:31.
%
%  UNITS.  Grad(Dot(B,B)) comes out in T^2/m, which is numerically identical to
%  mT^2/um:  1 T^2/m = (1e3 mT)^2 / (1e6 um) = 1 mT^2/um.  The gradient is
%  therefore never rescaled; only B is (T -> mT).
%
%  THE POINTS.  NQ equidistant points on the x_a axis of the ACTUATOR frame,
%  P_act = (s,0,0), |s| <= SPAN, mapped back to the measure frame only because that
%  is what indexes the lattice:
%      P_meas = P_act * R_act + [0 0 SPH_OFST]        (row-vector form)
%  b.b is a SCALAR and so frame-independent -- only the points carry the actuator
%  frame.  The gradient is a vector, projected at the very end,
%  d(b.b)/dx_a = G_meas . xhat_a with xhat_a = R_act(1,:) = the P1 pole axis.  The
%  field itself is never rotated.
%
%  INTERPOLATION.  The project's own regular-grid trilinear.  trilerp_ and the
%  lattice builder are reproduced here line for line from
%  Flux/Maxwell/function/conv_design_ws.m, which is the standing rule for this
%  engine (no shared loader; each caller carries its own copy).
%
%  Output: data/long2016_hexapole_halfcut/.mat/xa_fem_<tag>.mat
%
%  [ADDED 2026-09-20]

clearvars;  clc;

MODEL = 'long2016_hexapole_halfcut';   GEOM = 'tip40um';
MESH  = '0.02';                        % '0.1' | '0.02' -- solver mesh of the run to read
POLE  = 1;                             % x_a <- P1
SPAN  = 150;                           % |s| <= SPAN [um]
NQ    = 1001;                          % equidistant points on the axis
DEG   = 7;                             % polynomial order fitted to b.b (needs the B export)
XCHK  = true;                          % cross-check the gradient against central differences
WRITEPTS = true;                       % also write the sample points as an AEDT point list
NQ_EXTRA = 15001;                      % [] = off.  A SECOND, denser point list on the same
                                       %   axis, for the element-order test: 15001 points over
                                       %   +-150 um = 0.02 um spacing, about 128 samples per
                                       %   2.55 um slope-break interval (the 1001-point list
                                       %   gives only ~8, too few to tell a line from a
                                       %   parabola).  Written only -- nothing reads it here.

switch MESH
    case '0.1',  GRADFILE = fullfile('grad(0.1)','B_p1_grad.fld');  BFILE = 'B_p1_test.fld';  TAG = 'h0p1';
                 PTSOUT  = '';   DOTOUT = '';
    case '0.02', GRADFILE = 'B_p1_grad_0.02.fld';                   BFILE = '';               TAG = 'h0p02';
                 PTSOUT  = 'B_p1_grad_0.02_pts.fld';   % AEDT Grad(Dot(B,B)) AT the points
                 DOTOUT  = 'B_p1_dot_0.02_pts.fld';    % AEDT Dot(B,B)       AT the points
                                                       %   (overwritten 2026-09-21 with
                                                       %    Smooth(Grad(Smooth(Dot))) -- the
                                                       %    header decides how it is read)
                 BPTSOUT = 'B_p1_0.02_pts.fld';        % AEDT B_Vector       AT the points
    otherwise,   error('xa_grad_fem:mesh','MESH must be ''0.1'' or ''0.02''');
end

HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');

cfg  = model_config(MODEL, GEOM);
RACT = cfg.R_act;
xhat = RACT(POLE,:).';                 % the x_a unit vector, in measure components
fprintf(['mesh %s mm  ->  tag %s' newline], MESH, TAG);
fprintf(['xhat_a (measure) = [%.6f %.6f %.6f],  |xhat_a| = %.12f' newline], xhat, norm(xhat));

%% ---- the points, defined in the ACTUATOR frame --------------------------------
s  = linspace(-SPAN, SPAN, NQ).';
Pa = [s, zeros(NQ,2)];                             % [um], actuator frame
Pg = (Pa*1e-6) * RACT;                             % -> measure [m]
Pg(:,3) = Pg(:,3) + cfg.SPH_OFST;                  % -> global z

% The same points as an AEDT point list, so the Fields Calculator can evaluate
% Grad(Dot(B,B)) AT them instead of on a Cartesian grid we then interpolate.  It
% goes next to the .fld exports because that is where AEDT's file dialogs open and
% where its output will land.  Format: one "Unit=" header, then x y z per line, in
% the measure/global frame.  Independent of MESH, so either run rewrites the same
% file.
if WRITEPTS
    pts = fullfile(cfg.fld_dir, sprintf('xa_%d.pts', NQ));
    fid = fopen(pts, 'w');
    assert(fid > 0, 'xa_grad_fem:pts', 'cannot write %s', pts);
    fprintf(fid, 'Unit=mm\n');
    fprintf(fid, '%.10f %.10f %.10f\n', (Pg*1e3).');
    fclose(fid);
    nl = numel(strsplit(strtrim(fileread(pts)), newline));
    assert(nl == NQ+1, 'xa_grad_fem:ptsN', '%s has %d lines, expected %d', pts, nl, NQ+1);
    fprintf(['  [pts] %s  (%d points + header)' newline], pts, NQ);
end

% [ADDED 2026-09-21] the denser list.  Identical construction; only the count differs.
if ~isempty(NQ_EXTRA)
    se = linspace(-SPAN, SPAN, NQ_EXTRA).';
    Pe = ([se, zeros(NQ_EXTRA,2)]*1e-6) * RACT;
    Pe(:,3) = Pe(:,3) + cfg.SPH_OFST;
    pts2 = fullfile(cfg.fld_dir, sprintf('xa_%d.pts', NQ_EXTRA));
    fid = fopen(pts2, 'w');
    assert(fid > 0, 'xa_grad_fem:pts2', 'cannot write %s', pts2);
    fprintf(fid, 'Unit=mm\n');
    fprintf(fid, '%.10f %.10f %.10f\n', (Pe*1e3).');
    fclose(fid);
    % read back: the file itself must reproduce the axis to well under the spacing
    Mb = readmatrix(pts2, 'FileType','text', 'NumHeaderLines',1);
    assert(size(Mb,1) == NQ_EXTRA, 'xa_grad_fem:pts2N', ...
           '%s has %d points, expected %d', pts2, size(Mb,1), NQ_EXTRA);
    sb  = (Mb*1e-3 - [0 0 cfg.SPH_OFST]) * RACT(POLE,:).' * 1e6;    % back to s [um]
    off = max(abs(Mb*1e-3 - (Pe)), [], 'all');
    fprintf(['  [pts] %s  (%d points + header)' newline], pts2, NQ_EXTRA);
    fprintf(['        spacing %.4f um, round-trip s error %.2e um, write error %.2e m' newline], ...
            median(diff(sb)), max(abs(sb - se)), off);
    fprintf(['        %.1f samples per %.2f um break interval (the 1001-point list gives %.1f)' ...
             newline], 2.55/median(diff(sb)), 2.55, 2.55/(2*SPAN/(NQ-1)));
end

%% ---- gradient: interpolate the VECTOR, then project on the axis ----------------
gpath = fullfile(cfg.fld_dir, GRADFILE);
assert(isfile(gpath), 'xa_grad_fem:nofileG', 'gradient export not found: %s', gpath);
rg  = import_maxwell_fld(gpath);                   % [T^2/m] = [mT^2/um], measure frame
lat = build_lattice_(rg.x, rg.y, rg.z, [rg.bx, rg.by, rg.bz]);
fprintf(['  [lattice] %d x %d x %d, pitch %.1f/%.1f/%.1f um' newline], lat.n, lat.h*1e6);
[Gq, inG] = trilerp_(lat, Pg(:,1), Pg(:,2), Pg(:,3));
assert(all(inG), 'xa_grad_fem:outboxG', '%d of %d points fall outside the gradient box', ...
       nnz(~inG), NQ);
G_meas = Gq;                                       % NQ x 3 [mT^2/um], measure frame
g_xa   = G_meas * xhat;                            % d(b.b)/dx_a  [mT^2/um]
G_act  = G_meas * RACT.';                          % all three actuator components

fprintf([newline 'd(b.b)/dx_a on |s| <= %g um, %d points, P%d excited' newline], SPAN, NQ, POLE);
fprintf(['  range   %+.4f .. %+.4f mT^2/um    (at s=0  %+.4f)' newline], ...
        min(g_xa), max(g_xa), interp1(s,g_xa,0));
fprintf(['  the other two actuator components, for scale:' newline]);
fprintf(['    d/dy_a  %+.4f .. %+.4f     d/dz_a  %+.4f .. %+.4f mT^2/um' newline], ...
        min(G_act(:,2)), max(G_act(:,2)), min(G_act(:,3)), max(G_act(:,3)));

%% ---- the same gradient, evaluated BY AEDT at the very same points ---------------
%  The Fields Calculator was given xa_<NQ>.pts and asked for the same
%  Grad(Dot(B_Vector,B_Vector)), so this route skips BOTH the Cartesian export and
%  our trilinear interpolation.  It is not "uninterpolated" -- AEDT interpolates in
%  the tetrahedral mesh instead -- so the comparison measures one interpolation
%  against the other, on one identical set of points.
g_pts = [];   PCHK = struct('done',false);
if ~isempty(PTSOUT) && isfile(fullfile(cfg.fld_dir, PTSOUT))
    ppath = fullfile(cfg.fld_dir, PTSOUT);
    rp = import_maxwell_fld(ppath);
    assert(numel(rp.x) == NQ, 'xa_grad_fem:ptsN', ...
           '%s has %d points, expected %d', ppath, numel(rp.x), NQ);
    dev = max(max(abs([rp.x rp.y rp.z] - Pg)));
    assert(dev < 1e-9, 'xa_grad_fem:ptsXYZ', ...
           'the exported sample points differ from ours by up to %.3e m (order changed?)', dev);
    g_pts = [rp.bx rp.by rp.bz] * xhat;            % [mT^2/um]
    rel   = (g_xa - g_pts) ./ g_pts;
    PCHK  = struct('done',true, 'src',ppath, 'max_dev_m',dev, ...
                   'median_pct',100*median(abs(rel)), 'p90_pct',100*prctile(abs(rel),90), ...
                   'max_pct',100*max(abs(rel)), ...
                   'I_pts',trapz(s,g_pts), 'I_grid',trapz(s,g_xa));
    PCHK.I_rel_pct = 100*(PCHK.I_grid - PCHK.I_pts)/PCHK.I_pts;
    fprintf([newline 'AEDT evaluated at the points (no Cartesian grid, no trilinear)' newline]);
    fprintf(['  sample points match ours to %.1e m' newline], dev);
    fprintf(['  range   %+.4f .. %+.4f mT^2/um    (at s=0  %+.4f)' newline], ...
            min(g_pts), max(g_pts), interp1(s,g_pts,0));
    fprintf(['  grid-vs-points |rel|: median %.2f %%, p90 %.2f %%, max %.2f %%' newline], ...
            PCHK.median_pct, PCHK.p90_pct, PCHK.max_pct);
    fprintf(['  but the integrals agree: %.4f vs %.4f mT^2  (%+.3f %%)' newline], ...
            PCHK.I_grid, PCHK.I_pts, PCHK.I_rel_pct);
end

%% ---- b.b, also evaluated BY AEDT at the very same points ------------------------
%  A SCALAR export, so 4 columns (x y z value) and value in T^2 -- import_maxwell_fld
%  only reads 6-column vector files, hence the local reader.  Fitted with the same
%  degree-DEG polynomial as bb so the two smoothness numbers are on one scale.
bb_pts = [];  acof_pts = [];  p_pts = [];  res_pts = [];  DCHK = struct('done',false);
g_pts_sm = [];  SCHK = struct('done',false);
% [ADDED 2026-09-21] The role of this file is decided by its HEADER, not by its name.
%   'Dot(B_Vector, B_Vector)'                      -> scalar b.b at the points   (bb_pts)
%   '...Grad(...)' (e.g. Smooth(Grad(Smooth(Dot))))-> gradient vector at the points (g_pts_sm)
% B_p1_dot_0.02_pts.fld held the first on 2026-09-20 and was overwritten with the second on
% 2026-09-21; reading it by name would have put a gradient into the b.b slot.
DOTROLE = '';
if ~isempty(DOTOUT) && isfile(fullfile(cfg.fld_dir, DOTOUT))
    dq = fld_quantity_(fullfile(cfg.fld_dir, DOTOUT));
    if ~isempty(strfind(dq, 'Grad'))                                        %#ok<STREMP>
        DOTROLE = 'grad';
    elseif ~isempty(strfind(dq, 'Dot'))                                     %#ok<STREMP>
        DOTROLE = 'dot';
    end
    fprintf([newline '%s' newline '  header quantity: %s  -> role ''%s''' newline], ...
            DOTOUT, dq, DOTROLE);
end
if strcmp(DOTROLE, 'grad')
    spath = fullfile(cfg.fld_dir, DOTOUT);
    rs = import_maxwell_fld(spath);
    assert(numel(rs.x) == NQ, 'xa_grad_fem:smN', ...
           '%s has %d points, expected %d', spath, numel(rs.x), NQ);
    dev = max(max(abs([rs.x rs.y rs.z] - Pg)));
    assert(dev < 1e-9, 'xa_grad_fem:smXYZ', ...
           'the exported sample points differ from ours by up to %.3e m', dev);
    g_pts_sm = [rs.bx rs.by rs.bz] * xhat;         % [mT^2/um] = [T^2/m]
    SCHK = struct('done',true, 'src',spath, 'quantity',dq, 'max_dev_m',dev);
    fprintf(['d(b.b)/dx_a, smoothed, evaluated by AEDT at the points' newline]);
    fprintf(['  range   %+.4f .. %+.4f mT^2/um   (at s=0  %+.4f)' newline], ...
            min(g_pts_sm), max(g_pts_sm), interp1(s, g_pts_sm, 0));
    fprintf(['  monotonic: %s ; largest step between neighbours %.4f (typical %.4f)' newline], ...
            string(all(diff(g_pts_sm)>0)), max(abs(diff(g_pts_sm))), median(abs(diff(g_pts_sm))));
    if ~isempty(g_pts)
        rr = (g_pts_sm - g_pts) ./ g_pts;
        SCHK.rel_med_pct = 100*median(abs(rr));   SCHK.rel_max_pct = 100*max(abs(rr));
        SCHK.step_ratio = median(abs(diff(g_pts_sm))) / median(abs(diff(g_pts)));
        fprintf(['  vs the unsmoothed point export: median |rel| %.3f %%, max %.3f %%, ' ...
                 'neighbour step x%.3f' newline], ...
                SCHK.rel_med_pct, SCHK.rel_max_pct, SCHK.step_ratio);
    end
end
if strcmp(DOTROLE, 'dot')
    dpath = fullfile(cfg.fld_dir, DOTOUT);
    [Pd, vd] = read_scalar_fld_(dpath);
    assert(numel(vd) == NQ, 'xa_grad_fem:dotN', ...
           '%s has %d points, expected %d', dpath, numel(vd), NQ);
    dev = max(max(abs(Pd - Pg)));
    assert(dev < 1e-9, 'xa_grad_fem:dotXYZ', ...
           'the exported sample points differ from ours by up to %.3e m', dev);
    bb_pts   = vd * 1e6;                           % T^2 -> mT^2
    xt       = s / SPAN;
    Vd       = xt .^ (0:DEG);
    acof_pts = Vd \ bb_pts;
    p_pts    = Vd * acof_pts;
    res_pts  = p_pts - bb_pts;
    DCHK = struct('done',true, 'src',dpath, 'max_dev_m',dev, ...
                  'max_rel_pct',100*max(abs(res_pts./bb_pts)), ...
                  'rms_rel_pct',100*rms(res_pts./bb_pts));
    fprintf([newline 'b.b evaluated by AEDT at the points' newline]);
    fprintf(['  range   %.4f .. %.4f mT^2      (at s=0  %.4f)' newline], ...
            min(bb_pts), max(bb_pts), interp1(s,bb_pts,0));
    fprintf(['  monotonic: %s ; largest step between neighbours %.4f mT^2 (typical %.4f)' newline], ...
            string(all(diff(bb_pts)>0)), max(abs(diff(bb_pts))), median(abs(diff(bb_pts))));
    fprintf(['  degree-%d residual: %+.4f .. %+.4f mT^2  (max |rel| %.4f %%, rms %.4f %%)' newline], ...
            DEG, min(res_pts), max(res_pts), DCHK.max_rel_pct, DCHK.rms_rel_pct);
    if ~isempty(g_pts)                             % the two point exports must be consistent
        DCHK.D = bb_pts(end) - bb_pts(1);
        DCHK.I = trapz(s, g_pts);
        DCHK.ftc_rel_pct = 100*(DCHK.I - DCHK.D)/DCHK.D;
        fprintf(['  [FTC vs the point-wise gradient] b.b(+%g)-b.b(-%g) = %.4f  vs  ' ...
                 'trapz = %.4f mT^2  (%+.4f %%)' newline], SPAN, SPAN, DCHK.D, DCHK.I, DCHK.ftc_rel_pct);
    end
end

%% ---- B itself, evaluated BY AEDT at the very same points ------------------------
%  [ADDED 2026-09-21] Answers a specific question: inside one mesh element, does AEDT's
%  B vary linearly (a sloped straight line) or not at all (a horizontal line)?  First-order
%  nodal edge/vector potential elements give a B that is CONSTANT per tetrahedron; a nodal-
%  averaged (smoothed) output is linear instead.  The test is mechanical:
%    - piecewise constant -> long runs with d|b|/ds == 0 (flat)
%    - piecewise linear   -> no flats, the second difference is ~0 inside a run and spikes
%                            where the sample crosses a face
B_pts = [];  bmag_pts = [];  s_b = [];  SEGCHK = struct('done',false);
if ~isempty(BPTSOUT) && isfile(fullfile(cfg.fld_dir, BPTSOUT))
    bpath = fullfile(cfg.fld_dir, BPTSOUT);
    bq = fld_quantity_(bpath);
    assert(~isempty(strfind(bq, 'B_Vector')) && isempty(strfind(bq, 'Grad')) ...
           && isempty(strfind(bq, 'Dot')), 'xa_grad_fem:bptsQty', ...
           '%s holds ''%s'', not the plain B_Vector', bpath, bq);                %#ok<STREMP>
    rb = import_maxwell_fld(bpath);
    % [MODIFIED 2026-09-21] This export may use EITHER point list (1001 @ 0.3 um or 15001 @
    %   0.02 um), so its s axis is recovered from the coordinates in the file rather than
    %   compared against Pg.  Two guards replace that comparison: the points must lie on the
    %   x_a axis (y_a, z_a = 0) and be equally spaced.
    Pb = [rb.x rb.y rb.z];   Pb(:,3) = Pb(:,3) - cfg.SPH_OFST;
    Pab = (Pb * RACT.') * 1e6;                     % measure [m] -> actuator [um]
    s_b = Pab(:,1);
    dev = max(max(abs(Pab(:,2:3))));               % off-axis excursion [um]
    assert(dev < 1e-6, 'xa_grad_fem:bptsAxis', ...
           'the exported points leave the x_a axis by up to %.3e um', dev);
    assert(max(abs(diff(diff(s_b)))) < 1e-6, 'xa_grad_fem:bptsStep', ...
           'the exported points are not equally spaced');
    B_pts    = [rb.bx rb.by rb.bz] * 1e3;          % T -> mT
    bmag_pts = vecnorm(B_pts, 2, 2);               % |b| [mT]
    fprintf(['  %d points, spacing %.4f um, off-axis %.2e um' newline], ...
            numel(s_b), median(diff(s_b)), dev);
    s = s_b;                                       % local alias for the diagnostic below
    d1 = diff(bmag_pts);   d2 = diff(d1);
    flat = mean(abs(d1) < 0.02*median(abs(d1)));   % share of (near) zero-slope steps
    jmp  = find(abs(d2) > 8*median(abs(d2)));      % slope breaks = face crossings
    [~, kk] = max(diff(jmp));   i1 = jmp(kk)+2;   i2 = jmp(kk+1);
    lin = bmag_pts(i1) + (s(i1:i2)-s(i1)) * (bmag_pts(i2)-bmag_pts(i1)) / (s(i2)-s(i1));
    SEGCHK = struct('done',true, 'src',bpath, 'quantity',bq, 'max_dev_m',dev, ...
                    'flat_frac',flat, 'n_breaks',numel(jmp), ...
                    'break_gap_um',median(diff(s(jmp))), ...
                    'seg_s',[s(i1) s(i2)], 'seg_n',i2-i1+1, ...
                    'seg_dev_mT',max(abs(bmag_pts(i1:i2)-lin)), ...
                    'seg_span_mT',abs(bmag_pts(i2)-bmag_pts(i1)));
    SEGCHK.seg_dev_pct = 100*SEGCHK.seg_dev_mT/SEGCHK.seg_span_mT;
    fprintf([newline 'B evaluated by AEDT at the points  (%s)' newline], bq);
    fprintf(['  |b| range %.4f .. %.4f mT' newline], min(bmag_pts), max(bmag_pts));
    fprintf(['  zero-slope steps %.2f %% ; %d slope breaks, typical spacing %.2f um' newline], ...
            100*flat, numel(jmp), SEGCHK.break_gap_um);
    fprintf(['  longest run s = %.2f .. %.2f um (%d pts): deviation from a straight line ' ...
             '%.3e mT = %.2f %% of its %.3e mT change' newline], ...
            s(i1), s(i2), SEGCHK.seg_n, SEGCHK.seg_dev_mT, SEGCHK.seg_dev_pct, SEGCHK.seg_span_mT);
    s = linspace(-SPAN, SPAN, NQ).';           % restore the main axis
end

%% ---- where the trilinear interpolant kinks -------------------------------------
%  Trilinear interpolation is C0: its derivative jumps whenever the sample crosses a
%  cell face.  Along this axis the position is linear in s, so each face family gives
%  an exactly periodic set of s values -- solved here rather than detected, because
%  detection misses the weak ones.
%      P_meas(s) = s*1e-6 * xhat + [0 0 SPH_OFST]
%      face of axis d  <=>  (P_d - o_d)/h_d  is an integer
KINK = struct('s',{{}}, 'period',[], 'name',{{'x','y','z'}});
off  = [0 0 cfg.SPH_OFST];
for d = 1:3
    A = xhat(d)*1e-6 / lat.h(d);                   % d(face index)/ds   [1/um]
    B = (off(d) - lat.o(d)) / lat.h(d);
    if abs(A) < eps                                % the axis lies IN this family of faces
        KINK.s{d} = [];   KINK.period(d) = Inf;
        fprintf(['  [kink] %s: the whole axis lies on one face (no crossings)' newline], KINK.name{d});
        continue
    end
    m  = ceil(min(A*[-SPAN SPAN])+B) : floor(max(A*[-SPAN SPAN])+B);
    sk = (m - B)/A;
    KINK.s{d} = sort(sk(abs(sk) <= SPAN + 1e-9)).';
    KINK.period(d) = abs(lat.h(d)*1e6 / xhat(d));
    fprintf(['  [kink] %s faces: %d crossings, every %.4f um' newline], ...
            KINK.name{d}, numel(KINK.s{d}), KINK.period(d));
end

%% ---- B-dependent part: b.b, the polynomial fit, and the two cross-checks --------
bb = [];  B_meas = [];  B_act = [];  xt = [];  acof = [];  p = [];  res = [];
FTC = struct('done',false);   XC = struct('done',false);
if isempty(BFILE)
    fprintf([newline '[skip] mesh %s mm has no matching B export -- gradient only ' ...
             '(no b.b, no polynomial fit, no cross-checks)' newline], MESH);
else
    bpath = fullfile(cfg.fld_dir, BFILE);
    assert(isfile(bpath), 'xa_grad_fem:nofileB', 'B export not found: %s', bpath);
    rb = import_maxwell_fld(bpath);                % [T], measure frame
    assert(numel(rb.x) == numel(rg.x) && max(abs(rb.x-rg.x)) < 1e-9 ...
           && max(abs(rb.y-rg.y)) < 1e-9 && max(abs(rb.z-rg.z)) < 1e-9, ...
           'xa_grad_fem:grid', 'the two exports are not on the same lattice');
    latB = lat;   latB.V = zeros(size(lat.V));
    latB.V(lat.lin,:) = [rb.bx, rb.by, rb.bz];
    Bq     = trilerp_(latB, Pg(:,1), Pg(:,2), Pg(:,3));
    B_meas = 1e3 * Bq;                             % T -> mT, measure frame
    bb     = sum(B_meas.^2, 2);                    % b.b  [mT^2]  -- frame-invariant scalar
    B_act  = B_meas * RACT.';
    fprintf([newline 'b.b on the same points' newline]);
    fprintf(['  range   %.4f .. %.4f mT^2      (|b| %.4f .. %.4f mT)' newline], ...
            min(bb), max(bb), sqrt(min(bb)), sqrt(max(bb)));
    fprintf(['  at s=0  %.4f mT^2             (|b| %.4f mT)' newline], ...
            interp1(s,bb,0), sqrt(interp1(s,bb,0)));

    % -- least-squares polynomial in the NORMALISED coordinate --------------------
    %  p(xt) = sum_{k=0}^{DEG} a_k xt^k ,  xt = s/SPAN in [-1,1].  Normalising is not
    %  cosmetic: in raw um the degree-7 Vandermonde has cond ~1e15, on [-1,1] ~2e2.
    %  All a_k carry mT^2 because xt is unitless.  res = p - bb is by definition the
    %  LS residual: it averages to ~0 and changes sign.
    xt   = s / SPAN;
    Vp   = xt .^ (0:DEG);
    acof = Vp \ bb;
    p    = Vp * acof;
    res  = p - bb;
    fprintf([newline 'degree-%d least-squares fit of b.b in xt = x_a/%g' newline], DEG, SPAN);
    fprintf(['  cond(V) = %.3e   (same fit in raw um would be %.3e)' newline], ...
            cond(Vp), cond(s.^(0:DEG)));
    for k = 0:DEG, fprintf(['    a%d = %+14.6f mT^2' newline], k, acof(k+1)); end
    fprintf(['  residual p - b.b : %+.4f .. %+.4f mT^2   (max |rel| %.4f %%, rms |rel| %.4f %%)' newline], ...
            min(res), max(res), 100*max(abs(res./bb)), 100*rms(res./bb));

    % -- consistency: the gradient must integrate to the change in b.b ------------
    FTC = struct('done',true, 'I', trapz(s, g_xa), 'D', bb(end) - bb(1));
    FTC.rel_pct = 100*(FTC.I - FTC.D)/FTC.D;
    fprintf([newline '[FTC] integral of d(b.b)/dx_a = %.4f mT^2   vs   b.b(+%g)-b.b(-%g) = %.4f mT^2' newline], ...
            FTC.I, SPAN, SPAN, FTC.D);
    fprintf(['      relative difference %+.3f %%' newline], FTC.rel_pct);

    % -- unit / sign guard: central differences of b.b on the same lattice --------
    %  Kept in SI at the nodes (b.b in T^2, differenced per metre) so the result is
    %  in T^2/m = mT^2/um and lands on the scale of the exported gradient.
    if XCHK
        bbn  = rb.bx.^2 + rb.by.^2 + rb.bz.^2;     % [T^2] at the nodes
        lfd  = lat;   lfd.V = central_diff_(lat, bbn);
        Gf   = trilerp_(lfd, Pg(:,1), Pg(:,2), Pg(:,3));
        g_fd = Gf * xhat;
        rel  = (g_xa - g_fd) ./ max(abs(g_fd), eps);
        XC = struct('done',true, 'g_fd',g_fd, 'rel',rel, ...
                    'max_abs_rel_pct',100*max(abs(rel)), 'median_abs_rel_pct',100*median(abs(rel)));
        fprintf([newline '[XCHK] AEDT Grad(Dot(B,B))  vs  central differences of b.b (same solve)' newline]);
        fprintf(['  fd range  %+.4f .. %+.4f mT^2/um' newline], min(g_fd), max(g_fd));
        fprintf(['  relative difference: median %.2f %%, max %.2f %%' newline], ...
                XC.median_abs_rel_pct, XC.max_abs_rel_pct);
        if XC.median_abs_rel_pct > 25
            warning('xa_grad_fem:xchk', ['the two routes disagree by %.1f %% -- suspect a unit ' ...
                    'or frame error, not just tet discretisation'], XC.median_abs_rel_pct);
        end
    end
end

%% ---- save ---------------------------------------------------------------------
UNIT_BB = 'mT^2';
UNIT_G  = 'mT^2/um  (= T^2/m, numerically identical)';
SRC_G   = gpath;
SRC_B   = '';   if ~isempty(BFILE), SRC_B = fullfile(cfg.fld_dir, BFILE); end
SOLVE   = sprintf('P%d, solver mesh %s mm, export grid 0.02 mm', POLE, MESH);
out = fullfile(DAT, sprintf('xa_fem_%s.mat', TAG));
% [ADDED 2026-09-21] Do not destroy a quantity this run could not produce: the Dot export was
%   overwritten by the smoothed gradient, so bb_pts/res_pts have no source file any more.
%   Carry the stored ones over instead of writing empties over them.
if isempty(bb_pts) && isfile(out)
    Sold = load(out);
    if isfield(Sold,'bb_pts') && ~isempty(Sold.bb_pts)
        bb_pts = Sold.bb_pts;  acof_pts = Sold.acof_pts;  p_pts = Sold.p_pts;
        res_pts = Sold.res_pts;  DCHK = Sold.DCHK;
        fprintf([newline 'bb_pts carried over from the previous %s (its source file is gone)' ...
                 newline], out);
    end
end
save(out, 's','g_xa','G_meas','G_act','bb','B_meas','B_act', ...
     'Pa','Pg','SPAN','NQ','POLE','MODEL','GEOM','MESH','TAG','RACT','xhat', ...
     'xt','DEG','acof','p','res','KINK','g_pts','PCHK', ...
     'bb_pts','acof_pts','p_pts','res_pts','DCHK','g_pts_sm','SCHK', ...
     'B_pts','bmag_pts','s_b','SEGCHK', ...
     'UNIT_BB','UNIT_G','SRC_B','SRC_G','SOLVE','XC','FTC');
fprintf([newline 'wrote %s' newline], out);

%% ---- local: the quantity string in a .fld header --------------------------------
%  AEDT writes it on the "... data "<quantity>"" line, e.g.
%    Vector data "Smooth(Grad(Smooth(Dot(B_Vector, B_Vector))))"
function q = fld_quantity_(fpath)
    fid = fopen(fpath, 'r');
    assert(fid > 0, 'xa_grad_fem:openHdr', 'cannot open %s', fpath);
    q = '';
    for i = 1:8
        ln = fgetl(fid);
        if ~ischar(ln), break, end
        if numel(sscanf(ln, '%f')) >= 4, break, end      % data reached
        tok = regexp(ln, 'data\s*"([^"]*)"', 'tokens', 'once');
        if ~isempty(tok), q = tok{1}; end
    end
    fclose(fid);
end

%% ---- local: read a SCALAR .fld (4 columns: x y z value) -------------------------
%  AEDT writes "Scalar data ..." exports with one value column instead of three, so
%  import_maxwell_fld (fixed at 6 columns) cannot read them.  Header lines are
%  detected the same way: the first line that parses to at least 4 numbers is data.
function [P, v] = read_scalar_fld_(fpath)
    fid = fopen(fpath, 'r');
    assert(fid > 0, 'xa_grad_fem:openScalar', 'cannot open %s', fpath);
    nh = 0;
    for i = 1:8
        ln = fgetl(fid);
        if ~ischar(ln), break, end
        if numel(sscanf(ln, '%f')) >= 4, break, end
        nh = nh + 1;
    end
    frewind(fid);
    C = textscan(fid, '%f%f%f%f', 'HeaderLines', nh, 'CollectOutput', true, ...
                 'MultipleDelimsAsOne', true);
    fclose(fid);
    M = C{1};
    assert(~isempty(M), 'xa_grad_fem:emptyScalar', '%s parsed to nothing (nh=%d)', fpath, nh);
    P = M(:,1:3) * 1e-3;                            % mm -> m
    v = M(:,4);
    fprintf(['  Maxwell scalar .fld: %d points  (%s)' newline], size(M,1), fpath);
end

%% ---- local: lattice builder, reproduced from conv_design_ws.m ------------------
function lat = build_lattice_(x, y, z, V3)
    N  = numel(x);
    q  = @(vv) unique(round(vv*1e9))/1e9;                % quantise to nm
    xg = q(x);   yg = q(y);   zg = q(z);
    n  = [numel(xg) numel(yg) numel(zg)];
    assert(prod(n) == N, 'xa_grad_fem:notGrid', '.fld is not a full regular grid');
    o  = [xg(1) yg(1) zg(1)];
    h  = [(xg(end)-xg(1))/(n(1)-1), (yg(end)-yg(1))/(n(2)-1), (zg(end)-zg(1))/(n(3)-1)];
    ix = round((x-o(1))/h(1))+1;   iy = round((y-o(2))/h(2))+1;   iz = round((z-o(3))/h(3))+1;
    dev = max([max(abs(o(1)+(ix-1)*h(1)-x)), max(abs(o(2)+(iy-1)*h(2)-y)), ...
               max(abs(o(3)+(iz-1)*h(3)-z))]);
    assert(dev < 1e-9, 'xa_grad_fem:notUniform', 'grid is not uniform (max dev %.3e m)', dev);
    lin = sub2ind(n, ix, iy, iz);
    V = zeros(prod(n), 3);
    for cc = 1:3, V(lin,cc) = V3(:,cc);  end
    lat = struct('n',n, 'o',o, 'h',h, 'N_I',1, 'V',V, 'lin',lin);
end

%% ---- local: central differences of a nodal scalar, in lattice node order -------
%  The lattice stores nodes in sub2ind(n,ix,iy,iz) order, so an n-shaped array read
%  with (:) is already in that order and no re-indexing is needed on the way back.
function Gfd = central_diff_(lat, f)
    n = lat.n;   h = lat.h;
    A = zeros(n);   A(lat.lin) = f;
    Gfd = zeros(prod(n),3);
    for d = 1:3
        Bp = permute(A, [d setdiff(1:3,d)]);
        D  = zeros(size(Bp));
        D(2:end-1,:,:) = (Bp(3:end,:,:) - Bp(1:end-2,:,:)) / (2*h(d));
        D(1,:,:)       = (Bp(2,:,:)     - Bp(1,:,:))       / h(d);   % one-sided at the faces
        D(end,:,:)     = (Bp(end,:,:)   - Bp(end-1,:,:))   / h(d);
        D  = ipermute(D, [d setdiff(1:3,d)]);
        Gfd(:,d) = D(:);
    end
end

%% ---- local: trilerp, reproduced from conv_design_ws.m --------------------------
function [Vq, ok, loc] = trilerp_(lat, xq, yq, zq)
    nx = lat.n(1);   ny = lat.n(2);   nz = lat.n(3);
    fx = (xq - lat.o(1)) / lat.h(1);
    fy = (yq - lat.o(2)) / lat.h(2);
    fz = (zq - lat.o(3)) / lat.h(3);
    TOL = 1e-9;
    ok = fx >= -TOL & fx <= nx-1+TOL & fy >= -TOL & fy <= ny-1+TOL ...
       & fz >= -TOL & fz <= nz-1+TOL;
    fx = fx(ok);   fy = fy(ok);   fz = fz(ok);
    i0 = min(max(floor(fx), 0), nx-2);   tx = min(max(fx - i0, 0), 1);
    j0 = min(max(floor(fy), 0), ny-2);   ty = min(max(fy - j0, 0), 1);
    k0 = min(max(floor(fz), 0), nz-2);   tz = min(max(fz - k0, 0), 1);
    a = i0 + 1;   b = j0 + 1;   c = k0 + 1;
    sz = [nx ny nz];
    L  = @(da,db,dc) sub2ind(sz, a+da, b+db, c+dc);
    V  = lat.V;
    V000 = V(L(0,0,0),:);   V100 = V(L(1,0,0),:);
    V010 = V(L(0,1,0),:);   V110 = V(L(1,1,0),:);
    V001 = V(L(0,0,1),:);   V101 = V(L(1,0,1),:);
    V011 = V(L(0,1,1),:);   V111 = V(L(1,1,1),:);
    lp = @(A,Bb,t) A + t .* (Bb - A);
    B00 = lp(V000, V100, tx);   B01 = lp(V001, V101, tx);
    B10 = lp(V010, V110, tx);   B11 = lp(V011, V111, tx);
    C0  = lp(B00,  B10,  ty);   C1  = lp(B01,  B11,  ty);
    Vq  = lp(C0,   C1,   tz);
    loc = struct('i0',i0, 'j0',j0, 'k0',k0, 'tx',tx, 'ty',ty, 'tz',tz);
end
