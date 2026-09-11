function [Mbar_hat, gF_hat, gB_hat, out] = solve_force(H_hat, mgB, l_hat, UF, base)
%SOLVE_FORCE  Gauge the fitted force transfer matrix, then convert to B g.
%
%   [Mbar_hat, gF_hat, gB_hat, out] = SOLVE_FORCE(H_hat, mgB, l_hat, UF, base)
%     H_hat  : 6 x 6 fitted force transfer matrix (F H_I hat or F H_V hat)
%     mgB    : particle magnetisation constant [A.um^2/mT], from  m = mgB * b
%     l_hat  : fitted effective length [um]
%     UF     : force unit factor, default 1e3  (see UNITS below)
%     base   : 'current' (default) | 'voltage' -- LABELS AND UNITS ONLY
%
%     Mbar_hat : 6 x 6 gauged matrix, Mbar(1,1) = 5/6
%                  current -> K_I bar        voltage -> D bar
%     gF_hat   : force gain  F g   <- WHAT THE GAUGE COMPUTES
%                  current -> F g_I [pN/A^2] voltage -> F g_V [pN/mV^2]
%     gB_hat   : field gain  B g   <- DERIVED from gF_hat
%                  current -> B g_I [mT/A]   voltage -> B g_V [mT/mV]
%     out      : h11 / rt / gF / gB / scale / flipped columns / labels / checks
%
%   [RENAMED 2026-09-02] from solve_force_current, and made base-generic. The
%   ARITHMETIC IS IDENTICAL for the two bases -- every line below is unchanged
%   between current and voltage, including the scale sqrt(UF*mgB/(2*l)), which
%   contains neither a current nor a voltage. Only the printed names and the
%   units differ, which is all `base` selects.
%
%   THE GAUGE.  The fit only ever sees H_hat. Gauging it with Mbar(1,1) = 5/6:
%
%       H_hat  = sqrt(F g) * Mbar,   Mbar(1,1) = 5/6
%           =>  (5/6) * sqrt(F g) = h11,   h11 = H_hat(1,1)
%           =>  sqrt(F g) = 6*h11/5
%           =>  F g       = (6*h11/5)^2
%           =>  Mbar      = 5/(6*h11) * H_hat
%
%   B g is NOT an independent result of the fit. It follows from F g through
%
%       F g  = UF * mgB/(2*l) * (B g)^2
%       B g  = (6*h11/5) / sqrt(UF*mgB/(2*l))
%
%   and is only recoverable because mgB is supplied as a known constant. The
%   force data alone fixes the product mgB * (B g)^2, never either factor.
%
%   UNITS.  With lengths in um, fields in mT, excitations in A (current) or mV
%   (voltage) and mgB in A.um^2/mT, the force comes out in A.um.mT. A.m.T = N
%   exactly, so A.um.mT = 1e-9 N = 1 nN = 1e3 pN: UF = 1e3 converts to pN, the
%   project's force unit. It must be the SAME factor that was folded into c
%   when the measurement was generated. The scale sqrt(UF*mgB/(2*l)) carries
%   sqrt(pN)/mT in both bases.
%
%   COLUMN SIGNS.  f = (H u)' L (H u) is invariant under a sign flip of any
%   column of H_hat, so with one-at-a-time excitation the optimiser can land on
%   any of the 2^6 equivalent solutions. The signs are fixed here by the
%   physical requirement that self-excitation is positive, diag(Mbar) > 0 --
%   flipping column j flips Mbar(j,j) and nothing else about the fit quality.
%   Doing this before the gauge also guarantees h11 > 0, hence F g > 0, B g > 0.
%
%   The remaining structure of Mbar is NOT imposed; out.chk reports it as an
%   independent check. Read it per base: diag_pos and diag_dominant are
%   meaningful for both, but offdiag_neg and a near-zero rowsum_max are charge
%   neutrality properties of K_I bar and are NOT expected of D bar -- a voltage
%   run legitimately reports offdiag_neg = 0 and a rowsum of order 1. Those two
%   are inherited from whatever Mbar the flux calibration produced, so they say
%   something about that calibration, never about this gauge.
%
%   See also FITTING_FORCE, BUILD_L, BASE_MODEL.

    if nargin < 4 || isempty(UF),   UF   = 1e3;       end
    if nargin < 5 || isempty(base), base = 'current'; end
    validateattributes(H_hat, {'numeric'}, {'size',[6 6],'real','finite'}, mfilename, 'H_hat');
    validateattributes(mgB,   {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'mgB');
    validateattributes(l_hat, {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'l_hat');
    validateattributes(UF,    {'numeric'}, {'scalar','real','positive','finite'}, mfilename, 'UF');
    base = validatestring(base, {'current','voltage'}, mfilename, 'base');

    % ---- resolve the per-column sign ambiguity ---------------------------
    dsign   = sign(diag(H_hat)).';
    dsign(dsign == 0) = 1;                       % degenerate, leave untouched
    flipped = find(dsign < 0);
    H_fix   = H_hat .* dsign;                    % scales column j by dsign(j)

    % ---- gauge: this is where F g and Mbar come from ---------------------
    h11 = H_fix(1,1);
    assert(h11 > 0, 'solve_force:h11', ...
           'H_hat(1,1) = %g is not positive after the sign fix', h11);
    Mbar_hat = (5 / (6 * h11)) * H_fix;
    rt       = 6 * h11 / 5;                      % sqrt(F g), the factor on Mbar
    gF_hat   = rt^2;                             % F g

    % ---- derive B g from F g (needs mgB as a known constant) -------------
    sc     = sqrt(UF * mgB / (2 * l_hat));       % sqrt(F g) = sc * B g
    gB_hat = rt / sc;                            % B g

    % ---- physical structure checks (reported, not imposed) ---------------
    od         = Mbar_hat(~eye(6));
    [~, amax]  = max(abs(Mbar_hat), [], 2);
    chk = struct( ...
        'diag_pos',      all(diag(Mbar_hat) > 0), ...
        'diag_dominant', isequal(amax(:).', 1:6), ...
        'offdiag_neg',   all(od < 0), ...
        'rowsum_max',    max(abs(sum(Mbar_hat, 2))), ...
        'M11',           Mbar_hat(1,1));

    out = struct('h11', h11, 'rt', rt, 'gF', gF_hat, 'gB', gB_hat, 'scale', sc, ...
                 'UF', UF, 'base', base, 'labels', force_labels(base), ...
                 'flipped', flipped, 'nflip', numel(flipped), 'chk', chk);
end

% ---- the only base-dependent thing in this file -----------------------------
function lb = force_labels(base)
%FORCE_LABELS  Display names and units for one base. Plain names for the
%   console, LaTeX (no leading backslash escaping needed by the caller) for the
%   PDF. Kept here so solve_force stays the single source of the naming.
%   The *_m variants carry the measurement-side subscript; they are stored
%   rather than composed by the caller because appending ",m" has to go inside
%   the existing subscript braces, not after them.
    if strcmp(base, 'current')
        lb = struct( ...
            'u_name',   'I',            'u_unit',    'A', ...
            'u_tex',    'F',            'M_name',    'K_I bar', ...
            'M_tex',    '\bar{K}_{I}',  'M_tex_m',   '\bar{K}_{I,m}', ...
            'gF_name',  'F g_I',        'gF_tex',    '{}^{F}\hat{g}_{I}', ...
            'gF_tex_m', '{}^{F}\hat{g}_{I,m}', ...
            'gF_unit',  'pN/A^2',       'gF_texu',   'pN/A^{2}', ...
            'gB_name',  'B g_I',        'gB_tex',    '{}^{B}\hat{g}_{I}', ...
            'gB_tex_m', '{}^{B}\hat{g}_{I,m}', ...
            'gB_unit',  'mT/A',         'gB_texu',   'mT/A', ...
            'H_tex',    '{}^{F}\hat{H}_{I}');
    else
        lb = struct( ...
            'u_name',   'V',            'u_unit',    'mV', ...
            'u_tex',    'V',            'M_name',    'D bar', ...
            'M_tex',    '\bar{D}',      'M_tex_m',   '\bar{D}_{m}', ...
            'gF_name',  'F g_V',        'gF_tex',    '{}^{F}\hat{g}_{V}', ...
            'gF_tex_m', '{}^{F}\hat{g}_{V,m}', ...
            'gF_unit',  'pN/mV^2',      'gF_texu',   'pN/mV^{2}', ...
            'gB_name',  'B g_V',        'gB_tex',    '{}^{B}\hat{g}_{V}', ...
            'gB_tex_m', '{}^{B}\hat{g}_{V,m}', ...
            'gB_unit',  'mT/mV',        'gB_texu',   'mT/mV', ...
            'H_tex',    '{}^{F}\hat{H}_{V}');
    end
end
