function [P, info] = sample_axes_shells(R, Nr, R_act, opt)
%SAMPLE_AXES_SHELLS  Concentric shells of the six actuator-axis directions.
% =========================================================================
%   The design the user specified (2026-08-27): put the six pole directions
%   on Nr equally spaced shells, plus the centre.
%
%       directions (ACTUATOR frame)   (+1,0,0) (-1,0,0)   <- u axis, P1 / P2
%                                     (0,+1,0) (0,-1,0)   <- v axis, P3 / P4
%                                     (0,0,+1) (0,0,-1)   <- w axis, P5 / P6
%       shell radii                   r_k = (R/Nr)*k ,  k = 1 ... Nr
%       point count                   6*Nr + 1   (the +1 is the centre)
%
%       Nr = 1  ->  7      Nr = 3  ->  19     Nr = 10 ->  61
%       Nr = 2  -> 13      Nr = 4  ->  25     Nr = 166 -> 997
%
%   The outermost shell (k = Nr) sits exactly ON the sampling sphere r = R.
%
%   ── Why this layout ─────────────────────────────────────────────────
%     EVERY POLE GETS ITS OWN POINT, on every shell, and the +/- sides are
%     treated equally. That is what the 4-point designs cannot do at once:
%     taking +/- pairs of only two axes is exactly singular (rank 5, the
%     unused axis' charge pair is mirror-symmetric about the sample plane),
%     while taking only the + side of three axes is full rank but never
%     samples P2 / P4 / P6 -> their charge strengths stay poorly determined
%     and K_I_bar loses its physical structure.
%
%     Radial spacing also maximises the leverage on l_hat: moving along a
%     pole axis is the direction in which the distance to that charge
%     changes fastest, so the field falloff -- the only thing that carries
%     the length scale -- is sampled head on. (At the centre, by contrast,
%     S is exactly independent of l_hat.)
%
%   ── Known weakness ──────────────────────────────────────────────────
%     All points lie on three straight lines; there is no off-axis angular
%     coverage. Charge strengths and l_hat should come out clean, but the
%     components of e TRANSVERSE to each pole axis have little leverage
%     (moving a charge sideways is a second-order effect on the field along
%     its own axis). Expect single-parameter to do well and some eighteen-
%     parameter components to stay soft.
%
%   [MOVED 2026-08-30 user decision] This file now lives in
%   matlab/Flux/Maxwell/function/ (it used to sit in temp_code/scripts/).
%   It is the CANONICAL workspace sampler: main.m drives it through the
%   WS_AXSH_NR knob.  Contract -- design knobs in, point positions out:
%
%       in   Nr (shell count) and R (sampling radius);
%            R_act only fixes the frame, it is not a design parameter
%            (axsh_vs_R.m deliberately feeds a ROTATED R_act for its
%             rotation-invariance test, so it stays an explicit argument)
%       out  P (Np x 3, measure frame) -- positions only, NO field.
%            Fetching the field is conv_design_ws's job (opt.query mode).
%
%   Output is in the MEASURE frame (what conv_design_ws expects in
%   opt.query): P_meas = R_act' * P_act, since R_act's rows are the
%   actuator basis vectors written in measure coordinates.
%
%   Usage
%     cfg = model_config('long2016_hexapole_halfcut','tip40um');
%     P = sample_axes_shells(150e-6, 3, cfg.R_act);                  % 19 points
%     P = sample_axes_shells(150e-6, 3, cfg.R_act, struct('quiet',true));
%     P = sample_axes_shells(150e-6, 3, cfg.R_act, struct('center',false));  % 18
% =========================================================================
    if nargin < 1 || isempty(R),     R   = 1;        end
    if nargin < 2 || isempty(Nr),    Nr  = 1;        end
    if nargin < 3 || isempty(R_act), R_act = eye(3); end
    if nargin < 4 || isempty(opt),   opt = struct(); end
    QUIET = getdef_(opt, 'quiet',  false);
    ADDC  = getdef_(opt, 'center', true);      % include the r = 0 point
    assert(Nr >= 1 && Nr == round(Nr), 'Nr must be a positive integer');

    % six actuator-axis directions, in the charge order P1..P6 of Pc_base
    D = [ 1  0  0
         -1  0  0
          0  1  0
          0 -1  0
          0  0  1
          0  0 -1];

    h  = R / Nr;                               % shell spacing
    Pa = zeros(6*Nr, 3);
    for k = 1:Nr
        Pa((k-1)*6+1 : k*6, :) = D * (h*k);
    end
    if ADDC, Pa = [zeros(1,3); Pa]; end        % centre first, then shells outward

    P = Pa * R_act;                            % actuator -> measure frame

    info = struct('Nr',Nr, 'R',R, 'h',h, 'npts',size(P,1), ...
                  'radii',(1:Nr)*h, 'center',ADDC, 'P_act',Pa);

    if ~QUIET
        fprintf(['sample_axes_shells: Nr=%d -> %d points ' ...
                 '(6 axis directions x %d shells%s), r = %.1f ... %.1f um\n'], ...
                Nr, size(P,1), Nr, ternary_(ADDC,' + centre',''), h*1e6, R*1e6);
    end
end

% ----------------------------------------------------------------------------
function v = getdef_(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function v = ternary_(c, a, b)
    if c, v = a; else, v = b; end
end
