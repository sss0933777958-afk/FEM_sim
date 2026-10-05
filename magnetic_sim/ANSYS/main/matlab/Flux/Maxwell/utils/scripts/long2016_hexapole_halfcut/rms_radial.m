%% rms_radial.m -- why does the eighteen-parameter RMS converge more slowly (relatively)?
%  Hypothesis: RMS(Nc) falls with Nc mainly because of SAMPLING WEIGHT. The axes-shell ladder puts
%  6 of its 6*Nr+1 points on every shell, so few shells over-weight the outer radius, where the
%  residual is largest. A residual concentrated more strongly at the outer radius converges more
%  slowly in relative terms.
%  Test: fit both models once on a fine ladder (fixed parameters), get the residual-squared profile
%  e2(r), then predict eps(Nc) from the shell weights alone and compare with rms_ladder.mat.
%  Output: utils/data/<MODEL>/rms_radial.mat

MODEL = 'long2016_hexapole_halfcut';   GEOM = 'tip40um';   VARIANT = 'maxwell';
R = 150e-6;   l0 = 0.5e-3;   NRF = 200;                      % fine ladder for the profile
NC = [523 667];                                              % converged points (ten-level criterion)

here = fileparts(mfilename('fullpath'));   CAL = fileparts(fileparts(fileparts(here)));
addpath(fullfile(CAL, 'function'), fullfile(CAL, 'common_path'));
cfg = model_config(MODEL, GEOM);
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
opt = struct('model',MODEL, 'geom',GEOM, 'variant',VARIANT, 'frame','actuator', 'raw',raw);
[P, Bs] = conv_design_ws(NRF, R, opt);
r  = vecnorm(P, 2, 2);                                       % point radius [m]
rs = (0:NRF).' * R / NRF;                                    % shell radii (0 = centre)
sh = round(r / (R / NRF));                                   % shell index of every point

L  = load(fullfile(CAL, 'utils', 'data', MODEL, 'rms_ladder.mat'));
e2s = nan(NRF+1, 2);   q = nan(1,2);   outer = nan(1,2);   eps_pred = nan(numel(L.Nr), 2);
for m = 1:2
    [e, l_hat, J] = fitting(P, Bs, ad.Pc_base, l0, m == 2);
    S   = kernel(l_hat, pc(e, ad.Pc_base), P);
    res = S * ((S.' * S) \ (S.' * Bs)) - Bs;                 % 3Np x 6
    assert(abs(sum(res(:).^2) / J - 1) < 1e-8, 'residual does not reproduce J');
    e2  = sum(reshape(sum(res.^2, 2), 3, []), 1).';          % per point: 3 components x 6 excitations
    for k = 0:NRF, e2s(k+1, m) = mean(e2(sh == k)); end
    in  = rs > 0.2*R;                                        % power law e2 ~ r^q on the outer 80 %
    pq  = polyfit(log(rs(in)), log(e2s(in, m)), 1);   q(m) = pq(1);
    w   = rs.^2;                                             % share of the continuous mean from r > 0.8 R
    outer(m) = sum(e2s(rs > 0.8*R, m)) / sum(e2s(:, m));
    % prediction from shell weights alone (parameters frozen): RMS^2(Nr) ~ [e2(0) + 6*sum e2(R k/Nr)]/(6Nr+1)
    f = @(x) interp1(rs, e2s(:, m), x, 'linear');
    rp = arrayfun(@(nr) sqrt((f(0) + 6*sum(f((1:nr)*R/nr))) / (6*nr + 1)), L.Nr);
    eps_pred(:, m) = rp / rp(L.N == NC(m)) - 1;
end
eps_act = L.RMS ./ [L.RMS(L.N == NC(1), 1), L.RMS(L.N == NC(2), 2)] - 1;

for n = [13 25 43 55 103]
    k = find(L.N == n);
    fprintf('N=%3d | single eps actual %6.2f%% predicted %6.2f%% | eighteen actual %6.2f%% predicted %6.2f%%\n', ...
            n, 100*eps_act(k,1), 100*eps_pred(k,1), 100*eps_act(k,2), 100*eps_pred(k,2));
end
fprintf('radial power law e2 ~ r^q : single q = %.2f, eighteen q = %.2f\n', q);
fprintf('share of residual^2 (shell sum) at r > 0.8R : single %.1f%%, eighteen %.1f%%\n', 100*outer);
fprintf('centre residual^2 / outer-shell residual^2 : single %.4f, eighteen %.4f\n', ...
        e2s(1,1)/e2s(end,1), e2s(1,2)/e2s(end,2));

out = fullfile(CAL, 'utils', 'data', MODEL, 'rms_radial.mat');
save(out, 'rs', 'e2s', 'q', 'outer', 'eps_pred', 'eps_act', 'NRF', 'NC');
fprintf('saved %s\n', out);

% ---- same charge grid and kernel as function/fitting.m ----------------------
function Pc = pc(e17, Pc_base)
    E = zeros(3, 6);
    E(:,1) = e17(1:3);   E(:,2) = e17(4:6);   E(:,3) = e17(7:9);   E(:,4) = e17(10:12);
    E(:,5) = e17(13:15); E(1,6) = e17(16);    E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);
    Pc = Pc_base + E;
end
function S = kernel(l_hat, Pc, P)
    Np = size(P, 1);   pbar = P / l_hat;   S = zeros(3*Np, 6);
    for k = 1:6
        d = pbar - Pc(:,k).';   r3 = sum(d.^2, 2).^1.5;
        S(:,k) = reshape((d ./ r3).', 3*Np, 1);
    end
end
