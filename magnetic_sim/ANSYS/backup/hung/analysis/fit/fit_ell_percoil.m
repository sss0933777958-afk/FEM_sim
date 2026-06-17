%% fit_ell_percoil.m — Per-coil single-charge ell sweep (no K_I)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
%  Each coil: fit only the excited pole's charge position
%  c_k = ell_k * d_hat_k, C_k via VarPro

cnst = mt_constants();
cube_half = 100e-6;

alpha = cnst.alpha;
d_hat = zeros(3,6);
for i = 1:6
    theta = cnst.pole_angles(i)*pi/180;
    zs = sign(cnst.pole_tip_z_wp(i));
    d_hat(:,i) = [cos(theta)*sin(alpha); sin(theta)*sin(alpha); zs*cos(alpha)];
end

fprintf('=== Single-charge per-coil ell sweep (no K_I) ===\n\n');

ell_vals = zeros(1,6);
C_vals = zeros(1,6);
err_vals = zeros(1,6);

for k = 1:6
    dk = import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'), ...
        'wp', sprintf('coil%d',k));
    mask = abs(dk.x)<cube_half & abs(dk.y)<cube_half & abs(dk.z)<cube_half;
    px = dk.x(mask); py = dk.y(mask); pz = dk.z(mask);
    bfem = [dk.bx(mask); dk.by(mask); dk.bz(mask)];
    N = sum(mask);

    % Sweep ell
    ell_scan = linspace(300e-6, 3000e-6, 500);
    cost_scan = arrayfun(@(e) sc_cost(e, d_hat(:,k), px, py, pz, bfem), ell_scan);
    [~, im] = min(cost_scan);

    ell_opt = fminbnd(@(e) sc_cost(e, d_hat(:,k), px, py, pz, bfem), ...
        max(100e-6, ell_scan(im)-200e-6), min(3000e-6, ell_scan(im)+200e-6), ...
        optimset('TolX',1e-8,'Display','off'));

    % Evaluate at optimum
    [~, Ct, me] = sc_cost(ell_opt, d_hat(:,k), px, py, pz, bfem);

    ell_vals(k) = ell_opt;
    C_vals(k) = Ct;
    err_vals(k) = me;

    if cnst.pole_is_lower(k), ts='Lower'; else, ts='Upper'; end
    Ra = cnst.N_c / (cnst.mu_0 * abs(Ct));
    fprintf('Coil%d (P%d %s): ell=%7.1f um, C=%+.4e, R_a=%.3e, err=%.2f%%\n', ...
        k, k, ts, ell_opt*1e6, Ct, Ra, me);
end

fprintf('\n--- Summary ---\n');
lower = logical(cnst.pole_is_lower);
fprintf('Lower (P1,P3,P6) avg ell: %.1f um\n', mean(ell_vals(lower))*1e6);
fprintf('Upper (P2,P4,P5) avg ell: %.1f um\n', mean(ell_vals(~lower))*1e6);
fprintf('Lower avg |C|: %.4e\n', mean(abs(C_vals(lower))));
fprintf('Upper avg |C|: %.4e\n', mean(abs(C_vals(~lower))));
fprintf('Mean error: %.2f%%\n', mean(err_vals));

save(fullfile('..','..','data','single_charge_ell.mat'), ...
    'ell_vals','C_vals','err_vals','d_hat','cube_half');
fprintf('\nSaved to data/single_charge_ell.mat\n');

function [cost, C, mean_err] = sc_cost(e, dh, px, py, pz, bf)
    c = e * dh;
    dx = px-c(1); dy = py-c(2); dz = pz-c(3);
    r3 = (dx.^2+dy.^2+dz.^2).^(3/2);
    N = length(px);
    bu = 1e-7 * [-dx./r3; -dy./r3; -dz./r3];
    C = (bu'*bf) / (bu'*bu);
    res = C*bu - bf;
    cost = sum(res.^2);
    if nargout > 2
        em = sqrt(res(1:N).^2 + res(N+1:2*N).^2 + res(2*N+1:3*N).^2);
        fm = sqrt(bf(1:N).^2 + bf(N+1:2*N).^2 + bf(2*N+1:3*N).^2);
        mean_err = 100*mean(em./fm);
    end
end
