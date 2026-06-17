%% fit_ell_perlayer.m — Per-coil ell sweep + K_I fit with per-layer positions
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
%
%  Step 1: Each coil gets its own optimal ell (sweep + fminbnd)
%  Step 2: Average ell_lower and ell_upper → build per-layer positions
%  Step 3: Linear LS to fit K_I with these positions

cnst = mt_constants();
K_I_ideal = eye(6) - ones(6)/6;
cube_half = 100e-6;

alpha = cnst.alpha;
d_hat = zeros(3,6);
for i = 1:6
    theta = cnst.pole_angles(i)*pi/180;
    zs = sign(cnst.pole_tip_z_wp(i));
    d_hat(:,i) = [cos(theta)*sin(alpha); sin(theta)*sin(alpha); zs*cos(alpha)];
end

fprintf('=== Per-coil ell sweep + K_I fit ===\n\n');

%% Step 1: Per-coil ell sweep
ell_per_coil = zeros(1,6);
C_per_coil = zeros(1,6);
err_per_coil_A = zeros(1,6);

for k = 1:6
    dk = import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'), ...
        'wp', sprintf('coil%d',k));
    mask = abs(dk.x)<cube_half & abs(dk.y)<cube_half & abs(dk.z)<cube_half;
    p_wp = [dk.x(mask), dk.y(mask), dk.z(mask)];
    b_fem = [dk.bx(mask); dk.by(mask); dk.bz(mask)];
    N = sum(mask);

    Iv = zeros(6,1); Iv(k) = 1;
    w_k = K_I_ideal * Iv;

    % Coarse sweep
    ell_scan = linspace(400e-6, 2000e-6, 300);
    cost_scan = arrayfun(@(e) ell_cost(e, d_hat, w_k, cnst.k_m, p_wp, b_fem, N), ell_scan);
    [~, im] = min(cost_scan);

    % Fine search
    ell_opt = fminbnd(@(e) ell_cost(e, d_hat, w_k, cnst.k_m, p_wp, b_fem, N), ...
        max(100e-6, ell_scan(im)-200e-6), ell_scan(im)+200e-6, ...
        optimset('TolX',1e-8,'Display','off'));

    % Evaluate
    [~, C_opt, mean_err] = ell_cost(ell_opt, d_hat, w_k, cnst.k_m, p_wp, b_fem, N);
    ell_per_coil(k) = ell_opt;
    C_per_coil(k) = C_opt;
    err_per_coil_A(k) = mean_err;

    if cnst.pole_is_lower(k), ts='Lower'; else, ts='Upper'; end
    R_a = cnst.N_c / (cnst.mu_0 * abs(C_opt));
    fprintf('Coil%d (P%d %s): ell=%7.1f um, C=%+.4e, R_a=%.3e, err=%.2f%%\n', ...
        k, k, ts, ell_opt*1e6, C_opt, R_a, mean_err);
end

lower_mask = logical(cnst.pole_is_lower);
ell_lower = mean(ell_per_coil(lower_mask));
ell_upper = mean(ell_per_coil(~lower_mask));
fprintf('\nLower avg ell: %.1f um (P1,P3,P6)\n', ell_lower*1e6);
fprintf('Upper avg ell: %.1f um (P2,P4,P5)\n', ell_upper*1e6);

%% Step 2: Build positions with per-layer ell
pos = zeros(3,6);
for i = 1:6
    if cnst.pole_is_lower(i)
        pos(:,i) = ell_lower * d_hat(:,i);
    else
        pos(:,i) = ell_upper * d_hat(:,i);
    end
end

fprintf('\nPositions:\n');
for i = 1:6
    if cnst.pole_is_lower(i), ts='Lower'; else, ts='Upper'; end
    fprintf('  P%d (%s): [%+.1f, %+.1f, %+.1f] um, |c|=%.1f\n', ...
        i, ts, pos(1,i)*1e6, pos(2,i)*1e6, pos(3,i)*1e6, norm(pos(:,i))*1e6);
end

%% Step 3: Linear LS for K_I
d1 = import_ansys_data(fullfile('..','..','results','coil1','filleted'), 'wp', 'coil1');
mask = abs(d1.x)<cube_half & abs(d1.y)<cube_half & abs(d1.z)<cube_half;
N = sum(mask);
px = d1.x(mask); py = d1.y(mask); pz = d1.z(mask);

% Build A matrix
A_mat = zeros(3*N, 6);
for i = 1:6
    dx = px-pos(1,i); dy = py-pos(2,i); dz = pz-pos(3,i);
    r3 = (dx.^2+dy.^2+dz.^2).^(3/2);
    A_mat(1:N,i)       = cnst.k_m * (-dx./r3);
    A_mat(N+1:2*N,i)   = cnst.k_m * (-dy./r3);
    A_mat(2*N+1:3*N,i) = cnst.k_m * (-dz./r3);
end

fprintf('\nA matrix: %d x %d, cond = %.2e\n', size(A_mat,1), size(A_mat,2), cond(A_mat));

% Solve per-coil
Q = zeros(6,6);
err_KI = zeros(1,6);

fprintf('\n--- K_I linear LS ---\n');
fprintf('%-6s %10s\n', 'Coil', 'Vec Err%%');
for k = 1:6
    dk = import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'), ...
        'wp', sprintf('coil%d',k));
    b_fem = [dk.bx(mask); dk.by(mask); dk.bz(mask)];
    q_k = A_mat \ b_fem;
    Q(:,k) = q_k;

    b_model = A_mat * q_k;
    res = b_model - b_fem;
    em = sqrt(res(1:N).^2 + res(N+1:2*N).^2 + res(2*N+1:3*N).^2);
    fm = sqrt(b_fem(1:N).^2 + b_fem(N+1:2*N).^2 + b_fem(2*N+1:3*N).^2);
    err_KI(k) = 100*mean(em./fm);
    fprintf('Coil%d  %9.2f\n', k, err_KI(k));
end
fprintf('Mean:  %9.2f\n', mean(err_KI));

% Extract K_I and C_k
C_k = zeros(6,1);
K_I_fit = zeros(6,6);
for k = 1:6
    C_k(k) = Q(k,k) * 6/5;
    if abs(C_k(k)) > 1e-15
        K_I_fit(:,k) = Q(:,k) / C_k(k);
    end
end

R_a_k = cnst.N_c ./ (cnst.mu_0 * abs(C_k));

fprintf('\n--- Fitted K_I ---\n');
fprintf('         ');
for k=1:6, fprintf('  Coil%d  ',k); end; fprintf('\n');
for i=1:6
    fprintf('P%d  ',i);
    for k=1:6, fprintf(' %+.4f',K_I_fit(i,k)); end
    fprintf('\n');
end

fprintf('\nColumn sums: ');
for k=1:6, fprintf('%+.4f ',sum(K_I_fit(:,k))); end; fprintf('\n');

fprintf('\nC_k:  ');
for k=1:6, fprintf('%+.4e ',C_k(k)); end; fprintf('\n');
fprintf('R_a:  ');
for k=1:6, fprintf('%.3e ',R_a_k(k)); end; fprintf('\n');

% Save
sv.ell_per_coil = ell_per_coil;
sv.ell_lower = ell_lower; sv.ell_upper = ell_upper;
sv.pos = pos; sv.d_hat = d_hat;
sv.Q = Q; sv.K_I_fit = K_I_fit; sv.K_I_ideal = K_I_ideal;
sv.C_k = C_k; sv.R_a_k = R_a_k;
sv.err_per_coil_A = err_per_coil_A;
sv.err_KI = err_KI;
save(fullfile('..','..','data','KI_fit.mat'), '-struct', 'sv');
fprintf('\nSaved to data/KI_fit.mat\n');

%% Local function
function [cost, C, mean_err] = ell_cost(e, dh, w, km, pw, bf, N)
    pos = e * dh;
    bx=zeros(N,1); by=bx; bz=bx;
    for i=1:6
        dx=pw(:,1)-pos(1,i); dy=pw(:,2)-pos(2,i); dz=pw(:,3)-pos(3,i);
        r3=(dx.^2+dy.^2+dz.^2).^(3/2);
        bx=bx+(-w(i))*dx./r3; by=by+(-w(i))*dy./r3; bz=bz+(-w(i))*dz./r3;
    end
    bu = km*[bx;by;bz];
    C = (bu'*bf)/(bu'*bu);
    r = C*bu - bf;
    cost = sum(r.^2);
    if nargout > 2
        em = sqrt(r(1:N).^2+r(N+1:2*N).^2+r(2*N+1:3*N).^2);
        fm = sqrt(bf(1:N).^2+bf(N+1:2*N).^2+bf(2*N+1:3*N).^2);
        mean_err = 100*mean(em./fm);
    end
end
