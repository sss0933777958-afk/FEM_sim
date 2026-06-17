%% fit_KI_v1.m — Fit real K_I matrix from 6 single-coil FEM data
% VARIANT VERSION: see variants/README or analysis/README for context
%
%  Fix charge positions at ell * d_hat (sphere from Method [A])
%  For each coil k, solve linear least squares: b_k = A * q_k
%  q_ki = C_k * w_ki → extract K_I and C_k
%
%  No optimizer needed — pure linear algebra.

%% 1. Setup
cnst = mt_constants();
cube_half = 200e-6;

fprintf('=== Fit K_I: linear least squares ===\n\n');

% Load ell from Method [A]
A_fit = load('../../data/charge_model_fit.mat');
ell = A_fit.ell_opt;
fprintf('ell (from [A]): %.1f um\n', ell*1e6);

% Build d_hat and positions
alpha = cnst.alpha;
d_hat = zeros(3, 6);
for i = 1:6
    theta = cnst.pole_angles(i) * pi/180;
    z_sign = sign(cnst.pole_tip_z_wp(i));
    d_hat(:,i) = [cos(theta)*sin(alpha); sin(theta)*sin(alpha); z_sign*cos(alpha)];
end
pos = ell * d_hat;

fprintf('Charge positions (sphere):\n');
for i = 1:6
    fprintf('  P%d: [%+.1f, %+.1f, %+.1f] um, |c|=%.1f um\n', ...
        i, pos(1,i)*1e6, pos(2,i)*1e6, pos(3,i)*1e6, norm(pos(:,i))*1e6);
end

%% 2. Build basis matrix A (same for all coils, depends only on positions)
% Load first coil to get node coordinates
d1 = import_ansys_data(fullfile('..','..','results','coil1','filleted'), 'wp', 'coil1');
mask = abs(d1.x) < cube_half & abs(d1.y) < cube_half & abs(d1.z) < cube_half;
N = sum(mask);
px = d1.x(mask); py = d1.y(mask); pz = d1.z(mask);
fprintf('\nFitting cube: +/-%.0f um, %d nodes\n', cube_half*1e6, N);

% A matrix: 3N x 6
% Column i = k_m * [-(p - c_i) / |p - c_i|^3] for all N points, xyz stacked
A_mat = zeros(3*N, 6);
for i = 1:6
    dx = px - pos(1,i);
    dy = py - pos(2,i);
    dz = pz - pos(3,i);
    r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
    A_mat(1:N, i)       = cnst.k_m * (-dx ./ r3);
    A_mat(N+1:2*N, i)   = cnst.k_m * (-dy ./ r3);
    A_mat(2*N+1:3*N, i) = cnst.k_m * (-dz ./ r3);
end

%% 3. Solve for Q matrix (6x6): each column from one coil
Q = zeros(6, 6);  % Q(:,k) = q_k = C_k * w_k
residuals = zeros(1, 6);
err_per_coil = zeros(1, 6);

fprintf('\n--- Per-coil linear least squares ---\n');
fprintf('%-6s %10s %10s\n', 'Coil', 'Residual', 'Vec Err%%');

for k = 1:6
    coilname = sprintf('coil%d', k);
    dk = import_ansys_data(fullfile('..','..','results',coilname,'filleted'), 'wp', coilname);
    b_fem = [dk.bx(mask); dk.by(mask); dk.bz(mask)];

    % Linear least squares: A * q = b → q = A \ b
    q_k = A_mat \ b_fem;
    Q(:, k) = q_k;

    % Residual and error
    b_model = A_mat * q_k;
    res = b_model - b_fem;
    residuals(k) = norm(res);

    err_mag = sqrt(res(1:N).^2 + res(N+1:2*N).^2 + res(2*N+1:3*N).^2);
    fem_mag = sqrt(b_fem(1:N).^2 + b_fem(N+1:2*N).^2 + b_fem(2*N+1:3*N).^2);
    err_per_coil(k) = 100 * mean(err_mag ./ fem_mag);

    fprintf('Coil%d  %10.4e  %9.2f\n', k, residuals(k), err_per_coil(k));
end
fprintf('Mean:              %9.2f\n', mean(err_per_coil));

%% 4. Extract K_I and C_k from Q
% Q(:,k) = C_k * K_I(:,k)
% Normalization: K_I diagonal should be the dominant element
% C_k = Q(k,k) / K_I_ideal(k,k) = Q(k,k) / (5/6)
% Or simpler: C_k = sum of |Q(:,k)| with sign from diagonal

fprintf('\n--- Extracting K_I and C_k ---\n');
C_k = zeros(6, 1);
K_I_fit = zeros(6, 6);

for k = 1:6
    % C_k determined by: K_I diagonal = dominant → C_k = Q(k,k) / (5/6)
    % But we don't know the ideal diagonal. Use: C_k = norm of q_k / norm of ideal w_k
    % Simplest: C_k = Q(k,k) * 6/5 (assuming diagonal ≈ 5/6)
    C_k(k) = Q(k,k) * 6/5;
    if abs(C_k(k)) > 1e-15
        K_I_fit(:,k) = Q(:,k) / C_k(k);
    end
end

R_a_k = cnst.N_c ./ (cnst.mu_0 * abs(C_k));

fprintf('\nC_k and R_a:\n');
for k = 1:6
    if cnst.pole_is_lower(k), ts='Lower'; else, ts='Upper'; end
    fprintf('  Coil%d (%s): C_k = %+.4e, R_a = %.3e\n', k, ts, C_k(k), R_a_k(k));
end

fprintf('\n--- Fitted K_I matrix ---\n');
fprintf('         ');
for k = 1:6, fprintf('  Coil%d  ', k); end
fprintf('\n');
for i = 1:6
    fprintf('P%d  ', i);
    for k = 1:6
        fprintf(' %+.4f', K_I_fit(i,k));
    end
    fprintf('\n');
end

fprintf('\nColumn sums (should be ~0 for flux conservation):\n  ');
for k = 1:6, fprintf('%+.4f ', sum(K_I_fit(:,k))); end
fprintf('\n');

% Compare with ideal K_I
K_I_ideal = eye(6) - ones(6)/6;
K_I_diff = K_I_fit - K_I_ideal;
fprintf('\n--- K_I_fit - K_I_ideal (deviation from ideal) ---\n');
fprintf('         ');
for k = 1:6, fprintf('  Coil%d  ', k); end
fprintf('\n');
for i = 1:6
    fprintf('P%d  ', i);
    for k = 1:6
        fprintf(' %+.4f', K_I_diff(i,k));
    end
    fprintf('\n');
end
fprintf('\nMax |deviation|: %.4f\n', max(abs(K_I_diff(:))));

%% 5. Save
sv.Q = Q;
sv.K_I_fit = K_I_fit;
sv.K_I_ideal = K_I_ideal;
sv.C_k = C_k;
sv.R_a_k = R_a_k;
sv.pos = pos;
sv.ell = ell;
sv.d_hat = d_hat;
sv.cube_half = cube_half;
sv.N_nodes = N;
sv.err_per_coil = err_per_coil;
sv.err_mean = mean(err_per_coil);
save(fullfile('..','..','data','variants','KI_fit_v1.mat'), '-struct', 'sv');
fprintf('\nSaved to data/variants/KI_fit_v1.mat\n');
