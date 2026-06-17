%% fit_J_l250.m — Method [J] joint 6-coil fit for l=250 um variant
%  Adapts fit_J.m for l=250 um FEM data (filleted_l250_conv)
%  Initial positions: scale existing l=500 KI_fit.mat positions by 0.5
%  Output: data/J_ideal_fit_l250.mat (used to generate fitting_l250.tex table)

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
cnst = mt_constants();
K_I = eye(6) - ones(6)/6;  % IDEAL K_I
cube_half = 100e-6;

coil_data = struct();
for k=1:6
    dk = import_ansys_data( ...
        fullfile('..','..','results',sprintf('coil%d',k),'filleted_l250_conv'), ...
        'wp', sprintf('coil%d',k));
    mask = abs(dk.x)<cube_half & abs(dk.y)<cube_half & abs(dk.z)<cube_half;
    N_max = 5000;
    if sum(mask) > N_max
        idx_all = find(mask);
        rng(0);
        idx_sub = idx_all(randperm(length(idx_all), N_max));
        mask = false(size(mask));
        mask(idx_sub) = true;
    end
    coil_data(k).px   = dk.x(mask);
    coil_data(k).py   = dk.y(mask);
    coil_data(k).pz   = dk.z(mask);
    coil_data(k).bmag = sqrt(dk.bx(mask).^2 + dk.by(mask).^2 + dk.bz(mask).^2);
    coil_data(k).b_fem = [dk.bx(mask); dk.by(mask); dk.bz(mask)];
    coil_data(k).N    = sum(mask);
    Iv = zeros(6,1); Iv(k) = 1;
    coil_data(k).KI_w = K_I * Iv;
end

% Initial positions: scale l=500 KI_fit.mat pos by 0.5 (since l halved)
KI_data = load('../../data/KI_fit.mat');
pos_init = KI_data.pos * 0.5;

opts = optimset('TolX',1e-10,'TolFun',1e-22,'MaxIter',50000, ...
    'MaxFunEvals',50000,'Display','off');
cost_fn = @(x) cjoint(x, coil_data, 1e-7);

x0s = {pos_init(:)};
rng(42);
x0s{2} = pos_init(:) + 25e-6*randn(18,1);   % 25 um noise (scaled from 50um for l=500)
names = {'per-layer ell x 0.5', '+25um noise'};

fprintf('=== [J] IDEAL K_I on l=250 FEM data ===\n\n');
best_cost = inf; best_pos = []; best_ck = []; best_err = [];
for t = 1:2
    tic;
    [xo, fv] = fminsearch(cost_fn, x0s{t}, opts);
    elapsed = toc;
    [~, ck] = cjoint(xo, coil_data, 1e-7);
    pos = reshape(xo, 3, 6);
    epc = zeros(1, 6);
    for k = 1:6
        Nk = coil_data(k).N;
        pw = [coil_data(k).px, coil_data(k).py, coil_data(k).pz];
        bu = evcf(pos, coil_data(k).KI_w, 1e-7, pw);
        bm = ck(k)*bu;
        res = bm - coil_data(k).b_fem;
        em = sqrt(res(1:Nk).^2 + res(Nk+1:2*Nk).^2 + res(2*Nk+1:3*Nk).^2);
        epc(k) = mean(em ./ coil_data(k).bmag);
    end
    lower = logical(cnst.pole_is_lower);
    fprintf('Trial %d (%s): cost=%.6e, mean err=%.2f%%, %.0fs\n', ...
        t, names{t}, fv, 100*mean(epc), elapsed);
    fprintf('  Per-coil: '); fprintf('%.2f%% ', 100*epc); fprintf('\n');
    fprintf('  |c| Lower=%.1f, Upper=%.1f um\n', ...
        mean(vecnorm(pos(:,lower)))*1e6, mean(vecnorm(pos(:,~lower)))*1e6);
    fprintf('  C_k: '); fprintf('%.4e ', ck); fprintf('\n\n');
    if fv < best_cost
        best_cost = fv;
        best_pos  = pos;
        best_ck   = ck;
        best_err  = epc;
    end
end

%% Also compute ±50 um RMSE error (for table style matching fitting_l500.pdf)
fprintf('\n=== Computing Error over +/-50 um cube (like fitting_l500 table) ===\n');
rmse_cube50 = 0; sum_bmag = 0; n50_total = 0;
for k = 1:6
    dk = import_ansys_data( ...
        fullfile('..','..','results',sprintf('coil%d',k),'filleted_l250_conv'), ...
        'wp', sprintf('coil%d',k));
    mask50 = abs(dk.x)<50e-6 & abs(dk.y)<50e-6 & abs(dk.z)<50e-6;
    pw = [dk.x(mask50), dk.y(mask50), dk.z(mask50)];
    bf = [dk.bx(mask50); dk.by(mask50); dk.bz(mask50)];
    Iv = zeros(6,1); Iv(k) = 1;
    KI_w = K_I * Iv;
    bu = evcf(best_pos, KI_w, 1e-7, pw);
    bm = best_ck(k) * bu;
    res = bm - bf;
    Nk = sum(mask50);
    rmse_cube50 = rmse_cube50 + sum(res.^2);
    sum_bmag = sum_bmag + sum( sqrt(dk.bx(mask50).^2 + dk.by(mask50).^2 + dk.bz(mask50).^2) );
    n50_total = n50_total + Nk;
end
rmse_cube50 = sqrt(rmse_cube50 / (3*n50_total));
avg_bmag50  = sum_bmag / n50_total;
err_cube50  = 100 * rmse_cube50 / avg_bmag50;
fprintf('  +/-50um cube: RMSE=%.4e T, <|B|>=%.4e T, Error=%.3f%%\n', ...
    rmse_cube50, avg_bmag50, err_cube50);

%% Compute per-coil Ra from C_k
mu0 = 4*pi*1e-7;
N_c = 70;
Ra_all = N_c ./ (mu0 * abs(best_ck));
fprintf('\nPer-coil R_a [A/Wb]:\n');
for k = 1:6
    fprintf('  Coil%d (%s): %.4e\n', k, ternary(lower(k),'Lower','Upper'), Ra_all(k));
end
fprintf('  Lower avg: %.4e\n', mean(Ra_all(lower)));
fprintf('  Upper avg: %.4e\n', mean(Ra_all(~lower)));

%% Save
sv.pos          = best_pos;
sv.C_k          = best_ck;
sv.cost         = best_cost;
sv.err_per_coil = best_err;
sv.err_mean     = mean(best_err);
sv.err_cube50   = err_cube50;
sv.Ra_all       = Ra_all;
sv.Ra_lower     = mean(Ra_all(lower));
sv.Ra_upper     = mean(Ra_all(~lower));
sv.cube_half    = cube_half;
sv.method       = 'J-ideal-KI';
sv.K_I          = K_I;
sv.l_um         = 250;
save(fullfile('..','..','data','J_ideal_fit_l250.mat'), '-struct', 'sv');
fprintf('\nSaved: data/J_ideal_fit_l250.mat\n');

%% ---- helpers ----
function [c, Ck] = cjoint(x, cd, km)
    p = reshape(x, 3, 6);
    c = 0;
    Ck = zeros(6, 1);
    for k = 1:6
        pw = [cd(k).px, cd(k).py, cd(k).pz];
        bu = evcf(p, cd(k).KI_w, km, pw);
        Ck(k) = (bu' * cd(k).b_fem) / (bu' * bu);
        r = Ck(k) * bu - cd(k).b_fem;
        c = c + sum(r.^2);
    end
end

function b = evcf(p, w, km, pw)
    N = size(pw, 1);
    bx = zeros(N, 1); by = bx; bz = bx;
    for i = 1:6
        dx = pw(:,1) - p(1,i);
        dy = pw(:,2) - p(2,i);
        dz = pw(:,3) - p(3,i);
        r3 = (dx.^2 + dy.^2 + dz.^2).^(3/2);
        bx = bx + (-w(i))*dx./r3;
        by = by + (-w(i))*dy./r3;
        bz = bz + (-w(i))*dz./r3;
    end
    b = km * [bx; by; bz];
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
