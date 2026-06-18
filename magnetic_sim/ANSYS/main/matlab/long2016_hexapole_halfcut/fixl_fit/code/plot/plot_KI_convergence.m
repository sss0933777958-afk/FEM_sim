%% plot_KI_convergence.m — lsqnonlin convergence (cost%% vs iteration) per radius
%  Same K_I model/data as sweep_KI_radius.m, but gB initial guess = 100.
%  Records relative-RMS cost (%) at every iteration via an OutputFcn and overlays
%  the 10 ball radii (50..500 um) on one semilogy plot.

clear; clc; close all;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
fig_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];
I_actual = 0.6;
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

% load 6 coils once
C = struct('P',{},'Bn',{},'pj',{});
for k = 1:6
    cname = sprintf('coil%d', k);
    d   = import_ansys_data(fullfile(results_root, cname), 'wp', cname);
    air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize', false));
    zwp = d.z - cnst.SPH_OFST;
    C(k).P  = [d.x(air), d.y(air), zwp(air)];
    C(k).Bn = -[d.bx(air), d.by(air), d.bz(air)];
    C(k).pj = apdl_to_paper_idx(k);
end

R_list = (50:50:500)*1e-6;  Nmax = 4000;
Khat0  = eye(6) - ones(6)/6;  ell0 = 0.5e-3;
gB0    = 100;                              % <-- initial guess changed to 100
freemask = true(6); freemask(1,1) = false;
NIT    = 25;                               % fixed iteration count for every radius

hist = cell(1,numel(R_list));
for ri = 1:numel(R_list)
    R = R_list(ri);
    coil = struct('p',{},'bfem',{},'pj',{});
    for k = 1:6
        idx = find(sum(C(k).P.^2,2) < R^2);
        rng(0); if numel(idx)>Nmax, idx = idx(randperm(numel(idx),Nmax)); end
        Bn = C(k).Bn(idx,:);
        coil(k).p = C(k).P(idx,:);
        coil(k).bfem = [Bn(:,1);Bn(:,2);Bn(:,3)];
        coil(k).pj = C(k).pj;
    end
    x0 = [ell0*1e3; gB0; Khat0(freemask)];
    hist{ri} = fit_with_history(coil, dhat, I_actual, freemask, x0, NIT);
    fprintf('R=%3d um : %d iters, final cost = %.4e\n', ...
            round(R*1e6), numel(hist{ri})-1, hist{ri}(end));
end

%% plot
figure('Color','w','Position',[100 100 760 560]); hold on;
cmap = turbo(numel(R_list));
h = gobjects(1,numel(R_list));
for ri = 1:numel(R_list)
    h(ri) = semilogy(0:numel(hist{ri})-1, hist{ri}, '-o', 'Color', cmap(ri,:), ...
                     'LineWidth', 1.4, 'MarkerSize', 3, 'MarkerFaceColor', cmap(ri,:));
end
set(gca,'YScale','log'); grid on; box on;
xlabel('iteration'); ylabel('cost');
title('charge-model fit convergence');
legend(h, compose('R = %d \\mum', round(R_list*1e6)), 'Location','northeast', 'NumColumns',2);
set(gca,'FontSize',11);

outpng = fullfile(fig_dir, 'KI_cost_convergence_gB100.png');
exportgraphics(gcf, outpng, 'Resolution', 150);
fprintf('\nwrote %s\n', outpng);


%% ===== local functions =====
function hist = fit_with_history(coil, dhat, I, freemask, x0, NIT)
    hist = [];
    % tolerances = 0 + fixed MaxIterations -> every radius runs exactly NIT iters
    opts = optimoptions('lsqnonlin','Display','off', ...
        'MaxFunctionEvaluations',1e7,'MaxIterations',NIT, ...
        'FunctionTolerance',0,'StepTolerance',0,'OptimalityTolerance',0, ...
        'OutputFcn', @outfun);
    lsqnonlin(@(x) resid_all(x, coil, dhat, I, freemask), x0, [], [], opts);
    function stop = outfun(~, ov, state)
        stop = false;
        if strcmp(state,'iter')
            hist(end+1) = ov.resnorm;   % pure lsqnonlin cost J = sum||B_model-B_fem||^2
        end
    end
end

function r = resid_all(x, coil, dhat, I, freemask)
    [ell, gB, Khat] = unpack(x, freemask);
    r = [];
    for k = 1:numel(coil)
        pn = coil(k).p / ell;  N = size(pn,1);
        B  = zeros(3*N,1);
        w  = gB * Khat(:, coil(k).pj) * I;
        for i = 1:6
            dx = pn(:,1)-dhat(1,i); dy = pn(:,2)-dhat(2,i); dz = pn(:,3)-dhat(3,i);
            r3 = (dx.^2+dy.^2+dz.^2).^1.5;
            B  = B + w(i) * [dx./r3; dy./r3; dz./r3];
        end
        r = [r; B - coil(k).bfem]; %#ok<AGROW>
    end
end

function [ell, gB, Khat] = unpack(x, freemask)
    ell = x(1)*1e-3;  gB = x(2);
    Khat = zeros(6);  Khat(1,1) = 5/6;  Khat(freemask) = x(3:end);
end
