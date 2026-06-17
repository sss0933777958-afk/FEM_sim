%% fit_KI_full.m — hexapole point-charge model fit (per the document)
%  Model    : B(p) = gB * sum_i ( sum_j khat_ij I_j ) * (p/ell - dhat_i)/||p/ell - dhat_i||^3
%  Charges  : pc_i = ell * dhat_i  (dhat_i = unit pole-tip direction, paper P1..P6)
%  Fit vars : Khat_I^FEM (6x6, with Khat(1,1) fixed = 5/6), ell, gB
%  Init     : Khat0 = eye(6)-ones(6)/6 ,  ell0 = 0.5 mm ,  gB0 = 1
%  Cost     : J = sum_k || B_FEM(p_k) - B_model(p_k) ||^2  over all 6 coils
%  Data     : long2016_hexapole_halfcut, 6 coils baseline, 0.6 A.

clear; clc; close all;

%% 0. helpers + paths
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
out_dir      = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit';
if ~exist(out_dir,'dir'); mkdir(out_dir); end

%% 1. constants + pole-tip directions dhat_i
cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];      % coil k excites this paper pole
cube_half = 100e-6;                           % +/-100 um fitting cube
I_actual  = 0.6;                              % drive current [A]
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];   % 3x6 [m]
dhat = tip ./ vecnorm(tip);                   % 3x6 unit directions

%% 2. load 6-coil FEM B-field (WP frame, air nodes, cube)
coil = struct('p',{},'bfem',{},'bmag',{},'pj',{});
for k = 1:6
    cname = sprintf('coil%d', k);
    d   = import_ansys_data(fullfile(results_root, cname), 'wp', cname);
    air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize', false));
    zwp = d.z - cnst.SPH_OFST;
    m   = air & abs(d.x)<cube_half & abs(d.y)<cube_half & abs(zwp)<cube_half;
    coil(k).p    = [d.x(m), d.y(m), zwp(m)];
    % source convention: excited pole tip RADIATES outward (positive charge).
    % FEM winding gives the opposite sign, so negate B to match the document's
    % Q = +(N_c/mu0 R_a) K_I I convention -> gB and q come out positive.
    coil(k).bfem = -[d.bx(m); d.by(m); d.bz(m)];
    coil(k).bmag = sqrt(d.bx(m).^2 + d.by(m).^2 + d.bz(m).^2);
    coil(k).pj   = apdl_to_paper_idx(k);
end

%% 3. initial guess (per document)
Khat0 = eye(6) - ones(6)/6;     % diag 5/6, off-diag -1/6
ell0  = 0.5e-3;                 % 0.5 mm
gB0   = 1;
freemask = true(6); freemask(1,1) = false;       % Khat(1,1) fixed = 5/6

%% 4. minimize J over {Khat, ell, gB}  (lsqnonlin)
% pack x = [ell_mm; gB; Khat(free 35)]
x0   = [ell0*1e3; gB0; Khat0(freemask)];
opts = optimoptions('lsqnonlin','Display','iter', ...
    'MaxFunctionEvaluations',1e5,'MaxIterations',4e3, ...
    'FunctionTolerance',1e-20,'StepTolerance',1e-12);
xfit = lsqnonlin(@(x) resid_all(x, coil, dhat, I_actual, freemask), x0, [], [], opts);
[ell, gB, Khat] = unpack(xfit, freemask);
J = sum(resid_all(xfit, coil, dhat, I_actual, freemask).^2);

%% 5. report
fprintf('\nell = %.3f mm   gB = %.4e   J = %.4e\n\n', ell*1e3, gB, J);
fprintf('Khat_I^FEM =\n');
for i = 1:6, fprintf('  %9.4f', Khat(i,:)); fprintf('\n'); end

%% 5b. adopt k11 = 5/6 -> back out R_a and charges q_i
%   gB = N_c*(6/5)*k11 / (4*pi*R_a*ell^2)  ->  R_a = N_c*(6/5)*k11/(4*pi*|gB|*ell^2)
%   q_i = gB*ell^2*(sum_j khat_ij I_j)/k_m   (per coil)
k11 = 5/6;
k_m = cnst.k_m;                                  % mu_0/4pi = 1e-7
R_a = cnst.N_c*(6/5)*k11 / (4*pi*abs(gB)*ell^2);
fprintf('\nk11 = %.4f (assumed)    R_a = %.4e A/Wb\n', k11, R_a);

Q = zeros(6,6);                                  % q_i per coil [A*m]
for k = 1:6, Q(:,k) = gB*ell^2*(Khat(:,coil(k).pj)*I_actual)/k_m; end
fprintf('charges q_i [A*m] (rows P1..P6, cols coil1..6):\n');
for i = 1:6, fprintf('  %10.3e', Q(i,:)); fprintf('\n'); end

%% 6. save
save(fullfile(out_dir,'fit_KI_full.mat'), 'Khat','ell','gB','J', 'k11','R_a','Q', ...
     'Khat0','ell0','gB0','coil','cnst','cube_half','I_actual','apdl_to_paper_idx');


%% ===== local functions =====
function r = resid_all(x, coil, dhat, I, freemask)
    [ell, gB, Khat] = unpack(x, freemask);
    r = [];
    for k = 1:numel(coil)
        pn = coil(k).p / ell;
        N  = size(pn,1);
        B  = zeros(3*N,1);
        w  = gB * Khat(:, coil(k).pj) * I;         % 6x1 charge weights
        for i = 1:6
            dx = pn(:,1)-dhat(1,i); dy = pn(:,2)-dhat(2,i); dz = pn(:,3)-dhat(3,i);
            r3 = (dx.^2+dy.^2+dz.^2).^1.5;
            B  = B + w(i) * [dx./r3; dy./r3; dz./r3];
        end
        r = [r; B - coil(k).bfem];                 %#ok<AGROW>
    end
end

function [ell, gB, Khat] = unpack(x, freemask)
    ell = x(1)*1e-3;  gB = x(2);
    Khat = zeros(6);  Khat(1,1) = 5/6;  Khat(freemask) = x(3:end);
end
