%% sweep_KI_radius.m — K_I charge-model fit over BALL sampling regions
%  Sampling region = ball ||p|| <= R about the working point (WP), in WP frame.
%  R = 50:50:500 um (10 regions).  For each R: fit {Khat_I(6x6, (1,1)=5/6), ell, gB}
%  by lsqnonlin over all 6 coils, then emit one LaTeX script.
%  Model/convention identical to fit_KI_full.m (source sign = negate FEM B).
%
%  Outputs:
%   - .mat : magnetic_sim/ANSYS/main/MATLAB_data/long2016_hexapole_halfcut/charge_fit/fit_KI_ball/fit_KI_R<R>.mat  (per R)
%            + fit_KI_ball_sweep.mat (summary)
%   - .tex : magnetic_sim/ANSYS/main/doc/charge_model_fitting/long2016_hexapole_halfcut/scripts/
%            fit_KI_R<R>.tex (per R)  + fit_KI_ball_summary.tex (table)

clear; clc; close all;

%% paths
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut';
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fit_KI_ball';
tex_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\doc\charge_model_fitting\long2016_hexapole_halfcut\scripts';
if ~exist(data_dir,'dir'); mkdir(data_dir); end

%% constants + pole-tip directions
cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];
I_actual = 0.6;
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

%% load all 6 coils ONCE (full air-filtered WP-frame nodes)
fprintf('loading 6 coils ...\n');
C = struct('P',{},'Bn',{},'pj',{});
for k = 1:6
    cname = sprintf('coil%d', k);
    d   = import_ansys_data(fullfile(results_root, cname), 'wp', cname);
    air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize', false));
    zwp = d.z - cnst.SPH_OFST;
    C(k).P  = [d.x(air), d.y(air), zwp(air)];      % N x 3  [m], WP frame
    C(k).Bn = -[d.bx(air), d.by(air), d.bz(air)];  % N x 3  source-sign (negated)
    C(k).pj = apdl_to_paper_idx(k);
end
fprintf('  loaded (coil1 air nodes = %d)\n', size(C(1).P,1));

%% fit settings
R_list = (50:50:500) * 1e-6;        % ball radii [m]
Nmax   = 4000;                      % subsample cap per coil per region
Khat0  = eye(6) - ones(6)/6;
ell0   = 0.5e-3;  gB0 = 1;
freemask = true(6); freemask(1,1) = false;
opts = optimoptions('lsqnonlin','Display','off', ...
    'MaxFunctionEvaluations',1e5,'MaxIterations',4e3, ...
    'FunctionTolerance',1e-20,'StepTolerance',1e-12);

k_m = cnst.k_m;  k11 = 5/6;
summary = struct('R_um',{},'Ntot',{},'ell_mm',{},'gB',{},'Ra',{},'errpct',{});

%% sweep
for ri = 1:numel(R_list)
    R = R_list(ri); R_um = round(R*1e6);
    % build per-coil sampled sets inside the ball (subsampled, reproducible)
    coil = struct('p',{},'bfem',{},'pj',{});
    Ntot = 0; sumB2 = 0;
    for k = 1:6
        r2 = sum(C(k).P.^2, 2);
        idx = find(r2 < R^2);
        rng(0);                                  % reproducible subsample
        if numel(idx) > Nmax, idx = idx(randperm(numel(idx), Nmax)); end
        p  = C(k).P(idx,:);
        Bn = C(k).Bn(idx,:);
        coil(k).p    = p;
        coil(k).bfem = [Bn(:,1); Bn(:,2); Bn(:,3)];   % 3N x 1
        coil(k).pj   = C(k).pj;
        Ntot  = Ntot + numel(idx);
        sumB2 = sumB2 + sum(coil(k).bfem.^2);
    end

    % fit
    x0 = [ell0*1e3; gB0; Khat0(freemask)];
    xfit = lsqnonlin(@(x) resid_all(x, coil, dhat, I_actual, freemask), x0, [], [], opts);
    [ell, gB, Khat] = unpack(xfit, freemask);
    J = sum(resid_all(xfit, coil, dhat, I_actual, freemask).^2);
    errpct = 100*sqrt(J / sumB2);                % relative RMS error over region
    R_a = cnst.N_c*(6/5)*k11 / (4*pi*abs(gB)*ell^2);
    Q = zeros(6,6);
    for k = 1:6, Q(:,k) = gB*ell^2*(Khat(:,coil(k).pj)*I_actual)/k_m; end

    fprintf('R=%3d um | N=%6d | ell=%.3f mm | gB=%.4e | err=%.2f%%\n', ...
            R_um, Ntot, ell*1e3, gB, errpct);

    % save per-R mat
    save(fullfile(data_dir, sprintf('fit_KI_R%03d.mat', R_um)), ...
         'R','R_um','Khat','ell','gB','J','errpct','k11','R_a','Q','Ntot', ...
         'Khat0','ell0','gB0','I_actual','apdl_to_paper_idx','Nmax');

    % emit per-R LaTeX (three fitted quantities, ball region, error)
    write_tex(fullfile(tex_dir, sprintf('fit_KI_R%03d.tex', R_um)), ...
              R_um, Khat, ell, gB, errpct, Ntot);

    summary(ri) = struct('R_um',R_um,'Ntot',Ntot,'ell_mm',ell*1e3, ...
                         'gB',gB,'Ra',R_a,'errpct',errpct);
end

%% combined summary
save(fullfile(data_dir,'fit_KI_ball_sweep.mat'), 'summary','R_list','Nmax','I_actual');
write_summary_tex(fullfile(tex_dir,'fit_KI_ball_summary.tex'), summary);
fprintf('\ndone: %d per-R .tex + summary in\n  %s\n', numel(R_list), tex_dir);


%% ===== local functions =====
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

function write_tex(fname, R_um, Khat, ell, gB, errpct, Ntot)
    % presentation: flip upper-pole coil columns so all diagonals positive
    Khat(:, [2 4 5]) = -Khat(:, [2 4 5]);
    fid = fopen(fname,'w');
    fprintf(fid,'%% Auto-generated by sweep_KI_radius.m -- K_I fit, ball R=%d um\n', R_um);
    fprintf(fid,'\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n\n');
    fprintf(fid,['Fitting (sampling) region: ball $\\lVert\\mathbf{p}\\rVert \\le %d\\,' ...
                 '\\mu\\mathrm{m}$ about the working point (WP); %d sample nodes ' ...
                 '(6 coils).\n\n'], R_um, Ntot);
    fprintf(fid,['Sign convention: every excited pole shown as a source (upper-pole ' ...
                 'P2, P4, P5 coil columns flipped), so all diagonal entries are ' ...
                 'positive (emit) and all off-diagonal entries negative (receive).\n\n']);
    fprintf(fid,'\\[\n\\widehat{\\mathbf{K}}_I^{\\mathrm{FEM}} =\n\\begin{bmatrix}\n');
    for i = 1:6
        fprintf(fid,'%.4f', Khat(i,1));
        for j = 2:6, fprintf(fid,' & %.4f', Khat(i,j)); end
        if i < 6, fprintf(fid,' \\\\\n'); else, fprintf(fid,'\n'); end
    end
    fprintf(fid,'\\end{bmatrix}\n\\]\n\n');
    fprintf(fid,'\\[\n\\widehat{\\ell} = %.3f~\\mathrm{mm}, \\qquad \\widehat{g}_B = %s.\n\\]\n\n', ...
            ell*1e3, sci(gB));
    fprintf(fid,'Relative RMS field error over the region: $%.2f\\%%$.\n\n', errpct);
    fprintf(fid,'\\end{document}\n');
    fclose(fid);
end

function write_summary_tex(fname, S)
    fid = fopen(fname,'w');
    fprintf(fid,'%% Auto-generated by sweep_KI_radius.m -- ball-radius sweep summary\n');
    fprintf(fid,'\\documentclass{article}\n\\usepackage{amsmath,booktabs}\n\\begin{document}\n\n');
    fprintf(fid,'K$_I$ charge-model fit vs spherical sampling radius (about WP).\n\n');
    fprintf(fid,'\\begin{tabular}{rrrrr}\n\\toprule\n');
    fprintf(fid,'$R~[\\mu\\mathrm{m}]$ & $N$ & $\\widehat{\\ell}~[\\mathrm{mm}]$ & $\\widehat{g}_B$ & err [\\%%] \\\\\n\\midrule\n');
    for i = 1:numel(S)
        fprintf(fid,'%d & %d & %.3f & $%s$ & %.2f \\\\\n', ...
                S(i).R_um, S(i).Ntot, S(i).ell_mm, sci(S(i).gB), S(i).errpct);
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\n\\end{document}\n');
    fclose(fid);
end

function s = sci(x)
    if x == 0, s = '0'; return; end
    e = floor(log10(abs(x)));
    s = sprintf('%.3f\\times10^{%d}', x/10^e, e);
end
