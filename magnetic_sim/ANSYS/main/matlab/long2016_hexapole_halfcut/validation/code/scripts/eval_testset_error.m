%% eval_testset_error.m — generalization error of each radius-R fit on a COMMON test set
%  For every stored fit {Khat, gB, ell} at R = 150..500 um (from KI_trend_sweep_N1000.mat),
%  evaluate the point-charge model on ONE fixed test set and report
%       error = ( sum_n |B_model - B_FEM| / |B_FEM| * 100 ) / N      (mean per-point % )
%  Test set: real FEM air nodes in the shell 50 um <= ||p|| <= 500 um, 1000 samples
%  POOLED over the 6 coils (point,coil pairs), drawn once with rng(0) -> identical for
%  all fits (fair comparison). No interpolation (per plot-real-nodes rule).
%  Convention identical to sweep_KI_trend_N1000.m (source sign = negate FEM B, k11=5/6).
%
%  Outputs:
%    magnetic_sim/ANSYS/main/MATLAB_data/long2016_hexapole_halfcut/charge_fit/fitting_trend/testset_error_50_500.mat
%    magnetic_sim/ANSYS/main/doc/fitting_trend/testset_error_50_500.tex   (standalone, Overleaf-ready)
clear; clc; close all;

addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
data_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fitting_trend';
doc_dir  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\doc\fitting_trend';

%% constants + pole-tip directions
cnst = mt_constants();
apdl_to_paper_idx = [1, 3, 6, 5, 2, 4];
I_actual = 1;   % [MODIFIED] match FEM excitation (1A); was 0.6, kept in sync with sweep_KI_trend_N1000.m
tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
dhat = tip ./ vecnorm(tip);

%% load all 6 coils once (air-filtered, WP frame, source-sign negated)
fprintf('loading 6 coils ...\n');
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
fprintf('  loaded (coil1 air nodes = %d)\n', size(C(1).P,1));

%% ---- build ONE fixed test set: shell 50..500 um, STRATIFIED 1000 per coil (6000 total) ----
R_in  = 50e-6;   R_out = 500e-6;   Nper = 1000;     % per-coil budget
rng(0);                                              % reproducible draw for the whole set
Ptest = []; Bfem = []; kT = []; pjT = [];
for k = 1:6
    rr2 = sum(C(k).P.^2, 2);
    idx = find(rr2 >= R_in^2 & rr2 <= R_out^2);      % this coil's shell nodes
    sel = randperm(numel(idx), Nper);                % uniform draw of 1000 from this coil
    li  = idx(sel);
    Ptest = [Ptest; C(k).P(li,:)];   %#ok<AGROW>
    Bfem  = [Bfem;  C(k).Bn(li,:)];  %#ok<AGROW>
    kT    = [kT; k*ones(Nper,1)];    %#ok<AGROW>
    pjT   = [pjT; C(k).pj*ones(Nper,1)]; %#ok<AGROW>
    fprintf('  coil%d shell nodes = %d, drew %d\n', k, numel(idx), Nper);
end
Ntest = 6*Nper;                                      % 6000 total
Bfem_mag = vecnorm(Bfem,2,2);
fprintf('test-set = %d (stratified 1000/coil), |B_FEM| range: %.3e .. %.3e T\n', ...
        Ntest, min(Bfem_mag), max(Bfem_mag));

%% ---- load stored fits, select R = 150..500 um ----
S = load(fullfile(data_dir,'KI_trend_sweep_N1000.mat'));   % R_um, Khat_s, gB_s, ell_s, ...
keep = find(S.R_um >= 50 & S.R_um <= 500);
R_um = S.R_um(keep);
nFit = numel(keep);

err_test = zeros(1,nFit);   % the requested mean per-point % error
fprintf('\n  R[um]   ell[mm]   gB          test-err%%\n');
for fi = 1:nFit
    ri   = keep(fi);
    ell  = S.ell_s(ri);   gB = S.gB_s(ri);   Khat = S.Khat_s(:,:,ri);

    % evaluate model on test set, grouped by coil (same pj/weights)
    Bmod = zeros(Ntest,3);
    for k = 1:6
        m = (kT == k);
        if ~any(m); continue; end
        PN = Ptest(m,:) / ell;                 % normalized positions
        w  = gB * Khat(:, apdl_to_paper_idx(k)) * I_actual;   % 6x1 weights
        B  = zeros(sum(m),3);
        for i = 1:6
            D  = PN - dhat(:,i).';
            r3 = (sum(D.^2,2)).^1.5;
            B  = B + w(i) * (D ./ r3);
        end
        Bmod(m,:) = B;
    end

    relerr = vecnorm(Bmod - Bfem, 2, 2) ./ Bfem_mag;   % per-point |dB|/|B_FEM|
    err_test(fi) = mean(relerr) * 100;                 % (sum(...)*100)/N
    fprintf('  %4d   %.3f    %.4e   %.3f\n', R_um(fi), ell*1e3, gB, err_test(fi));
end

%% ---- figure: error vs R ----
fig_dir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\figures\long2016_hexapole_halfcut';
err_mean = mean(err_test);
f = figure('Color','w','Position',[100 100 760 480]);
plot(R_um, err_test, '-o', 'LineWidth',1.6, 'MarkerSize',5, ...
     'MarkerFaceColor',[0.2 0.45 0.85], 'Color',[0.2 0.45 0.85]); hold on;
yline(err_mean, '--', sprintf('mean %.2f%%', err_mean), ...
     'Color',[0.9 0.3 0.2], 'LineWidth',1.5, 'FontSize',11, ...
     'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom');
grid on; box on;
xlabel('Fit sampling radius  R  [\mum]','FontSize',12);
ylabel('error  =  mean_n |B_{model}-B_{FEM}| / |B_{FEM}| \times100  [%]','FontSize',12);
title(sprintf('Generalization error on fixed 50-500 \\mum test set (N=%d nodes)', Ntest),'FontSize',12);
set(gca,'FontSize',11);
exportgraphics(f, fullfile(fig_dir,'testset_error_vs_R.png'), 'Resolution',150);
fprintf('wrote testset_error_vs_R.png to\n  %s\n', fig_dir);

%% save
save(fullfile(data_dir,'testset_error_50_500.mat'), ...
     'R_um','err_test','Ntest','Nper','R_in','R_out','I_actual', ...
     'apdl_to_paper_idx');
fprintf('\nsaved testset_error_50_500.mat\n');

%% ---- write standalone LaTeX (Overleaf-ready) ----
[~,imin] = min(err_test);
fid = fopen(fullfile(doc_dir,'testset_error_50_500.tex'),'w');
fprintf(fid, '%% Auto-generated by eval_testset_error.m\n');
fprintf(fid, '\\documentclass[a4paper,11pt]{article}\n');
fprintf(fid, '\\usepackage{amsmath,amssymb,booktabs,geometry,longtable}\n');
fprintf(fid, '\\geometry{margin=2.2cm}\n');
fprintf(fid, '\\begin{document}\n\n');
fprintf(fid, '\\section*{Charge-model generalization error on a common 50--500\\,$\\mu$m test set}\n');
fprintf(fid, ['\\noindent Model \\texttt{long2016\\_hexapole\\_halfcut}, 6 coils, $I=0.6$\\,A, source sign $=$ FEM $\\mathbf B$ negated. ' ...
    'For each radius $R$ the fit $\\{\\hat{\\mathbf K}_I\\,(\\hat k_{11}{=}5/6),\\hat g_B,\\hat\\ell\\}$ obtained over the ball $\\lVert\\mathbf p\\rVert\\le R$ ' ...
    '(\\texttt{KI\\_trend\\_sweep\\_N1000.mat}) is evaluated on ONE fixed test set: %d real FEM air (point,coil) samples drawn stratified (\\texttt{rng(0)}, %d uniformly per coil) ' ...
    'from the shell $50\\,\\mu\\mathrm{m}\\le\\lVert\\mathbf p\\rVert\\le500\\,\\mu\\mathrm{m}$. The same %d test samples are used for every $R$.\\\\[4pt]\n'], Ntest, Nper, Ntest);
fprintf(fid, ['\\noindent\\textbf{Error metric.}\\quad' ...
    '$\\displaystyle \\mathrm{error}\\;\\triangleq\\;\\frac{1}{N}\\sum_{n=1}^{N}\\frac{\\lVert\\mathbf B_{\\mathrm{model}}(\\mathbf p_n)-\\mathbf B_{\\mathrm{FEM}}(\\mathbf p_n)\\rVert}{\\lVert\\mathbf B_{\\mathrm{FEM}}(\\mathbf p_n)\\rVert}\\times100,\\qquad N=%d.$\\\\[6pt]\n'], Ntest);
fprintf(fid, '\\noindent Point-charge field: $\\mathbf B(\\mathbf p)=\\hat g_B\\sum_{i=1}^{6}\\Big(\\sum_j \\hat k_{ij}I_j\\Big)\\dfrac{\\mathbf p/\\hat\\ell-\\hat{\\mathbf d}_i}{\\lVert\\mathbf p/\\hat\\ell-\\hat{\\mathbf d}_i\\rVert^{3}}.$\n\n');
fprintf(fid, '\\paragraph{Result.} Minimum test error $=%.3f\\%%$ at $R=%d\\,\\mu$m.\n\n', err_test(imin), R_um(imin));
fprintf(fid, '\\begin{longtable}{rr}\n');
fprintf(fid, '\\caption{Mean per-point relative error on the fixed 50--500\\,$\\mu$m test set vs.\\ the fit''s sampling radius $R$.}\\\\\n');
fprintf(fid, '\\toprule\n$R_{\\mathrm{fit}}$ [$\\mu$m] & error [\\%%] \\\\\n\\midrule\n\\endhead\n');
for fi = 1:nFit
    fprintf(fid, '%d & %.3f \\\\\n', R_um(fi), err_test(fi));
end
fprintf(fid, '\\bottomrule\n\\end{longtable}\n\n');
fprintf(fid, '\\end{document}\n');
fclose(fid);
fprintf('wrote testset_error_50_500.tex to\n  %s\n', doc_dir);
