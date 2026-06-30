%% make_R040_charge_figs.m
%  Produce the 4 equivalent-charge figures at R=40um (new graded all-source gB50 fit):
%    P1_circuit_charge_R40_zoom.png   (2D circuit + charge)
%    P2_circuit_charge_R40_zoom.png   (2D circuit + charge)
%    P1_charge_only_R40_3D.png        (3D charge position)
%    P2_charge_only_R40_3D.png        (3D charge position)
%  The plot functions only read 'ell' from fit_KI_R040.mat, so we materialise that file
%  from the R=40um entry of KI_convergence_gB50.mat (ell_hat = 0.844 mm). Field arrows in
%  the 2D plots come from the old coil1/coil5 _all data (the new graded sims only dumped the
%  2mm WP sphere); the R-dependent quantity here is the charge position ell*dhat.
clear; clc; close all;
plotdir = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\plot';
addpath(plotdir);

S  = load('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\Calibration_using_FEM_modeling\fix_dir\data\KI_convergence_gB50.mat');
ri = find(round(S.R_list*1e6)==40,1);
f  = S.fitP(ri);
ell = f.ell; gB = f.gB; Khat = f.Khat; I_actual = S.I_actual;     %#ok<NASGU>
R = 40e-6; Nmax = S.Neq; apdl_to_paper_idx = [1 3 6 5 2 4];         %#ok<NASGU>
matf = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\Calibration_using_FEM_modeling\fix_dir\data\fit_KI_R040.mat';
save(matf, 'ell','gB','Khat','I_actual','R','Nmax','apdl_to_paper_idx');
fprintf('wrote %s  (ell = %.3f mm)\n', matf, ell*1e3);

plot_P1_circuit_charge('zoom', true, 40);   % 2D P1 circuit + charge
plot_P2_circuit_charge(true, 40);           % 2D P2 circuit + charge
plot_pole_circuit_charge_3d(1, false, 40);  % 3D P1 charge position
plot_pole_circuit_charge_3d(2, false, 40);  % 3D P2 charge position
