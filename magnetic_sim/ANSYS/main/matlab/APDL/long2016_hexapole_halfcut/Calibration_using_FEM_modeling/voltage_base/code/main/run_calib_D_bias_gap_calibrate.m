function run_calib_D_bias_gap_calibrate()
% run_calib_D_bias_gap_calibrate — 算 gap_calibrate 的 **18-param(bias) voltage 校正**。
%   等同 voltage_base/main.m 但 USE_BIAS=true，存到 **calib_D_gap_calibrate_bias.mat**（不覆蓋舊 single 檔）。
%   目的：補齊缺的 18-param D̄/Vmat/Pc 供電壓力模型用（見 plot_circle_force_voltage_18p.m）。
%   重用 current_base/voltage_base 既有函式（load_coils_actuator/fit_varpro/build_A/make_Pc/
%   build_sensor_geometry/extract_Vmat_interp）。

    USE_BIAS = true;  VARIANT = 'gap_calibrate';  N_I = 6;
    R_select = 150e-6; I_actual = 1; S_hall = 130; ell0 = 0.5e-3;
    n_uniform = 10000; AXIAL_TOL = 0.10e-3; dataset = 'all';

    CAL = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
           'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
    addpath(fullfile(CAL,'current_base','code','main_function'));   % load_coils_actuator/select_ball/build_S_matrix/fitting
    addpath(fullfile(CAL,'voltage_base','code','main_function'));   % build_sensor_geometry/build_V_matrix/solve_D_bar_gain
    model = 'long2016_hexapole_halfcut';
    results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
    mesh_csv_dir = fullfile(fileparts(results_root),'csv',VARIANT);
    mat_dir = fullfile(CAL,'voltage_base','data');

    cnst = mt_constants();  apdl_to_paper_idx = [1,3,6,5,2,4];

    D = load_coils_actuator(model, cnst, apdl_to_paper_idx, dataset, VARIANT);
    [P, Bstack, npts] = select_ball(D, R_select);
    s_sink = ones(1,6); for j=1:6, if ismember(apdl_to_paper_idx(j),[1 3 6]), s_sink(j)=-1; end; end
    Bstack = (-Bstack).*s_sink;                              % all-source

    [ell_hat, e_hat, Pc, J] = fitting(P, Bstack, D.Pc_base, ell0, USE_BIAS);
    fprintf('BIAS fit: ell=%.2f µm, ||e_hat||=%.4e, J=%.4e\n', ell_hat*1e6, norm(e_hat), J);

    A  = build_S_matrix(ell_hat, Pc, P);  M = A.'*A;  Dv = M\(A.'*Bstack);
    ell_hat = ell_hat*1e6;

    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);
    [Vmat, exc_sign] = build_V_matrix(results_root, cnst, apdl_to_paper_idx, ...
                           sensor_pos, sensor_n, S_hall, mesh_csv_dir, n_uniform, [], AXIAL_TOL, VARIANT);

    [D_bar, ghat_V_B, H_V] = solve_D_bar_gain(Dv, Vmat);

    [~, paper_to_apdl] = sort(apdl_to_paper_idx);
    Vmat_p = Vmat(:, paper_to_apdl);  Dv_p = Dv(:, paper_to_apdl);
    E36 = zeros(3,6);
    E36(:,1)=e_hat(1:3); E36(:,2)=e_hat(4:6); E36(:,3)=e_hat(7:9); E36(:,4)=e_hat(10:12); E36(:,5)=e_hat(13:15);
    E36(1,6)=e_hat(16); E36(2,6)=e_hat(17); E36(3,6)=e_hat(1)-e_hat(4)+e_hat(8)-e_hat(11)+e_hat(15);

    Dmat = H_V;  g_V = ghat_V_B;  Pc_base = D.Pc_base;  R_act = D.R_act;
    out = fullfile(mat_dir, 'calib_D_gap_calibrate_bias.mat');
    save(out, 'Dmat','D_bar','g_V','Dv_p','Vmat_p','exc_sign','ell_hat','e_hat','E36','Pc','Pc_base','R_act', ...
              'J','S_hall','R_select','npts','VARIANT','USE_BIAS','apdl_to_paper_idx','paper_to_apdl', ...
              'sensor_pos','sensor_n','n_uniform');
    fprintf('已存 %s  (D_bar(1,1)=%.4f, ell=%.1f µm, g_V=%.4e mT/mV)\n', out, D_bar(1,1), ell_hat, g_V);
end
