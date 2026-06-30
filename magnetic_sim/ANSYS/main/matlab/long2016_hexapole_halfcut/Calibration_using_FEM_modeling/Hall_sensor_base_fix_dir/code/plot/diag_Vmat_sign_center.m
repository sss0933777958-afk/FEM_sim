function diag_Vmat_sign_center(variant)
% DIAG_VMAT_SIGN_CENTER  用「圓柱底面中心單點內插」法算 6×6 all-source Vmat，印 sign 表。
%   variant (選填)：FEM 變體子夾，預設 'standard'（baseline）；可傳 'gap200um_mueq' 對照。
%   對照舊法（粗網格圓柱平均 / 1 節點），看 P2-under-P1 等 off-diag sign 是否改變。
    if nargin < 1 || isempty(variant), variant = 'standard'; end
    TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
            'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir'];
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(fullfile(TREE,'code','function'));
    rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
    mesh_csv_dir = fullfile(rr,'mesh','standard','csv');   % 連接性恆用 baseline 網格（gap200um 同網格，csv/ 子夾）

    cnst = mt_constants();
    apdl_to_paper_idx = [1,3,6,5,2,4];
    S_hall = 130;
    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);
    [Vmat, ~, Bn] = extract_Vmat_interp_center(rr, cnst, apdl_to_paper_idx, ...
                                               sensor_pos, sensor_n, S_hall, mesh_csv_dir, variant);

    colpole = {'P1','P3','P6','P5','P2','P4'};   % APDL coil1..6 → paper 激發極
    paper_of_col = [1 3 6 5 2 4];
    fprintf('\n=== [%s] 單點內插 all-source Vmat sign（列=sensor P1..P6，欄=激發極）===\n', variant);
    fprintf('          ');  for c=1:6, fprintf(' exc%-3s', colpole{c}); end; fprintf('\n');
    for i=1:6
        fprintf('  senP%d :', i);
        for j=1:6, fprintf('   %+d ', sign(Vmat(i,j))); end
        fprintf('\n');
    end
    fprintf('\nP2 sensor under P1 激發 = Vmat(2,1) = %+.3e V  (all-source 單點 B·n+ = %+.3e T)  → %s\n', ...
            Vmat(2,1), Bn(2,1), tern(Vmat(2,1)<0,'負(匯進P2)','正(回流)'));
    fprintf('\n3 對 opposing pair（upper←lower）單點 B·n+：\n');
    pairs = [2 1; 4 2; 5 3];   % (sensor row, exc col)：P2←P1(col1=P1)、P4←P3(col2=P3)、P5←P6(col3=P6)
    plab  = {'P2←P1','P4←P3','P5←P6'};
    for k=1:3
        fprintf('   %s : Vmat=%+.3e V  ｜ all-source B·n+=%+.3e T\n', plab{k}, ...
                Vmat(pairs(k,1),pairs(k,2)), Bn(pairs(k,1),pairs(k,2)));
    end
    fprintf('\noff-diagonal 為正者：\n');
    for i=1:6, for j=1:6
        if paper_of_col(j)~=i && Vmat(i,j)>0
            fprintf('   senP%d ← exc%s : %+.3e\n', i, colpole{j}, Vmat(i,j));
        end
    end, end
end
function o=tern(c,a,b), if c,o=a; else,o=b; end, end
