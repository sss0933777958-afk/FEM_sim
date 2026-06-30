function diag_Vmat_sign()
% DIAG_VMAT_SIGN  印 baseline 6×6 all-source Vmat 的 sign 表，找出哪些 off-diagonal 為正。
%   目的：確認「P2-under-P1 讀正」是不是唯一異常，還是 P4/P5（其他上極）也有。
%   列=sensor 極 P1..P6；欄=激發（APDL coil1..6 = paper P1,P3,P6,P5,P2,P4）。
    TREE = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
            'long2016_hexapole_halfcut\Calibration_using_FEM_modeling\Hall_sensor_base_fix_dir'];
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    addpath(fullfile(TREE,'code','function'));
    rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';

    cnst = mt_constants();
    apdl_to_paper_idx = [1,3,6,5,2,4];
    S_hall = 130;
    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);
    Vmat = extract_Vmat(rr, cnst, apdl_to_paper_idx, sensor_pos, sensor_n, S_hall, 'standard');

    colpole = {'P1','P3','P6','P5','P2','P4'};   % APDL coil1..6 → paper 激發極
    fprintf('\n=== baseline all-source Vmat sign（列=sensor P1..P6，欄=激發極）===\n');
    fprintf('          ');  for c=1:6, fprintf(' exc%-3s', colpole{c}); end; fprintf('\n');
    for i=1:6
        fprintf('  senP%d :', i);
        for j=1:6
            s = sign(Vmat(i,j));
            if i==find(strcmp(colpole,sprintf('P%d',i))), tag='[D]'; else tag='   '; end %#ok<NASGU>
            fprintf('   %+d ', s);
        end
        fprintf('\n');
    end
    % 重點：P2 sensor(列2) under P1 激發(欄1)
    fprintf('\nP2 sensor under P1 激發 = Vmat(2,1) = %+.3e  → %s\n', ...
            Vmat(2,1), tern(Vmat(2,1)<0,'負(匯進P2)','正(回流)'));
    % 列出所有「off-diagonal 為正」的 (sensor, 激發) 對
    paper_of_col = [1 3 6 5 2 4];
    fprintf('\n所有 off-diagonal 為正的項（理論上 all-source off-diag 應幾乎全負）：\n');
    for i=1:6, for j=1:6
        if paper_of_col(j)~=i && Vmat(i,j)>0
            fprintf('   senP%d ← exc%s : %+.3e\n', i, colpole{j}, Vmat(i,j));
        end
    end, end
end
function o=tern(c,a,b), if c,o=a; else,o=b; end, end
