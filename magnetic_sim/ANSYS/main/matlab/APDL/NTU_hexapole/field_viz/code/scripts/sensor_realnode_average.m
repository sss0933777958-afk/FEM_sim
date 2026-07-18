function sensor_realnode_average()
% SENSOR_REALNODE_AVERAGE  NTU sensor 真實節點平均（不內插）→ PDF + .mat
%   直接對 sensor 加密球的實際 FEM 節點取無權重平均（非取樣/內插；per plot-real-nodes）。
%   板下 S1(n+ -z) / 板上 S2(n+ +z)；比較 upper_assembly vs 單極 baseline。
%   ⟨|B|⟩=mean(bsum)、⟨B·n+⟩=mean(sgn*bz)（正=沿外法線 n+、離板），單位 mT（unit-reference）。
%   讀 sim deck 已抽出的球節點 <model>_bfield_sensor_S{1,2}.dat（import_ansys_data 已修 Fortran 負號吃空白）。
%   輸出：results/sensor_realnode_average.pdf + data/sensor_realnode_average.mat。
%   ★ 獨立腳本：不動既有 1000-點內插 sensor_average.m。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';
    fv   = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\NTU_hexapole\field_viz';

    models = {'upper_assembly','singlepole'};
    sens(1) = struct('name','S1','pos',[13.15 -0.99 -0.41],'sgn',-1,'ntxt','-z','where','below plate');
    sens(2) = struct('name','S2','pos',[13.15 -0.99  0.66],'sgn',+1,'ntxt','+z','where','above plate');

    %% ---- 真實節點平均 ----
    R = struct();
    for m = 1:numel(models)
        mdl = models{m};
        for k = 1:2
            d = import_ansys_data(fullfile(root,mdl,'coil1'), ['sensor_' sens(k).name], mdl); % SI [m],[T]
            bmag = d.bsum * 1e3;                 % |B| [mT]
            bnp  = sens(k).sgn * d.bz * 1e3;     % B.n+ [mT]（沿外法線）
            r.N   = numel(d.bz);
            r.Bm  = mean(bmag);  r.Bms = std(bmag);
            r.Bn  = mean(bnp);   r.Bns = std(bnp);
            R.(mdl).(sens(k).name) = r;
            fprintf('%-15s %s (n+=%s): N=%d  <|B|>=%.4f (sd %.4f) mT  <B.n+>=%.4f (sd %.4f) mT\n', ...
                mdl, sens(k).name, sens(k).ntxt, r.N, r.Bm, r.Bms, r.Bn, r.Bns);
        end
    end

    %% ---- 存 .mat 到 field_viz/data/（per matlab-output-layout）----
    datadir = fullfile(fv,'data');
    if ~exist(datadir,'dir'), mkdir(datadir); end
    save(fullfile(datadir,'sensor_realnode_average.mat'), 'R', 'sens', 'models');
    fprintf('saved: %s\n', fullfile(datadir,'sensor_realnode_average.mat'));

    %% ---- 出 PDF（xelatex；比照 sensor_average.m）----
    outdir = fullfile(fv,'results');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    tex = fullfile(outdir,'sensor_realnode_average.tex');
    U = R.upper_assembly;  P = R.singlepole;
    fid = fopen(tex,'w');
    fprintf(fid,'\\documentclass[11pt]{article}\n\\usepackage[margin=1in]{geometry}\n');
    fprintf(fid,'\\usepackage{amsmath,booktabs}\n\\begin{document}\n\\pagestyle{empty}\n');
    fprintf(fid,'\\begin{center}{\\large\\bfseries NTU hexapole upper assembly --- Hall-sensor real-node average}\\\\[2pt]\n');
    fprintf(fid,'1\\,A $\\times$ 50 turns, coil on post0 \\quad$|$\\quad average over the actual FEM nodes of each sensor sphere ($R=150\\,\\mu$m)\\end{center}\n');
    fprintf(fid,'\\vspace{6pt}\n\\begin{center}\n\\begin{tabular}{llcccc}\n\\toprule\n');
    fprintf(fid,'Sensor & center $(x,y,z)$ [mm] & $\\langle|B|\\rangle$ upper [mT] & $\\langle B\\cdot n_+\\rangle$ upper [mT] & $\\langle|B|\\rangle$ single [mT] & ratio\\\\\n\\midrule\n');
    for k = 1:2
        u = U.(sens(k).name);  p = P.(sens(k).name);
        fprintf(fid,'%s (%s) & $(%.2f,\\,%.2f,\\,%+.2f)$ & $%.4f\\;(\\pm%.4f)$ & $%+.4f\\;(\\pm%.4f)$ & $%.4f$ & $%.2f\\times$\\\\\n', ...
            sens(k).name, sens(k).where, sens(k).pos(1),sens(k).pos(2),sens(k).pos(3), ...
            u.Bm,u.Bms, u.Bn,u.Bns, p.Bm, u.Bm/p.Bm);
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{center}\n');
    fprintf(fid,'\\vspace{4pt}\\noindent\\small Real-node average (no interpolation): unweighted mean over the actual FEM nodes inside each sensor sphere ($N=%d$ for S1, $%d$ for S2); $(\\pm)$ = std over those nodes. $B\\cdot n_+$: sign along outward normal $n_+$ (away from the pole plate; S1 $n_+=-z$ below, S2 $n_+=+z$ above). ratio = $\\langle|B|\\rangle_{\\text{upper}}/\\langle|B|\\rangle_{\\text{single}}$ (single = single-pole baseline, same real-node method).\n', U.S1.N, U.S2.N);
    fprintf(fid,'\\end{document}\n');
    fclose(fid);

    xelatex = 'C:\Users\Kuo\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe';
    old = cd(outdir);
    system(sprintf('"%s" -interaction=nonstopmode -halt-on-error "%s"', xelatex, tex));
    cd(old);
    pdf = fullfile(outdir,'sensor_realnode_average.pdf');
    if exist(pdf,'file')
        for e = {'aux','log'}
            f = fullfile(outdir,['sensor_realnode_average.' e{1}]);
            if exist(f,'file'), delete(f); end
        end
        delete(tex);
        fprintf('saved: %s\n', pdf);
    else
        error('xelatex 編譯失敗（.tex 留在 results/ 供查）');
    end
end
