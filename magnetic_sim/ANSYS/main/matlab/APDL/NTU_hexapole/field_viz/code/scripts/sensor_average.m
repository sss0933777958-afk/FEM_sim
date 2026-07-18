function sensor_average(model, cen_mm, tag)
% SENSOR_AVERAGE  NTU：兩 sensor 有效範圍內 1000 點面積平均（重心內插）→ PDF
%   取樣圓柱（比照 extract_Vmat_interp）：底面=sensor 中心、沿 n+ 高 H=0.1mm、半徑 R=0.15mm，
%   面積均勻取 N=1000 點；對局部真實 FEM 節點建 Delaunay → 重心內插 B(p)=Σλ_iB_i（=FEM 線性 shape fn）。
%   輸出 ⟨B·n+⟩、⟨|B|⟩（mT，unit-reference）→ results/sensor_average_<model><tag>.pdf。
%   model  ：'singlepole'（預設）或 'upper_assembly'。
%   cen_mm ：可選 2x3 sensor 中心 [mm]（覆蓋預設 (13.15,-0.99,∓)）；n+ 慣例維持 S1=-z / S2=+z。
%   tag    ：可選輸出檔名後綴（區隔不同 sensor 位置，例 '_x8p15'）。
    if nargin<1 || isempty(model), model='singlepole'; end
    if nargin<2 || isempty(cen_mm), cen_mm=[13.15 -0.99 -0.41; 13.15 -0.99 0.66]; end
    if nargin<3, tag=''; end
    rng(1);
    N = 1000;  R = 0.15e-3;  H = 0.10e-3;          % 取樣點數 / 圓柱半徑 / 高 [m]
    OVER = 2000;                                   % 過取樣（凸包外 NaN 剔除後仍夠 1000）

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    rr = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';
    d  = import_ansys_data(fullfile(rr,model,'coil1'),'all',model);   % SI [m],[T]
    X = [d.x d.y d.z];  Bn3 = [d.bx d.by d.bz];

    S(1) = struct('name','S1','c',cen_mm(1,:)*1e-3,'n',[0 0 -1],'ntxt','-z', ...
                  'nf',0,'Bn',0,'Bns',0,'Bm',0,'Bms',0);
    S(2) = struct('name','S2','c',cen_mm(2,:)*1e-3,'n',[0 0 +1],'ntxt','+z', ...
                  'nf',0,'Bn',0,'Bns',0,'Bm',0,'Bms',0);

    for k=1:2
        c=S(k).c(:); ni=S(k).n(:);
        loc = vecnorm(X-c.',2,2) < 0.30e-3;                 % 局部節點（球內+一點外）建 Delaunay
        P=X(loc,:); BP=Bn3(loc,:);
        TR=delaunayTriangulation(P);
        Tb=null(ni.'); t1=Tb(:,1); t2=Tb(:,2);              % 圓柱局部基底 ⊥ n+
        a=H*rand(OVER,1); r=R*sqrt(rand(OVER,1)); th=2*pi*rand(OVER,1);  % 面積均勻
        pts = c.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';
        ti=pointLocation(TR,pts); good=find(~isnan(ti));
        good=good(1:min(N,numel(good)));                    % 取前 N=1000 個有效點
        bc=cartesianToBarycentric(TR,ti(good),pts(good,:));
        conn=TR.ConnectivityList(ti(good),:);
        Bint=zeros(numel(good),3);
        for cc=1:4, Bint=Bint+bc(:,cc).*BP(conn(:,cc),:); end
        Bmag=vecnorm(Bint,2,2)*1e3;  Bproj=(Bint*ni)*1e3;   % mT
        S(k).nf=numel(good);
        S(k).Bn=mean(Bproj); S(k).Bns=std(Bproj);
        S(k).Bm=mean(Bmag);  S(k).Bms=std(Bmag);
        fprintf('%s (n+=%s): %d pts  <B.n+>=%.4f (sd %.4f) mT  <|B|>=%.4f (sd %.4f) mT\n', ...
                S(k).name,S(k).ntxt,S(k).nf,S(k).Bn,S(k).Bns,S(k).Bm,S(k).Bms);
    end

    %% ---- 出 PDF（xelatex）----
    outdir='G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\NTU_hexapole\field_viz\results';
    base=['sensor_average_' model tag];
    mtxt=strrep(model,'_','\_');
    tex=fullfile(outdir,[base '.tex']);
    fid=fopen(tex,'w');
    fprintf(fid,'\\documentclass[11pt]{article}\n\\usepackage[margin=1in]{geometry}\n');
    fprintf(fid,'\\usepackage{amsmath,booktabs}\n\\begin{document}\n\\pagestyle{empty}\n');
    fprintf(fid,'\\begin{center}{\\large\\bfseries NTU %s --- Hall-sensor field average (interpolated)}\\\\[2pt]\n',mtxt);
    fprintf(fid,'1\\,A $\\times$ 50 turns \\quad$|$\\quad %d-point area average over sampling cylinder ($R=150\\,\\mu$m, $H=100\\,\\mu$m along $n_+$)\\end{center}\n',N);
    fprintf(fid,'\\vspace{6pt}\n\\begin{center}\n\\begin{tabular}{lccccc}\n\\toprule\n');
    fprintf(fid,'Sensor & center $(x,y,z)$ [mm] & $n_+$ & $\\langle B\\cdot n_+\\rangle$ [mT] & $\\langle|B|\\rangle$ [mT] & pts\\\\\n\\midrule\n');
    for k=1:2
        fprintf(fid,'%s & $(%.2f,\\,%.2f,\\,%+.2f)$ & $%s$ & $%+.4f\\;(\\pm%.4f)$ & $%.4f\\;(\\pm%.4f)$ & %d\\\\\n', ...
            S(k).name, S(k).c(1)*1e3,S(k).c(2)*1e3,S(k).c(3)*1e3, S(k).ntxt, ...
            S(k).Bn,S(k).Bns, S(k).Bm,S(k).Bms, S(k).nf);
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{center}\n');
    fprintf(fid,'\\vspace{4pt}\\noindent\\small $(\\pm)$ = std over the %d sampled points. $B\\cdot n_+$: sign along outward normal $n_+$ (away from pole plate). Barycentric interpolation on the local Delaunay of real FEM nodes.\n',N);
    fprintf(fid,'\\end{document}\n');
    fclose(fid);

    xelatex='C:\Users\Kuo\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe';
    old=cd(outdir);
    system(sprintf('"%s" -interaction=nonstopmode -halt-on-error "%s"',xelatex,tex));
    cd(old);
    pdf=fullfile(outdir,[base '.pdf']);
    if exist(pdf,'file')
        delete(fullfile(outdir,[base '.aux']));
        delete(fullfile(outdir,[base '.log']));
        delete(tex);
        fprintf('saved: %s\n',pdf);
    else
        error('xelatex 編譯失敗（.tex 留在 results/ 供查）');
    end
end
