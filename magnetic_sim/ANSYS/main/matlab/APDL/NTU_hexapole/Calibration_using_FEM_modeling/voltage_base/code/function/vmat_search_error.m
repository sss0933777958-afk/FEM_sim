function vmat_search_error(variant, target_pct)
%VMAT_SEARCH_ERROR  NTU：sensing search-area 內找「V-matrix 相差 ~target% 且自激發電壓為正」的點 + Frobenius。
% -------------------------------------------------------------------------
% 問題：Hall sensor 可在其 search area（外圓 r_max=2.865mm）內任意合法位置擺放。每顆 sensor 各在自己的
%   search 圓盤內找一點，使其那一列相對偏差 ≈ target%（預設 50%）**且**自激發電壓（該列自己欄）為正
%   （激發極自己的電壓為正的號誌慣例）。用這 6 點組出 V，算相對 Frobenius 誤差（≈target%）。
%   V_pdf = data/calib_D_<variant>.mat 的 Vmat_p（已套自激發全正號誌）。target_pct=[] → 改找「最大偏差」。
%
% 號誌：自激發全正 = 對 raw sensor 讀值套 col_sign（翻下極激發欄 P1/P3/P6），與 main.m 一致
%   （exc_sign·col_sign = +1 → 等同 raw；G 與 V 同翻，D̄ 不變）。
%
% 兩階段：①搜尋（快，點值）每 sensor 建 1 次 Delaunay+pointLocation（6 激發共用權重）、候選=圓柱底面中心
%   （n+=+z、往上長）、點值 B·n+ 內插算列 → 依 target 選點（自激發>0 的候選中挑相對偏差最接近 target%）；
%   ②精算（6 點各底面中心、1000 點圓柱平均、build_V_matrix + col_sign）→ 相對 Frobenius。輸出 PDF。
% -------------------------------------------------------------------------
    if nargin<1 || isempty(variant),   variant   = 'full_assembly'; end
    if nargin<2,                        target_pct = 50;            end   % [] → 找最大偏差

    %% ---- paths ----
    CAL = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
           'NTU_hexapole\Calibration_using_FEM_modeling'];
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');      % import_ansys_data
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common'); % ansys_path
    addpath(fullfile(CAL,'current_base','code','main_function'));   % NTU mt_constants
    addpath(fullfile(CAL,'voltage_base','code','main_function'));   % build_sensor_geometry/build_V_matrix/emit_tex
    results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';
    here  = fileparts(mfilename('fullpath'));
    base  = fileparts(fileparts(here));                    % .../voltage_base

    %% ---- 載 V_pdf + 慣例 ----
    S = load(fullfile(base,'data',sprintf('calib_D_%s.mat',variant)));
    V_pdf         = S.Vmat_p;                              % 6×6 [mV]，paper 序，自激發全正
    S_hall        = S.S_hall;
    n_uniform     = S.n_uniform;
    apdl_to_paper = S.apdl_to_paper_idx;                   % [1 3 6 5 2 4]
    paper_to_apdl = S.paper_to_apdl;                       % [1 5 2 6 4 3]
    cnst = mt_constants();
    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);
    AXIAL_TOL = cnst.SENSOR_H;
    col_sign = ones(1,6);                                  % 自激發全正：翻下極激發欄（coil 序）
    for j=1:6, if ismember(apdl_to_paper(j),[1 3 6]), col_sign(j)=-1; end; end

    %% ---- 候選偏移（R_valid 圓盤，同 sensing_area_count 網格）----
    r_max = 2.865e-3;  r = 0.15e-3;  R_valid = r_max - r;  dd = 0.05e-3;
    deltas = dd:dd:R_valid;  phis = (0:0.1:360-0.1)*pi/180;
    [DL,PH] = meshgrid(deltas, phis);
    offx = [0; DL(:).*cos(PH(:))];  offy = [0; DL(:).*sin(PH(:))];  Ncand = numel(offx);
    LOCz = 0.12e-3;
    % 搜尋圓「以極中線 (13.15,0) 為圓心」（pole-local），sensor 標稱在 (13.15,−0.99)＝offset (0,−0.99)（比照 mt_constants SX/SY）。
    SXc = 13.15e-3;  SYnom = -0.99e-3;   % 圓心 x、sensor 標稱 y（offset）。offx/offy = 相對圓心之偏移。
    if isempty(target_pct), fprintf('候選 Ncand=%d；目標=最大偏差（自激發>0）\n', Ncand);
    else, fprintf('候選 Ncand=%d；目標=相對偏差最接近 %.0f%%（自激發>0）\n', Ncand, target_pct); end

    %% ---- 載 coil1：建每 sensor 的 Delaunay + pointLocation 權重 ----
    Bn_reg=cell(1,6); bc=cell(1,6); conn=cell(1,6); good=cell(1,6); idx_reg=cell(1,6); Preg=cell(1,6);
    CC=cell(1,6); RI=cell(1,6);
    d1 = import_ansys_data(fullfile(results_root, variant, 'coil1'), 'all', variant);
    X1 = [d1.x d1.y d1.z];  B1 = [d1.bx d1.by d1.bz];
    for i = 1:6
        thi = cnst.pole_angles(i)*pi/180;  Ri = [cos(thi) -sin(thi); sin(thi) cos(thi)];
        cc  = [Ri*[SXc; 0]; sensor_pos(3,i)];              % 搜尋圓心 world（pole-local (13.15,0)）
        CC{i}=cc; RI{i}=Ri;  ni = sensor_n(:,i);
        sel = (hypot(X1(:,1)-cc(1), X1(:,2)-cc(2)) <= r_max) & (abs(X1(:,3)-cc(3)) <= LOCz);
        idx_reg{i} = find(sel);  Preg{i} = X1(sel,:);
        TR = delaunayTriangulation(Preg{i});
        candxy = (Ri * [SXc+offx.'; offy.']).';            % Ncand×2 world（offset 在 pole-local）
        cand = [candxy, cc(3)+zeros(Ncand,1)];
        ti = pointLocation(TR, cand);  g = ~isnan(ti);
        bc{i}=nan(Ncand,4); conn{i}=nan(Ncand,4); good{i}=g;
        bc{i}(g,:)   = cartesianToBarycentric(TR, ti(g), cand(g,:));
        conn{i}(g,:) = TR.ConnectivityList(ti(g),:);
        Bn_reg{i} = zeros(numel(idx_reg{i}),6);  Bn_reg{i}(:,1) = B1(sel,:)*ni;
        fprintf('  sensor P%d：region %d 節點，凸包內候選 %d/%d\n', i, numel(idx_reg{i}), nnz(g), Ncand);
    end
    for kc = 2:6
        dk = import_ansys_data(fullfile(results_root, variant, sprintf('coil%d',kc)), 'all', variant);
        Xk = [dk.x dk.y dk.z];  Bk = [dk.bx dk.by dk.bz];
        for i = 1:6
            assert(max(abs(Xk(idx_reg{i},:) - Preg{i}), [], 'all') < 1e-9, sprintf('coil%d 順序不一致 P%d', kc, i));
            Bn_reg{i}(:,kc) = Bk(idx_reg{i},:) * sensor_n(:,i);
        end
    end

    %% ---- 每 sensor：候選列（自激發全正）→ 依 target 選點（自激發>0）----
    found_pos = zeros(3,6);  info = zeros(6,7);  dev_sel = zeros(6,1);  self_sel = zeros(6,1);
    for i = 1:6
        BnC = nan(Ncand,6);  g = good{i};
        for kc = 1:6
            reg = Bn_reg{i}(:,kc);  val = nan(Ncand,1);
            val(g) = sum(bc{i}(g,:) .* reg(conn{i}(g,:)), 2);  BnC(:,kc) = val;
        end
        ROWc = S_hall * (BnC*1e3);                          % 自激發全正號誌 = raw 讀值 [mV]（exc·col=+1）
        ROWp = ROWc(:, paper_to_apdl);                     % paper 欄
        dev_rel = vecnorm(ROWp - V_pdf(i,:), 2, 2) / norm(V_pdf(i,:)) * 100;   % 相對 %
        selfV   = ROWp(:, i);                               % 自激發（該列自己欄）電壓
        ok = all(~isnan(ROWp),2) & all(sign(ROWp) == sign(V_pdf(i,:)), 2);  % 凸包內 + 整列號誌與 V_pdf 一致（對角正、off-diag 負）
        if isempty(target_pct)
            score = dev_rel;  score(~ok) = -inf;  [~, im] = max(score);          % 最大偏差
        else
            score = abs(dev_rel - target_pct);  score(~ok) = inf;  [~, im] = min(score);  % 最接近 target
        end
        found_pos(:,i) = [RI{i}*[SXc+offx(im); offy(im)]; CC{i}(3)];   % 找到點 world（offset 在 pole-local）
        dev_sel(i) = dev_rel(im);  self_sel(i) = selfV(im);
        info(i,:) = [hypot(offx(im),offy(im))*1e3, mod(atan2(offy(im),offx(im))*180/pi,360), ...
                     offx(im)*1e3, offy(im)*1e3, found_pos(:,i).'*1e3];   % dx/dy = 相對圓心之 pole-local 偏移
        fprintf('  P%d：δ=%.3f mm, φ=%.1f°, 列相對偏差=%.2f%%, 自激發=%.1f mV\n', ...
                i, info(i,1), info(i,2), dev_sel(i), self_sel(i));
    end

    %% ---- 精算 V（6 點各底面中心、1000 點圓柱平均、col_sign 自激發全正）----
    fprintf('精算 V（build_V_matrix，6 個點）...\n');
    V_found = build_V_matrix(results_root, cnst, apdl_to_paper, found_pos, sensor_n, ...
                             S_hall, [], n_uniform, [], AXIAL_TOL, variant);
    V_found = V_found .* col_sign;                          % 自激發全正號誌（與 main.m 一致）
    V_found_p = V_found(:, paper_to_apdl);
    dV = V_found_p - V_pdf;
    errF_abs = norm(dV,'fro');  errF_rel = errF_abs / norm(V_pdf,'fro') * 100;
    allpos = all(diag(V_found_p) > 0);
    fprintf('\n===== 相對 Frobenius = %.3f %% （絕對 %.3f mV）；自激發全正=%d =====\n', errF_rel, errF_abs, allpos);

    %% ---- 兩張圖（空白基底、只畫 sensor1 P1）→ figures/shared/ ----
    pole = {'P1','P2','P3','P4','P5','P6'};
    figdir = fullfile(base,'figures','shared');  if ~exist(figdir,'dir'); mkdir(figdir); end
    if isempty(target_pct), ftag='max dev'; else, ftag=sprintf('~%.0f%%',target_pct); end
    % 圖1：sensor1 原位置（標稱 = 圓心偏移 (0,−0.99)）
    draw_search_fig([], r_max, R_valid, [0, SYnom*1e3], {'P1'}, ...
        fullfile(figdir,'search_pts_initial.png'), 'sensor1 (P1) original position');
    % 圖2：sensor1 找到的新位置
    draw_search_fig([], r_max, R_valid, info(1,3:4), {'P1'}, ...
        fullfile(figdir,'search_pts_found.png'), sprintf('sensor1 (P1) found (%s) position', ftag));
    fprintf('已輸出 %s\\search_pts_{initial,found}.png\n', figdir);

    %% ---- 出 PDF → results/eighteen_param/ ----
    out_dir = fullfile(base,'results','eighteen_param');  if ~exist(out_dir,'dir'); mkdir(out_dir); end
    tex_path = fullfile(out_dir, sprintf('vmat_search_error_%s.tex', variant));
    pdf_path = fullfile(out_dir, sprintf('vmat_search_error_%s.pdf', variant));
    if isempty(target_pct), tgt = 'max deviation'; else, tgt = sprintf('$\\approx%.0f\\%%$ relative deviation', target_pct); end
    T = emit_tex();
    fid = fopen(tex_path,'w');
    fprintf(fid,'%% Auto-generated by vmat_search_error.m\n');
    fprintf(fid,'\\documentclass[11pt]{article}\n\\usepackage[margin=1in]{geometry}\n\\usepackage{amsmath}\n\\begin{document}\n');
    fprintf(fid,'\\section*{NTU voltage\\_base — sensing search-area V-matrix deviation (%s, self-excitation $>0$)}\n', tgt);
    fprintf(fid,'\\noindent Search area: $r_{\\max}=%.3f$ mm, sensor $r=%.2f$ mm, valid-center disk $R=%.3f$ mm; grid $d\\varphi=0.1^\\circ$, $d\\delta=0.05$ mm ($N=%d$). Candidate $=$ cylinder bottom-face center ($n^+=+z$). Sign: self-excitation positive (flip sink P1/P3/P6 excitation column on $G$ and $V$). Each sensor picks the position whose row keeps the SAME sign pattern as $V_{\\mathrm{pdf}}$ (diagonal $>0$, off-diagonals $<0$) and per-row relative deviation closest to the target.\\\\[4pt]\n', ...
            r_max*1e3, r*1e3, R_valid*1e3, Ncand);
    fprintf(fid,'\\noindent\\textbf{Selected sensor positions (one per sensor):}\\\\[2pt]\n');
    fprintf(fid,'\\begin{tabular}{c|ccccccc}\n');
    fprintf(fid,'Sensor & $\\delta$ (mm) & $\\varphi$ ($^\\circ$) & $dx$ (mm) & $dy$ (mm) & world $(x,y,z)$ (mm) & row dev (\\%%) & self V (mV)\\\\\\hline\n');
    for i=1:6
        fprintf(fid,'%s & %.3f & %.1f & %+.3f & %+.3f & $(%.2f,\\,%.2f,\\,%.2f)$ & %.2f & %.1f\\\\\n', ...
                pole{i}, info(i,1), info(i,2), info(i,3), info(i,4), info(i,5), info(i,6), info(i,7), dev_sel(i), self_sel(i));
    end
    fprintf(fid,'\\end{tabular}\\\\[10pt]\n');
    T.mat(fid, 'V_{\mathrm{pdf}}~[\mathrm{mV}]',   V_pdf,     pole, '');
    T.mat(fid, 'V_{\mathrm{found}}~[\mathrm{mV}]', V_found_p, pole, '');
    T.mat(fid, '\Delta V=V_{\mathrm{found}}-V_{\mathrm{pdf}}~[\mathrm{mV}]', dV, pole, '');
    fprintf(fid,'\\[\n\\frac{\\lVert \\Delta V\\rVert_F}{\\lVert V_{\\mathrm{pdf}}\\rVert_F}\\times100\\%% = %.3f\\%%\\qquad(\\lVert \\Delta V\\rVert_F = %.3f~\\mathrm{mV})\n\\]\n', errF_rel, errF_abs);
    fprintf(fid,'\\noindent All selected self-excitation voltages $V_{\\mathrm{found}}(i,i)>0$: %s.\n', mat2str(allpos));
    fprintf(fid,'\\end{document}\n');
    fclose(fid);

    xelatex = 'C:\Users\Kuo\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe';
    old = cd(out_dir);
    [st,so] = system(sprintf('"%s" -interaction=nonstopmode -halt-on-error "%s"', xelatex, tex_path));
    cd(old);
    if st~=0 || ~exist(pdf_path,'file'); fprintf('%s\n',so); error('xelatex 編譯失敗'); end
    for e = {'.tex','.aux','.log','.out'}, f = fullfile(out_dir, sprintf('vmat_search_error_%s%s', variant, e{1})); if exist(f,'file'); delete(f); end, end
    fprintf('已存 %s\n', pdf_path);
end

% ============================================================================
function draw_search_fig(base_pts, r_max, R_valid, marks, labels, outpng, ttl)
% sensing_area 基底（外圓 r_max + 有效中心區 R_valid + 橘色取樣點）+ 圈出 marks（黑圈）。
    tt = linspace(0,2*pi,400);
    fig = figure('Color','w','Position',[100 100 760 800]);  ax = axes(fig);  hold(ax,'on');
    plot(ax, r_max*1e3*cos(tt),   r_max*1e3*sin(tt),   'k-',  'LineWidth',2.5);
    plot(ax, R_valid*1e3*cos(tt), R_valid*1e3*sin(tt), '--',  'Color',[0.2 0.5 0.9], 'LineWidth',1.8);
    if ~isempty(base_pts)
        plot(ax, base_pts(:,1), base_pts(:,2), '.', 'Color',[0.85 0.33 0.10], 'MarkerSize',4);   % 取樣基底
    end
    for k = 1:size(marks,1)
        plot(ax, marks(k,1), marks(k,2), 'o', 'MarkerSize',15, 'MarkerEdgeColor','k', 'LineWidth',2.6, 'MarkerFaceColor','none');
        plot(ax, marks(k,1), marks(k,2), '+', 'Color','k', 'MarkerSize',11, 'LineWidth',1.8);
        if ~isempty(labels), text(ax, marks(k,1)+0.12, marks(k,2)+0.12, labels{k}, 'FontWeight','bold','FontSize',13,'Color','k'); end
    end
    axis(ax,'equal'); box(ax,'on'); grid(ax,'off');
    lim = r_max*1e3*1.08;  xlim(ax,[-lim lim]); ylim(ax,[-lim lim]);
    set(ax,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]);
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end));  yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel(ax,'x (mm)','FontWeight','bold'); ylabel(ax,'y (mm)','FontWeight','bold');
    title(ax, ttl, 'FontWeight','bold','FontSize',14);
    ax.Toolbar.Visible='off';
    exportgraphics(fig, outpng, 'Resolution',150);  close(fig);
end
