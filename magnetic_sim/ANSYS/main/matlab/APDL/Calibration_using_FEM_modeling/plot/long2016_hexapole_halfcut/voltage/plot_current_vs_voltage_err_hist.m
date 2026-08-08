function plot_current_vs_voltage_err_hist()
%% plot_current_vs_voltage_err_hist.m -- 電荷(電流)模型 vs 電壓(sensor)模型 場誤差疊圖
% =========================================================================
% Current_base（電荷/電流模型）：B = S·(^Bĝ_I·K̄·F)  → 電荷 Gc = ^Bĝ_I·K̄·F。
% Voltage_base（電壓/sensor 模型）：B = S·(Ĥ_V·V)     → 電荷 Gv = Ĥ_V·Vmat。
% 兩者數學上恆等：Ĥ_V = Dv·Vᵀ(VVᵀ)⁻¹ ⇒ Ĥ_V·V = Dv = Gc（6×6 可逆 Vmat）→ 場/誤差/CV 完全相同
% （電壓重建 ≡ 電流重建），疊圖必完全重合。
%
% [MODIFIED 2026-07-14] VARIANT='no_gap'。no_gap 無 Hall-sensor 校正（無 calib_D、無 refined
%   sensor tet CSV），故 Vmat 改「自算」：用 no_gap coilN 'all' 空氣節點的局部 Delaunay，於各
%   sensor 圓柱內均勻撒 1000 點重心內插 B、⟨B·n+⟩×S_hall 平均（同 extract_Vmat_interp 取樣慣例，
%   只是內插來源＝air 節點而非 refined tets）。legend = mean + CV(%)。
% 誤差 = 逐點逐激發 |B_model−B_FEM|（向量差大小，mT）；bin 0.005；選項①粗體框圖。
% =========================================================================
    clc;
    VARIANT='no_gap'; R_select=150e-6; ell0=0.5e-3;   % ell0 [m]（fit_bias 在 SI）
    CAL = ['G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\' ...
           'long2016_hexapole_halfcut\Calibration_using_FEM_modeling'];
    OUTDIR = fullfile(CAL,'voltage_base','figures','shared');
    results_root = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data';
    % [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live config。 原註：mt_constants/import_ansys_data/filter_iron_nodes
    CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
    addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');
    addpath(fullfile(CAL,'current_base','code','main_function'));               % load_coils_actuator/select_ball/fitting/build_S_matrix/solve_KI_bar_gain
    addpath(fullfile(CAL,'voltage_base','code','main_function')); % build_sensor_geometry
    model='long2016_hexapole_halfcut'; cnst=mt_constants(); apdl_to_paper_idx=[1,3,6,5,2,4];
    if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end
    S_hall = 130;   % mV/mT（尺度不影響誤差：Ĥ_V 會重標；僅為 Vmat 單位一致）

    % ---- no_fix bias 擬合（一次）→ A, Bstack, Dv ----
    D = load_coils_actuator(model, cnst, apdl_to_paper_idx, 'all', VARIANT);
    [P, Bstack, npts] = select_ball(D, R_select);
    s_sink=ones(1,6); for j=1:6, if ismember(apdl_to_paper_idx(j),[1 3 6]), s_sink(j)=-1; end; end
    Bstack = (-Bstack).*s_sink;     % flip-sink all-source（Bstack 已 mT）
    [ell, e_hat, Pc] = fitting(P, Bstack, D.Pc_base, ell0, true);   % 18-param bias；回 (ℓ̄,ê,Pc)
    A  = build_S_matrix(ell, Pc, P);
    Dv = (A.'*A)\(A.'*Bstack);

    % ---- Current_base：Gc = ^Bĝ_I·K̄·F ----
    [K_bar, ghat_I_B] = solve_KI_bar_gain(ell, Pc, P, Bstack, D.F);
    Gc = ghat_I_B * K_bar * D.F;
    err_c = vecnorm(reshape(A*Gc - Bstack, 3, npts*6), 2, 1).';

    % ---- Voltage_base：自算 no_gap Vmat（air 節點內插）→ Ĥ_V → Gv ----
    [sensor_pos, sensor_n] = build_sensor_geometry(cnst);
    [Vmat, ngood] = vmat_nogap_interp(results_root, cnst, apdl_to_paper_idx, sensor_pos, sensor_n, S_hall, VARIANT);
    fprintf('Vmat 內插命中：min %d / max %d / 1000（每 sensor×激發）\n', min(ngood(:)), max(ngood(:)));
    H_V = Dv * Vmat.' / (Vmat * Vmat.');                   % Ĥ_V = Dv·Vᵀ(VVᵀ)⁻¹
    Gv  = H_V * Vmat;                                       % = Dv（恆等）
    err_v = vecnorm(reshape(A*Gv - Bstack, 3, npts*6), 2, 1).';

    fprintf('自檢：max|Gc−Gv|=%.2e ；max|Gc−Dv|=%.2e ；max|Gv−Dv|=%.2e\n', ...
            max(abs(Gc(:)-Gv(:))), max(abs(Gc(:)-Dv(:))), max(abs(Gv(:)-Dv(:))));
    ec=err_c; ev=err_v;
    cvc=100*std(ec)/mean(ec); cvv=100*std(ev)/mean(ev);
    fprintf('Current_base: mean=%.4f mT CV=%.1f%% ；Voltage_base: mean=%.4f mT CV=%.1f%%\n', ...
            mean(ec), cvc, mean(ev), cvv);
    Jc = sum(ec.^2); Jv = sum(ev.^2);   % cost J = ‖A·G−Bstack‖² [mT²]
    fprintf('COST J [mT^2]: current_base=%.6e | hall(voltage)_base=%.6e | |diff|=%.3e (relative %.2e)\n', ...
            Jc, Jv, abs(Jc-Jv), abs(Jc-Jv)/Jc);

    %% ---- 疊圖（選項①粗體框圖、count、bin 0.005、legend mean+CV%）----
    XMAX=0.8; edges=0:0.005:XMAX;
    fig=figure('Color','w','Position',[100 100 820 580]); ax=axes(fig); hold(ax,'on');
    histogram(ax, ec, edges, 'FaceColor',[0.85 0.33 0.10], 'FaceAlpha',0.55, 'EdgeColor','w');
    histogram(ax, ev, edges, 'FaceColor',[0.20 0.40 0.70], 'FaceAlpha',0.55, 'EdgeColor','w');
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]); box(ax,'on'); grid(ax,'off');
    xlim(ax,[0 XMAX]); xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end)); yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel(ax,'|B_{model} - B_{FEM}| (mT)','FontWeight','bold'); ylabel(ax,'Count','FontWeight','bold');
    lc=sprintf('Current\\_base  (mean %.3f mT, CV %.1f%%)', mean(ec), cvc);
    lv=sprintf('Voltage\\_base  (mean %.3f mT, CV %.1f%%)', mean(ev), cvv);
    lg=legend(ax,{lc,lv},'Location','northeast','Interpreter','tex'); set(lg,'FontSize',14,'FontWeight','bold','Box','on');

    out=fullfile(OUTDIR,sprintf('current_vs_voltage_err_hist_%s.png',VARIANT));
    exportgraphics(fig,out,'Resolution',600); fprintf('已輸出 %s\n', out);
end

% =========================================================================
function [Vmat, ngood] = vmat_nogap_interp(results_root, cnst, idx, sensor_pos, sensor_n, S_hall, variant)
% no_gap Vmat：coilN 'all' 空氣節點局部 Delaunay，各 sensor 圓柱內均勻 1000 點重心內插 → ⟨B·n+⟩×S_hall。
% 圓柱 R=0.15mm × H=0.10mm 沿 n+；鄰域 R_nbr=1.0mm（涵蓋取樣圓柱、含足夠 air 節點）。all-source 翻下極激發欄。
    n_uniform=1000; sensor_r=0.15e-3; axial_tol=0.10e-3; R_nbr=1.0e-3;
    rng(0);
    samp = cell(1,6);
    for i=1:6
        ci = sensor_pos(:,i) + [0;0;cnst.SPH_OFST];       % sensor 中心（ANSYS 框）
        ni = sensor_n(:,i);
        t1 = [-ni(2); ni(1); 0]; if norm(t1)<1e-9, t1=[1;0;0]; end
        t1 = t1/norm(t1); t2 = cross(ni,t1);              % 圓柱橫切基底 ⊥ n+
        a  = axial_tol*rand(n_uniform,1);
        r  = sensor_r*sqrt(rand(n_uniform,1));            % √U → 面積均勻
        th = 2*pi*rand(n_uniform,1);
        samp{i} = ci.' + a.*ni.' + (r.*cos(th)).*t1.' + (r.*sin(th)).*t2.';
    end
    Vmat = zeros(6,6); ngood = zeros(6,6);
    for kc = 1:6
        cn = sprintf('coil%d', kc);
        ds = import_ansys_data(fullfile(results_root, variant, cn),'all',cn);
        air = filter_iron_nodes(ds.x, ds.y, ds.z, cnst, struct('visualize',false));   % 只留空氣節點（避免跨鋼/空氣界面內插）
        Xa = [ds.x(air), ds.y(air), ds.z(air)];
        Ba = 1e3*[ds.bx(air), ds.by(air), ds.bz(air)];    % Tesla → mT
        for i = 1:6
            ci = (sensor_pos(:,i) + [0;0;cnst.SPH_OFST]).';  ni = sensor_n(:,i);
            loc = find(vecnorm(Xa - ci, 2, 2) < R_nbr);
            TR  = delaunayTriangulation(Xa(loc,:));
            ti  = pointLocation(TR, samp{i});  good = ~isnan(ti);
            bc  = cartesianToBarycentric(TR, ti(good), samp{i}(good,:));
            conn= TR.ConnectivityList(ti(good),:);
            Bp  = zeros(nnz(good),3);
            for c = 1:4, Bp = Bp + bc(:,c).*Ba(loc(conn(:,c)),:); end
            Vmat(i,kc) = S_hall * mean(Bp*ni);            % ⟨B·n+⟩[mT]×S_hall[mV/mT]
            ngood(i,kc)= nnz(good);
        end
    end
    exc = ones(1,6); for j=1:6, if ismember(idx(j),[1 3 6]), exc(j)=-1; end; end   % all-source 翻下極激發欄
    Vmat = Vmat .* exc;
end
