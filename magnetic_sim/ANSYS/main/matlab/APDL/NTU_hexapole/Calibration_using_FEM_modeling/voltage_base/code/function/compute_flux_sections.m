function T = compute_flux_sections(H, MODE, ZCUT_IN)
%COMPUTE_FLUX_SECTIONS  NTU：導柱 / 磁極板兩個截面的磁通量 Φ（coil1 = P1 自激發）。
%   Φ = Σ (B·n̂)·ΔA —— 截面切成邊長 H 的微小正方形，**每格中心一點**，該點的 B 由周圍
%   FEM 節點做**重心內插**（Delaunay + pointLocation + cartesianToBarycentric），再乘格面積加總。
%
%   兩個截面（法線都是軸向 → B·n̂ 退化成取單一分量）：
%     ① 導柱 @ z = 15.00：圓盤 r2.5 @ (20.5,0)、n̂ = +z → Φ = Σ B_z·ΔA
%     ② 磁極 @ x = 13.15：矩形 y∈[-w,+w] × z∈[8.80,9.05]、n̂ = +x → Φ = Σ B_x·ΔA
%        （w 由 9 段輪廓內插得 2.9942mm，不寫死）
%   z=15 只有導柱（襯套 flange 從 15.60 起）→ 該截面兩 variant 幾何相同。
%
%   ⚠ 三個踩過的坑（本 session 實測，改動前先看）：
%     1. **格子必須鋪滿**：用 `z0+H/2:H:z1` 在 H=0.10 會給 3 排 = 0.30mm（板厚只有 0.25）→ 面積多算
%        +20%。故用 N = round(L/H)、實際格距 L/N。
%     2. **Delaunay 只能用鋼節點建**：B 在鋼/空氣界面不連續（B_x 對板面是切向 → 差 μr 倍），
%        混著建會被空氣值拉低且隨格距抖動（板僅 0.25mm 厚、整片都在界面附近，受害最重）。
%     3. **收斂已驗**：H 0.05 → 0.01（格數 ×25）Φ 只變 0.3%（導柱）/ 0.8%（磁極）→ 0.05 足夠。
%
%   [ADDED] MODE = 'sections'（預設，上述兩截面）| 'box'（襯套 flange 段 Gauss 盒三塊，見 §box）
%
%   §box —— 繞 P1 導柱軸 (20.5,0) 的 Gauss 盒（範圍圖：code/plot/plot_flux_box.m）：
%     ① 底面 z=15.50、r<=2.5      鋼（導柱），n̂=+z → Σ B_z·ΔA
%     ② 底面 z=15.50、2.5<r<=7.5  空氣，        n̂=+z → Σ B_z·ΔA
%     ③ 側壁 r=7.52、z=15.60~16.60 空氣（貼 flange 外緣），n̂=r̂(θ) → Σ (B_x cosθ + B_y sinθ)·ΔA
%   ⚠ 頂面 z=16.60 未含 → 這三塊**尚未封閉**，不可期待 ①+②+③ = 0。
%   ⚠ 取樣面刻意離開鋼/空氣界面（底面 15.50 而非 flange 底 15.60、側壁 7.52 而非 7.50）：
%      界面上 B 切向不連續、FEM 節點只有一份，貼著取會被另一側材料汙染。偏移量 ~0.02–0.10mm
%      « 該區 tet 尺度，對 Φ 影響遠小於 1%（同一 tet 內線性）。
%   ⚠ ①② 與 §sections 的導柱截面**不同 z**（15.50 vs 15.00）→ 值本來就不會相同（漏磁隨 z 變）。
%
%   T = compute_flux_sections(H, MODE)  H 預設 0.05 [mm]；回 table（variant × section × Φ[µWb] × 面積 × 格數）

    if nargin < 1 || isempty(H), H = 0.05; end
    if nargin < 2 || isempty(MODE), MODE = 'sections'; end
    if nargin < 3 || isempty(ZCUT_IN), ZCUT_IN = 15.00; end   % [ADDED] 導柱截面 z（可換，見下）

    here  = fileparts(mfilename('fullpath'));
    vbase = fileparts(fileparts(here));
    addpath(here);                                                       % ntu_pole_profile
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % import_ansys_data
    droot = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\NTU_hexapole\data';

    % ---- 截面幾何 ----
    ZCUT = ZCUT_IN;  PC = [20.5, 0];  PR = 2.5;        % ① 導柱圓盤（z 由參數給，預設 15.00）
    XCUT = 13.15;  ZB = 8.80;  ZT = 9.05;              % ② 磁極矩形
    prof = ntu_pole_profile();
    lo = prof(prof(:,2)<0,:);  [~,i] = sort(lo(:,1));  lo = lo(i,:);
    W  = abs(interp1(lo(:,1), lo(:,2), XCUT));         % 板半寬 @ XCUT
    TOL = 1e-6;

    variants = {'full_assembly_sleeve','full_assembly'};
    vlab     = {'sleeve mu_r=280','sleeve mu_r=1 (none)'};

    if strcmpi(MODE,'box')
        T = run_box(H, droot, variants, vlab, PC, PR);
        return;
    end

    rows = {};
    fprintf('\n格距 H = %.3f mm\n', H);
    fprintf('%-22s | %-24s | %12s | %10s | %8s\n', 'variant','section','Phi [uWb]','A [mm^2]','cells');
    fprintf('%s\n', repmat('-',1,88));
    for v = 1:2
        d = import_ansys_data(fullfile(droot,variants{v},'coil1'),'all',variants{v});
        X = [d.x d.y d.z]*1e3;  B3 = [d.bx d.by d.bz]*1e3;               % mm / mT
        % 鋼件遮罩（★ Delaunay 只用鋼節點）
        isPost  = hypot(X(:,1)-PC(1), X(:,2)-PC(2)) <= PR+TOL & X(:,3)>=9.05-TOL & X(:,3)<=26.80+TOL;
        isPlate = inpolygon(X(:,1),X(:,2),prof(:,1),prof(:,2)) & X(:,3)>=ZB-TOL & X(:,3)<=ZT+TOL;

        % ---- ① 導柱 @ z=ZCUT，n̂=+z：格心落在 r<PR 內（圓邊界誤差 O(H)）----
        sp = isPost & abs(X(:,3)-ZCUT) < 1.5;
        TR = delaunayTriangulation(X(sp,:));  BP = B3(sp,:);
        g = -PR+H/2 : H : PR;  [GX,GY] = ndgrid(g,g);  in = GX.^2+GY.^2 <= PR^2;
        P = [PC(1)+GX(in), PC(2)+GY(in), ZCUT*ones(nnz(in),1)];
        [phi1, n1] = integrate(TR, BP, P, 3, H*H);                        % 取 B 的第 3 分量 = B·(+z)
        A1 = nnz(in)*H*H;
        rows(end+1,:) = {vlab{v}, sprintf('post  @ z=%.2f (n=+z)',ZCUT), phi1, A1, n1};  %#ok<AGROW>
        fprintf('%-22s | %-24s | %+12.4f | %10.4f | %8d\n', vlab{v}, sprintf('post  @ z=%.2f',ZCUT), phi1, A1, n1);

        % ---- ② 磁極 @ x=XCUT，n̂=+x：★ 強制整除鋪滿 ----
        sq = isPlate & abs(X(:,1)-XCUT) < 1.5;
        TR2 = delaunayTriangulation(X(sq,:));  BP2 = B3(sq,:);
        Nz = round((ZT-ZB)/H);  dz = (ZT-ZB)/Nz;
        Ny = round(2*W/H);      dy = 2*W/Ny;
        [QY,QZ] = ndgrid(-W+dy/2 : dy : W-dy/2,  ZB+dz/2 : dz : ZT-dz/2);
        P2 = [XCUT*ones(numel(QY),1), QY(:), QZ(:)];
        [phi2, n2] = integrate(TR2, BP2, P2, 1, dy*dz);                   % 取第 1 分量 = B·(+x)
        A2 = Ny*Nz*dy*dz;
        rows(end+1,:) = {vlab{v}, sprintf('pole  @ x=%.2f (n=+x)',XCUT), phi2, A2, n2};  %#ok<AGROW>
        fprintf('%-22s | %-24s | %+12.4f | %10.4f | %8d\n', vlab{v}, sprintf('pole  @ x=%.2f',XCUT), phi2, A2, n2);
    end
    fprintf('%s\n', repmat('-',1,88));
    fprintf('理論面積：post %.4f mm^2 / pole %.4f mm^2\n', pi*PR^2, 2*W*(ZT-ZB));
    T = cell2table(rows, 'VariableNames', {'variant','section','Phi_uWb','A_mm2','cells'});
end

% ============================================================================
function T = run_box(H, droot, variants, vlab, PC, PR)
% §box：flange 段 Gauss 盒三塊（①② 底面 n̂=+z / ③ 側壁 n̂=r̂）
    R_FL = 7.5;  ZB = 15.60;  ZM = 16.60;
    ZDISC = 15.50;          % 底面取樣 z（往空氣側離開 flange 底 15.60 共 0.10mm）
    R_SIDE = R_FL + 0.02;   % 側壁取樣 r（往空氣側離開 flange 外緣 7.50 共 0.02mm）
    TOL = 1e-6;

    rows = {};
    fprintf('\n格距 H = %.3f mm | Gauss box @ post axis (%.1f, %.1f)\n', H, PC);
    fprintf('%-22s | %-30s | %12s | %10s | %8s\n','variant','face','Phi [uWb]','A [mm^2]','cells');
    fprintf('%s\n', repmat('-',1,94));
    for v = 1:2
        d = import_ansys_data(fullfile(droot,variants{v},'coil1'),'all',variants{v});
        X = [d.x d.y d.z]*1e3;  B3 = [d.bx d.by d.bz]*1e3;             % mm / mT
        rr = hypot(X(:,1)-PC(1), X(:,2)-PC(2));

        % ---- 共用：底面笛卡兒格心（同 §sections 的分格法，值可對照）----
        g = -R_FL+H/2 : H : R_FL;  [GX,GY] = ndgrid(g,g);  RG = hypot(GX,GY);
        in1 = RG <= PR;              % ① 導柱段
        in2 = RG > PR & RG <= R_FL;  % ② 空氣段

        % ---- ① 導柱鋼節點（z=ZDISC 在導柱 9.05~26.80 內部，無界面問題）----
        sp = rr <= PR+TOL & X(:,3)>=9.05-TOL & X(:,3)<=26.80+TOL & abs(X(:,3)-ZDISC)<1.5;
        P1 = [PC(1)+GX(in1), PC(2)+GY(in1), ZDISC*ones(nnz(in1),1)];
        [p1, n1] = integrate(delaunayTriangulation(X(sp,:)), B3(sp,:), P1, 3, H*H);
        rows(end+1,:) = {vlab{v}, sprintf('(1) post  z=%.2f r<=%.1f',ZDISC,PR), p1, nnz(in1)*H*H, n1};   %#ok<AGROW>

        % ---- ② 空氣節點（r>2.5、z<15.60 → 該區全空氣，兩 variant 皆然）----
        sa = rr > PR+TOL & rr < R_FL+1.0 & X(:,3)>=13.5 & X(:,3)<=ZB-TOL;
        P2 = [PC(1)+GX(in2), PC(2)+GY(in2), ZDISC*ones(nnz(in2),1)];
        [p2, n2] = integrate(delaunayTriangulation(X(sa,:)), B3(sa,:), P2, 3, H*H);
        rows(end+1,:) = {vlab{v}, sprintf('(2) air   z=%.2f %.1f<r<=%.1f',ZDISC,PR,R_FL), p2, nnz(in2)*H*H, n2}; %#ok<AGROW>

        % ---- ③ 側壁：極座標鋪滿（弧長 = H），n̂ = r̂(θ) 逐格不同 ----
        Nz  = round((ZM-ZB)/H);        dz = (ZM-ZB)/Nz;
        Nth = round(2*pi*R_FL/H);      dth = 2*pi/Nth;
        [TH,ZZ] = ndgrid(dth/2 : dth : 2*pi-dth/2,  ZB+dz/2 : dz : ZM-dz/2);
        th = TH(:);
        P3 = [PC(1)+R_SIDE*cos(th), PC(2)+R_SIDE*sin(th), ZZ(:)];
        ss = rr > R_FL-TOL & rr < R_FL+1.5 & X(:,3)>=ZB-1.0 & X(:,3)<=ZM+1.0;   % flange 外的空氣
        dA3 = R_FL*dth*dz;
        [p3, n3] = integrate_radial(delaunayTriangulation(X(ss,:)), B3(ss,:), P3, th, dA3);
        rows(end+1,:) = {vlab{v}, sprintf('(3) side  r=%.2f z=%.2f~%.2f',R_SIDE,ZB,ZM), p3, Nth*Nz*dA3, n3}; %#ok<AGROW>

        for k = size(rows,1)-2 : size(rows,1)
            fprintf('%-22s | %-30s | %+12.4f | %10.4f | %8d\n', rows{k,1}, rows{k,2}, rows{k,3}, rows{k,4}, rows{k,5});
        end
        fprintf('%-22s | %-30s | %+12.4f |            |\n', vlab{v}, '(1)+(2) disc r<=7.5', p1+p2);
        fprintf('%-22s | %-30s | %+12.4f |   <-- NOT closed (no top face)\n', vlab{v}, '(1)+(2)+(3)', p1+p2+p3);
        fprintf('%s\n', repmat('-',1,94));
    end
    fprintf('理論面積：(1) %.4f / (2) %.4f / (3) %.4f mm^2\n', pi*PR^2, pi*(R_FL^2-PR^2), 2*pi*R_FL*(ZM-ZB));
    T = cell2table(rows, 'VariableNames', {'variant','section','Phi_uWb','A_mm2','cells'});
end

% ============================================================================
function [phi, n] = integrate_radial(TR, BP, P, th, dA)
% 側壁專用：法線 n̂ = r̂(θ) = (cosθ, sinθ, 0) 逐格不同 → B·n̂ = Bx·cosθ + By·sinθ
    ti = pointLocation(TR, P);  ok = ~isnan(ti);
    bc = cartesianToBarycentric(TR, ti(ok), P(ok,:));
    cn = TR.ConnectivityList(ti(ok),:);
    Bi = zeros(nnz(ok),3);
    for k = 1:4, Bi = Bi + bc(:,k).*BP(cn(:,k),:); end
    Bn  = Bi(:,1).*cos(th(ok)) + Bi(:,2).*sin(th(ok));
    phi = sum(Bn) * dA * 1e-3;              % mT·mm² → µWb
    n   = nnz(ok);
end

% ============================================================================
function [phi, n] = integrate(TR, BP, P, comp, dA)
% 每格中心一點 → 重心內插 B → 取第 comp 分量（= B·n̂，法線為軸向）→ ×格面積加總。
%   單位：B [mT] × A [mm^2] = 1e-3 T × 1e-6 m^2 = 1e-9 Wb = 1e-3 µWb
    ti = pointLocation(TR, P);  ok = ~isnan(ti);
    bc = cartesianToBarycentric(TR, ti(ok), P(ok,:));
    cn = TR.ConnectivityList(ti(ok),:);
    Bi = zeros(nnz(ok),3);
    for k = 1:4, Bi = Bi + bc(:,k).*BP(cn(:,k),:); end
    phi = sum(Bi(:,comp)) * dA * 1e-3;      % → µWb
    n   = nnz(ok);
end
