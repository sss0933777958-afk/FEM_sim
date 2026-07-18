function plot_circle_force_voltage_18p()
% plot_circle_force_voltage_18p — 電壓基礎力模型(18-param) 沿 R=50µm 圓的力。
%   F⃗ = ᶠĝ_V·( Vᵀ D̄ᵀ [Lx;Ly;Lz] D̄ V )，ᶠĝ_V=1.875e-7；Lx=∂/∂x(SᵀS)，S_i(p̄)=(p̄−Pc_i)/|p̄−Pc_i|³。
%   = ᶠĝ_V·∇(|B̃|²)，B̃=Σ_i q_i S_i(p̄)，q=D̄·V（coil1@1A 的 sensor 電壓）。
%   資料：gap_calibrate 18-param（calib_D_gap_calibrate_bias.mat）：D̄、Vmat、Pc（actuator 正規化）、ℓ̂、R_act。
%   圓在 measure/WP frame(z=0)，轉 actuator 算 p̄=R_act·P/ℓ̂，力直接為 measure 分量（∇ 對 measure 座標）。
%   解析電荷模型 → 平滑力（對照 current_base 的 raw FEM |B|² 力）。單位：力 pN、長度 µm。
%   輸出：voltage_base/figures/eighteen_param/circle_force_voltage_coil1.png。

    here  = fileparts(mfilename('fullpath'));
    vbase = fileparts(fileparts(here));
    figdir = fullfile(vbase,'figures','eighteen_param'); if ~exist(figdir,'dir'); mkdir(figdir); end
    outpng = fullfile(figdir,'circle_force_voltage_coil1.png');

    S = load(fullfile(vbase,'data','calib_D_gap_calibrate_bias.mat'));
    Dbar = S.D_bar;  Vmat = S.Vmat_p;  Pc = S.Pc;  ell = S.ell_hat;  R_act = S.R_act;  % ell µm; Pc 3x6 norm
    V = Vmat(:,1);                       % coil1 = P1 激發 sensor 電壓 [mV]
    q = Dbar * V;                        % 6x1 電荷向量（Pc 欄序）
    fprintf('ℓ̂=%.1f µm | q=D̄V = [%s]\n', ell, sprintf('%.3g ', q));

    % R=50µm 圓（measure/WP frame, z=0）
    Nn = 720; th = linspace(0,2*pi,Nn).'; Rc = 50;
    Pm = [Rc*cos(th), Rc*sin(th), zeros(Nn,1)];   % µm, measure frame

    % Φ(P)=|Σ q_i S_i(p̄)|²，p̄=R_act·P/ℓ̂；F=1.875e-7·∇_P Φ（中央差分，δ µm → measure 分量）
    PRE = 1.875e-7; dl = 0.5;
    ex=[dl 0 0]; ey=[0 dl 0]; ez=[0 0 dl];
    Fx = PRE*(phi(Pm+ex,R_act,ell,Pc,q) - phi(Pm-ex,R_act,ell,Pc,q))/(2*dl);
    Fy = PRE*(phi(Pm+ey,R_act,ell,Pc,q) - phi(Pm-ey,R_act,ell,Pc,q))/(2*dl);
    Fz = PRE*(phi(Pm+ez,R_act,ell,Pc,q) - phi(Pm-ez,R_act,ell,Pc,q))/(2*dl);
    idx=(1:Nn).';
    fprintf('|F|: Fx[%.3g,%.3g] Fy[%.3g,%.3g] Fz[%.3g,%.3g] pN\n', ...
        min(Fx),max(Fx),min(Fy),max(Fy),min(Fz),max(Fz));

    f=figure('Color','w','Position',[100 100 1150 900]);
    lbl={'F_x','F_y','F_z'}; Fc={Fx,Fy,Fz};
    for k=1:3
        subplot(3,1,k); plot(idx,Fc{k},'k-','LineWidth',1.4); grid on;
        set(gca,'FontSize',14,'FontWeight','bold','LineWidth',1.5);
        ylabel([lbl{k} ' (pN)'],'FontWeight','bold'); xlim([1 Nn]);
        if k==1; title('coil1 @ 1A    voltage-base 18-param force on R=50 \mum circle','FontWeight','bold'); end
        if k==3; xlabel('point index (along R=50\mum circle)','FontWeight','bold'); end
    end
    exportgraphics(f,outpng,'Resolution',150); close(f);
    fprintf('wrote %s\n', outpng);
end

function v = phi(Pmat, R_act, ell, Pc, q)
% Φ = |Σ_i q_i S_i(p̄)|²，p̄=R_act·P/ℓ̂（Pmat µm, measure frame）
    pb = (R_act * Pmat.') / ell;         % 3xN 正規化 actuator
    N = size(Pmat,1); Btil = zeros(3,N);
    for i=1:6
        Di = pb - Pc(:,i);               % 3xN
        r3 = sum(Di.^2,1).^1.5;          % 1xN
        Btil = Btil + q(i)*(Di./r3);     % 3xN
    end
    v = sum(Btil.^2,1).';                % Nx1
end
