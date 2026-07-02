function plot_ellipsoids_yaxis_3d()
% PLOT_ELLIPSOIDS_YAXIS_3D  沿 y 軸的致動橢球 + U 主軸（fix-ℓ, ref 平面 z=0）。
%   從 ref 平面沿 y 軸取 7 點 (0,y,0)（y=−150:50:150µm），每點對 T=S_i(p)·Ĥ_I 做 SVD：
%     致動橢球 = 單位電流球經 T 的像，半軸 σ_k 沿 U(:,k)（U=場端主軸）。
%   畫：每點橢球（半透明）+ U 三個單位向量箭頭（k=1 紅/2 綠/3 藍）。只留橢球+U（無球/圓/線）。
%   ★ model-derived。橢球需 daspect([1 1 1]) 才不變形。輸出 fix_dir/figures/ellipsoids_yaxis_3d.png。

    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');  % mt_constants
    cnst = mt_constants();
    here   = fileparts(mfilename('fullpath'));
    fixdir = fileparts(fileparts(here));
    figdir = fullfile(fixdir,'figures');   if ~exist(figdir,'dir'); mkdir(figdir); end

    S = load(fullfile(fixdir,'data','fit_fixl_R150um_gap200um_mueq.mat'),'ell','gB','Khat');
    ell_m  = S.ell*1e-6;   Hhat_I = S.gB*S.Khat;
    tip  = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];
    dhat = tip ./ vecnorm(tip);

    yv = -150:50:150;   ny = numel(yv);                            % 7 點 [µm]
    Uc = cell(1,ny);  SVc = cell(1,ny);
    for i = 1:ny
        p  = [0; yv(i); 0]*1e-6;
        Dk = p/ell_m - dhat;
        [U,Sig,~] = svd((Dk ./ (vecnorm(Dk).^3)) * Hhat_I);
        Uc{i} = U;  SVc{i} = diag(Sig);
    end
    s1max = max(cellfun(@(s) s(1), SVc));
    sc    = 20 / s1max;                                            % 最大半軸 ≈ 20µm
    Lu    = 26;                                                    % U 箭頭長 [µm]
    col   = [0.85 0.10 0.10; 0.10 0.60 0.20; 0.10 0.35 0.85];     % k=1/2/3 紅/綠/藍

    [ex,ey,ez] = sphere(24);  sz0 = size(ex);
    fig = figure('Color','w','Position',[60 60 1400 520]); hold on;
    for i = 1:ny
        c = [0; yv(i); 0];  U = Uc{i};  sv = SVc{i};
        %% 橢球（半軸 σ_k·sc 沿 U）
        E = U * diag(sv*sc) * [ex(:) ey(:) ez(:)].';              % 3×Npts
        surf(reshape(E(1,:),sz0)+c(1), reshape(E(2,:),sz0)+c(2), reshape(E(3,:),sz0)+c(3), ...
             'FaceColor',[0.45 0.62 0.88], 'FaceAlpha',0.30, 'EdgeColor','none');
        %% U 三個單位向量
        for k = 1:3
            quiver3(c(1),c(2),c(3), U(1,k),U(2,k),U(3,k), Lu, ...
                    'Color',col(k,:), 'LineWidth',2.5, 'MaxHeadSize',0.6);
        end
    end
    plot3(0,0,0,'k+','MarkerSize',12,'LineWidth',2.0);

    %% 框 / 風格（box on + daspect：橢球不變形）
    grid off; box on; daspect([1 1 1]);
    xlim([-32 32]); ylim([-175 175]); zlim([-32 32]);
    view(120, 25);
    set(gca,'FontSize',15,'FontWeight','bold','LineWidth',2,'TickLength',[.012 .012]);
    set(gca,'XTick',-30:30:30,'YTick',-150:75:150,'ZTick',-30:30:30);
    xlabel('x_m (\mum)','FontWeight','bold');
    ylabel('y_m (\mum)','FontWeight','bold');
    zlabel('z_m (\mum)','FontWeight','bold');
    ax = gca; ax.Toolbar.Visible = 'off';
    hold off;

    out = fullfile(figdir,'ellipsoids_yaxis_3d.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('ny=%d, sc=%.3f µm/(mT/A)；saved %s\n', ny, sc, out);
    close(fig);
end
