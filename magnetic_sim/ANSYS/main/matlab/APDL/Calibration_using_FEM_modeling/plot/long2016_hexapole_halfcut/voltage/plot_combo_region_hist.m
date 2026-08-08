%% plot_combo_region_hist.m -- 50 組任意 signed 電流的 voltage_base 相對預測誤差直方圖 + mean/CV
% 情境：sensor 貼好、voltage_base 模型(Ĥ_V)校正好後，對不同激發電流 I 的 WP 場預測能力。
% voltage 預測場 = S·Ĥ_V·V（Ĥ_V=[mT/mV]、V=sensor 電壓向量=V_mat·I）；self-consistent Ĥ_V·V≡G·I
%   → 數值上 S·Ĥ_V·V = S·G·I，故誤差 = ‖res·I‖，res=S·G−b。
% x = 相對預測誤差 = ‖S·Ĥ_V·V − b·I‖ / ‖b·I‖  (%)（誤差佔實際場的百分比；與電流大小無關）。
% y = count/50×100%。標示 50 組的 mean 與 CV=std/mean（不用外部 ref）。
% 50 組 signed 電流 {±1,±2,±3}^6。SET=1→rng(7) 原始 50 組；SET=2→排除 set1 後另抽 50 組（rng8，與 set1 不重複）。風格①、nb=18。
% 輸出 figures/eighteen_param/combo_region_hist.png（SET1）或 combo_region_hist_set2.png（SET2）。
clear; clc;
here=fileparts(mfilename('fullpath')); CAL=fileparts(fileparts(fileparts(here))); MODELDIR=fileparts(CAL);
% [MODIFIED 2026-08-08] 脫離 backup（規則 no-backup-data）→ live config。
CALROOT = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\APDL\Calibration_using_FEM_modeling';
addpath(fullfile(CALROOT,'function'), fullfile(CALROOT,'utils'), fullfile(CALROOT,'common_path'));
addpath(fullfile(MODELDIR,'common'));
addpath(fullfile(CAL,'current_base','code','main_function'));
addpath(fullfile(CAL,'voltage_base','code','main_function'));
cnst=mt_constants(); map=[1,3,6,5,2,4]; R=150e-6;
S=load(fullfile(CAL,'voltage_base','data','calib_D_graded.mat'));
D=load_coils_actuator('long2016_hexapole_halfcut',cnst,map,'all','graded');
[P,b,~]=select_ball(D,R); ss=ones(1,6); for j=1:6, if ismember(map(j),[1 3 6]),ss(j)=-1;end; end
b=(-b).*ss; A=build_S_matrix(S.ell_hat*1e-6,S.Pc,P); G=(A'*A)\(A'*b); res=A*G-b;   % res=S·G−b (all-source)

SET=0;    % 0=全 46656 組 {±1,±2,±3}^6 枚舉；1=原始50組(rng7)；2=排除 set1 另抽50組(rng8)；12=兩批合併100組
lv=[-3 -2 -1 1 2 3];
if SET==0
    % [-3,+3] 連續實數電流「所有組合」→ uniform 密集取樣代表全體（paper 序），分塊向量化
    rng(42); Np=200000; pool=-3+6*rand(Np,6);     % 200000×6 任意實數 ∈[-3,3]
    relpct=zeros(Np,1); blk=2000;
    for i0=1:blk:Np
        idx=i0:min(i0+blk-1,Np);
        Ip=(ss'.*pool(idx,map)');                 % 6×n：套 all-source 番號→coil 框
        relpct(idx)=(100*vecnorm(res*Ip)./vecnorm(b*Ip))';   % ‖S·Ĥ_V·V−b·I‖/‖b·I‖ (%)
    end
    MEANv=mean(relpct); CV=std(relpct)/MEANv*100;                               % CV=std/mean
    fprintf('SET 0 | [-3,3] 實數 %d 組 相對預測誤差: min=%.3f mean=%.3f median=%.3f max=%.3f %% | CV=std/mean=%.1f%%\n', ...
            Np, min(relpct),MEANv,median(relpct),max(relpct), CV);            % 太多不逐組列印
else
    rng(7); pool1=zeros(0,6);                                                    % set1
    while size(pool1,1)<50, c=lv(randi(6,1,6)); if ~ismember(c,pool1,'rows'), pool1(end+1,:)=c; end; end
    rng(8); pool2=zeros(0,6);                                                    % set2（排除 set1）
    while size(pool2,1)<50
        c=lv(randi(6,1,6));
        if ~ismember(c,pool1,'rows') && ~ismember(c,pool2,'rows'), pool2(end+1,:)=c; end
    end
    switch SET, case 1, pool=pool1; case 2, pool=pool2; otherwise, pool=[pool1;pool2]; end
    Np=size(pool,1);
    relpct=zeros(Np,1);
    for n=1:Np, Ip=(ss.*pool(n,map))'; relpct(n)=100*norm(res*Ip)/norm(b*Ip); end   % ‖S·Ĥ_V·V−b·I‖/‖b·I‖ (%)
    MEANv=mean(relpct); CV=std(relpct)/MEANv*100;                               % CV=std/mean
    fprintf('SET %d | %d組 相對預測誤差: min=%.3f mean=%.3f max=%.3f %% | CV=std/mean=%.1f%%\n', ...
            SET, Np, min(relpct),MEANv,max(relpct), CV);
    fprintf('=== SET %d 的 %d 組電流組合（paper 序 [I_P1 I_P2 I_P3 I_P4 I_P5 I_P6]） ===\n', SET, Np);
    for n=1:Np, fprintf('%3d: [%2d %2d %2d %2d %2d %2d]  rel=%.3f%%\n', n, pool(n,:), relpct(n)); end
end

% ---- 直方圖 ----
fig_dir=fullfile(CAL,'voltage_base','figures','eighteen_param'); if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
if SET==0
    % [-3,3] 實數密集取樣：Dmodel 風格（鋼藍 bar+白邊+右上方框）+ 紅虛線 mean、y=百分比、nb=180
    nb=180; edges=linspace(min(relpct),max(relpct),nb+1); ctr=(edges(1:end-1)+edges(2:end))/2;
    N=histcounts(relpct,edges); y=N/Np*100;                                     % y=百分比（current_sets/total×100%）
    fig=figure('Color','w','Position',[100 100 900 600]); ax=axes(fig); hold(ax,'on');
    bar(ax,ctr,y,1,'FaceColor',[0.20 0.40 0.70],'EdgeColor','w','LineWidth',0.2);
    hN =plot(ax,nan,nan,'LineStyle','none');                                    % dummy：total 佔 legend 一列
    hM =xline(ax,MEANv,'--','Color',[0.85 0.10 0.10],'LineWidth',2.2);          % mean 紅虛線（進 legend）
    hCV=plot(ax,nan,nan,'LineStyle','none');                                    % dummy：CV 佔 legend 一列
    hold(ax,'off'); xlim(ax,[0 max(relpct)*1.03]); ylim(ax,[0 max(y)*1.12]);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]); box(ax,'on'); grid(ax,'off');
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end)); yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel(ax,'$\left|\,S_i\,{}^{B}\hat{H}_{V}\,V - b_i\,I\,\right|\,/\,\left|\,b_i\,I\,\right|\ \ (\%)$','Interpreter','latex','FontSize',18);
    ylabel(ax,'current sets count / total current sets  (%)','FontWeight','bold');
    legend(ax,[hN hM hCV],{sprintf('total current sets = %d',Np),sprintf('mean = %.2f%%',MEANv),sprintf('CV = std/mean = %.0f%%',CV)}, ...
           'Location','northeast','FontSize',13,'FontWeight','bold','Box','on');   % 紅虛線標在 mean 列
    out=fullfile(fig_dir,'combo_region_hist_allreal.png'); exportgraphics(fig,out,'Resolution',300);
    fprintf('已輸出 %s\n', out);
else
    % 50/100 組（原樣式：深藍半透明 + legend mean/CV、y=count/N×100%）
    nb=18; if SET==12, nb=20; end; edges=linspace(min(relpct),max(relpct),nb+1); ctr=(edges(1:end-1)+edges(2:end))/2;
    N=histcounts(relpct,edges); y=N/Np*100;
    fig=figure('Color','w','Units','inches','Position',[1 1 8.8 5.2]); ax=axes(fig); hold(ax,'on');
    bar(ax,ctr,y,1,'FaceColor',[0.00 0.20 0.60],'EdgeColor','w','LineWidth',0.3,'FaceAlpha',0.55);
    hM=xline(ax,MEANv,'--','Color',[0.85 0.10 0.10],'LineWidth',2.2);              % mean 紅虛線（進 legend）
    hCV=plot(ax,nan,nan,'LineStyle','none');                                        % dummy handle：CV 佔 legend 一列
    hold(ax,'off'); xlim(ax,[0 max(relpct)*1.05]); ylim(ax,[0 max(y)*1.20]);
    set(ax,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickLength',[.018 .018]); box(ax,'on'); grid(ax,'off');
    xt=get(ax,'XTick'); set(ax,'XTick',xt(1:2:end)); yt=get(ax,'YTick'); set(ax,'YTick',yt(1:2:end));
    xlabel(ax,'$\left|\,S\,{}^{B}\hat{H}_{V}\,V - b\,I\,\right|\,/\,\left|\,b\,I\,\right|\ \ (\%)$','Interpreter','latex','FontSize',18); ylabel(ax,sprintf('count / %d  (%%)',Np),'FontWeight','bold');
    legend(ax,[hM hCV],{sprintf('mean = %.2f%%',MEANv),sprintf('CV = std/mean = %.0f%%',CV)}, ...
           'Location','northeast','FontSize',13,'FontWeight','bold','Box','on');
    fname='combo_region_hist.png'; if SET==2, fname='combo_region_hist_set2.png'; elseif SET==12, fname='combo_region_hist_all100.png'; end
    out=fullfile(fig_dir,fname); exportgraphics(fig,out,'Resolution',150);
    fprintf('已輸出 %s\n', out);
end
