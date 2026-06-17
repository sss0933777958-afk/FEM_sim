%% fit_J_50um.m — [J] ideal K_I on +-50um (same range as Long)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
%  + RMSE figure + LaTeX table output
cnst = mt_constants();
K_I = eye(6) - ones(6)/6;
I_apdl = [1;-1;1;-1;-1;1];
cube_half = 50e-6;  % SAME as Long

% Load FEM data +-50um
coil_data = struct();
for k=1:6
    dk=import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'),'wp',sprintf('coil%d',k));
    mask = abs(dk.x)<cube_half & abs(dk.y)<cube_half & abs(dk.z)<cube_half;
    coil_data(k).px=dk.x(mask);coil_data(k).py=dk.y(mask);coil_data(k).pz=dk.z(mask);
    coil_data(k).bmag=sqrt(dk.bx(mask).^2+dk.by(mask).^2+dk.bz(mask).^2);
    coil_data(k).b_fem=[dk.bx(mask);dk.by(mask);dk.bz(mask)];
    coil_data(k).N=sum(mask);
    Iv=zeros(6,1);Iv(k)=1;coil_data(k).KI_w=K_I*Iv;
end
N = coil_data(1).N;
fprintf('+-50um: %d nodes per coil\n', N);

% Optimize
KI_data=load('../../data/KI_fit.mat');
opts=optimset('TolX',1e-10,'TolFun',1e-22,'MaxIter',50000,'MaxFunEvals',50000,'Display','off');
[xo,~]=fminsearch(@(x)cjoint(x,coil_data,1e-7),KI_data.pos(:),opts);
[~,C_k]=cjoint(xo,coil_data,1e-7);
pos=reshape(xo,3,6);

% Compute per-coil RMSE and relative error
pxs=coil_data(1).px;pys=coil_data(1).py;pzs=coil_data(1).pz;
all_sq_err=[];
all_rel=[];
for k=1:6
    bf=coil_data(k).b_fem;wk=coil_data(k).KI_w;
    bxk=zeros(N,1);byk=bxk;bzk=bxk;
    for i=1:6,dx=pxs-pos(1,i);dy=pys-pos(2,i);dz=pzs-pos(3,i);r3=(dx.^2+dy.^2+dz.^2).^(3/2);
    bxk=bxk+(-wk(i))*dx./r3;byk=byk+(-wk(i))*dy./r3;bzk=bzk+(-wk(i))*dz./r3;end
    bmod=[C_k(k)*1e-7*bxk;C_k(k)*1e-7*byk;C_k(k)*1e-7*bzk];
    res=bmod-bf;
    sq=res(1:N).^2+res(N+1:2*N).^2+res(2*N+1:3*N).^2;
    em=sqrt(sq); fm=sqrt(bf(1:N).^2+bf(N+1:2*N).^2+bf(2*N+1:3*N).^2);
    all_sq_err=[all_sq_err;sq];
    all_rel=[all_rel;100*em./fm];
end
rmse=sqrt(mean(all_sq_err))*1e3;
me_rel=mean(all_rel);
Ntotal=length(all_sq_err);
all_abs_err=sqrt(all_sq_err)*1e3;

fprintf('\n=== Results (+-50um, same as Long) ===\n');
fprintf('RMSE = %.4f mT\n', rmse);
fprintf('Mean relative error = %.2f%%\n', me_rel);
fprintf('Nodes: %d per coil, %d total\n', N, Ntotal);

% Print positions
fprintf('\n%-4s %-6s %10s %12s\n','Pole','Type','|l| [um]','Tip->Charge');
for i=1:6
    if cnst.pole_is_lower(i),ts='Lower';else,ts='Upper';end
    ell_i=norm(pos(:,i))*1e6;
    fprintf('P%d   %-6s %10.1f %12.1f\n',i,ts,ell_i,ell_i-500);
end
lower=logical(cnst.pole_is_lower);
fprintf('Lower avg: %.1f um\nUpper avg: %.1f um\n',...
    mean(vecnorm(pos(:,lower)))*1e6,mean(vecnorm(pos(:,~lower)))*1e6);

% === RMSE Figure ===
% (a) Alternating quiver
d1=import_ansys_data(fullfile('..','..','results','coil1','filleted'),'wp','coil1');
mask_q=abs(d1.x)<cube_half&abs(d1.y)<cube_half&abs(d1.z)<cube_half;
Nall=sum(mask_q);px_all=d1.x(mask_q);py_all=d1.y(mask_q);pz_all=d1.z(mask_q);
rng(42);idx_500=randperm(Nall,min(500,Nall));
px=px_all(idx_500);py=py_all(idx_500);pz=pz_all(idx_500);Nq=length(px);

bxf=zeros(Nq,1);byf=bxf;bzf=bxf;
bxm=zeros(Nq,1);bym=bxm;bzm=bxm;
for k=1:6
    dk=import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'),'wp',sprintf('coil%d',k));
    Fx=scatteredInterpolant(dk.x,dk.y,dk.z,dk.bx,'natural');
    Fy=scatteredInterpolant(dk.x,dk.y,dk.z,dk.by,'natural');
    Fz=scatteredInterpolant(dk.x,dk.y,dk.z,dk.bz,'natural');
    bxf=bxf+I_apdl(k)*Fx(px,py,pz);byf=byf+I_apdl(k)*Fy(px,py,pz);bzf=bzf+I_apdl(k)*Fz(px,py,pz);
    Iv=zeros(6,1);Iv(k)=1;wk=K_I*Iv;
    bxk=zeros(Nq,1);byk=bxk;bzk=bxk;
    for i=1:6,dx=px-pos(1,i);dy=py-pos(2,i);dz=pz-pos(3,i);r3=(dx.^2+dy.^2+dz.^2).^(3/2);
    bxk=bxk+(-wk(i))*dx./r3;byk=byk+(-wk(i))*dy./r3;bzk=bzk+(-wk(i))*dz./r3;end
    bxm=bxm+I_apdl(k)*C_k(k)*1e-7*bxk;bym=bym+I_apdl(k)*C_k(k)*1e-7*byk;bzm=bzm+I_apdl(k)*C_k(k)*1e-7*bzk;
end
sc=500;

fig=figure('Position',[50 50 1800 750]);
subplot('Position',[0.05 0.08 0.42 0.82]);
h1=quiver3(px*1e6,py*1e6,pz*1e6,bxf*sc,byf*sc,bzf*sc,0,'b-','LineWidth',1.0);hold on;
h2=quiver3(px*1e6,py*1e6,pz*1e6,bxm*sc,bym*sc,bzm*sc,0,'r-','LineWidth',1.0);
xlabel('x (\mum)','FontSize',14);ylabel('y (\mum)','FontSize',14);zlabel('z (\mum)','FontSize',14);
xlim([-50 50]);ylim([-50 50]);zlim([-50 50]);
legend([h1 h2],'FEM','Model','Location','northeast','FontSize',13);
title('(a)','FontSize',20,'FontWeight','bold');
view(135,25);grid on;set(gca,'FontSize',13);box on;

ax2=subplot('Position',[0.55 0.08 0.35 0.78]);
scatter(1:Ntotal,all_abs_err,3,'b','filled','MarkerFaceAlpha',0.4);hold on;
yline(rmse,'r--','LineWidth',2);
xlabel('Index of Points (6 coils)','FontSize',14);
ylabel('|{\bfB}_{model} - {\bfB}_{FEM}| [mT]','FontSize',14,'Interpreter','tex');
title('(b)','FontSize',20,'FontWeight','bold');
grid on;set(gca,'FontSize',13);ylim([0 min(max(all_abs_err)*1.05,0.2)]);

annotation('textbox',[0.905 0.12 0.09 0.04],...
    'String',['RMSE = ' num2str(rmse,'%.4f') ' mT'],...
    'Color','r','FontSize',11,'FontWeight','bold',...
    'EdgeColor','none','BackgroundColor','none',...
    'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','none');
annotation('textbox',[0.55 0.88 0.42 0.06],...
    'String','RMSE = \surd( mean( |{\bfB}_{model} - {\bfB}_{FEM}|^2 ) )',...
    'FontSize',12,'EdgeColor','k','BackgroundColor','w',...
    'HorizontalAlignment','center','FitBoxToText','on','Interpreter','tex','Margin',4);

set(fig,'PaperPositionMode','auto');
print(fig,fullfile('..','..','figures','analytic','fitting_J_idealKI_50um_RMSE.png'),'-dpng','-r200');
fprintf('\nFigure saved.\n');

% Save data
sv.pos=pos;sv.C_k=C_k;sv.rmse=rmse;sv.me_rel=me_rel;sv.cube_half=cube_half;sv.N=N;
save(fullfile('..','..','data','J_idealKI_50um.mat'),'-struct','sv');

function [c,Ck]=cjoint(x,cd,km)
    p=reshape(x,3,6);c=0;Ck=zeros(6,1);
    for k=1:6,pw=[cd(k).px,cd(k).py,cd(k).pz];
    bu=evcf(p,cd(k).KI_w,km,pw);Ck(k)=(bu'*cd(k).b_fem)/(bu'*bu);
    r=Ck(k)*bu-cd(k).b_fem;c=c+sum(r.^2);end
end
function b=evcf(p,w,km,pw)
    N=size(pw,1);bx=zeros(N,1);by=bx;bz=bx;
    for i=1:6,dx=pw(:,1)-p(1,i);dy=pw(:,2)-p(2,i);dz=pw(:,3)-p(3,i);
    r3=(dx.^2+dy.^2+dz.^2).^(3/2);
    bx=bx+(-w(i))*dx./r3;by=by+(-w(i))*dy./r3;bz=bz+(-w(i))*dz./r3;end
    b=km*[bx;by;bz];
end
