%% plot_J_quiver.m — same style as B6x figure (alternating + random subsample)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
cnst = mt_constants();
K_I = eye(6) - ones(6)/6;
I_apdl = [1;-1;1;-1;-1;1];

% Re-run [J] ideal K_I to get positions
cube_half = 100e-6;
coil_data = struct();
for k=1:6
    dk=import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'),'wp',sprintf('coil%d',k));
    mask = abs(dk.x)<cube_half & abs(dk.y)<cube_half & abs(dk.z)<cube_half;
    N_max=5000;
    if sum(mask)>N_max, idx_all=find(mask); rng(0); idx_sub=idx_all(randperm(length(idx_all),N_max)); mask=false(size(mask)); mask(idx_sub)=true; end
    coil_data(k).px=dk.x(mask); coil_data(k).py=dk.y(mask); coil_data(k).pz=dk.z(mask);
    coil_data(k).bmag=sqrt(dk.bx(mask).^2+dk.by(mask).^2+dk.bz(mask).^2);
    coil_data(k).b_fem=[dk.bx(mask);dk.by(mask);dk.bz(mask)];
    coil_data(k).N=sum(mask);
    Iv=zeros(6,1); Iv(k)=1;
    coil_data(k).KI_w = K_I * Iv;
end
KI_data = load('../../data/KI_fit.mat');
opts=optimset('TolX',1e-10,'TolFun',1e-22,'MaxIter',50000,'MaxFunEvals',50000,'Display','off');
[xo,~]=fminsearch(@(x) cjoint(x,coil_data,1e-7), KI_data.pos(:), opts);
[~,C_k]=cjoint(xo,coil_data,1e-7);
pos=reshape(xo,3,6);

% (a) Alternating superposition on RANDOM subsampled FEM nodes in +-50um
d1=import_ansys_data(fullfile('..','..','results','coil1','filleted'),'wp','coil1');
cube_q=50e-6;
mask_q = abs(d1.x)<cube_q & abs(d1.y)<cube_q & abs(d1.z)<cube_q;
Nall = sum(mask_q);
px_all=d1.x(mask_q); py_all=d1.y(mask_q); pz_all=d1.z(mask_q);

% Random subsample 500 points
rng(42);
idx_500 = randperm(Nall, min(500, Nall));
px=px_all(idx_500); py=py_all(idx_500); pz=pz_all(idx_500);
Nq=length(px);

% FEM alternating at these random points
bxf=zeros(Nq,1);byf=bxf;bzf=bxf;
for k=1:6
    dk=import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'),'wp',sprintf('coil%d',k));
    Fx=scatteredInterpolant(dk.x,dk.y,dk.z,dk.bx,'natural');
    Fy=scatteredInterpolant(dk.x,dk.y,dk.z,dk.by,'natural');
    Fz=scatteredInterpolant(dk.x,dk.y,dk.z,dk.bz,'natural');
    bxf=bxf+I_apdl(k)*Fx(px,py,pz);
    byf=byf+I_apdl(k)*Fy(px,py,pz);
    bzf=bzf+I_apdl(k)*Fz(px,py,pz);
end

% Model alternating using [J] positions + per-coil C_k
bxm=zeros(Nq,1);bym=bxm;bzm=bxm;
for k=1:6
    Iv=zeros(6,1);Iv(k)=1;wk=K_I*Iv;
    bxk=zeros(Nq,1);byk=bxk;bzk=bxk;
    for i=1:6,dx=px-pos(1,i);dy=py-pos(2,i);dz=pz-pos(3,i);r3=(dx.^2+dy.^2+dz.^2).^(3/2);
    bxk=bxk+(-wk(i))*dx./r3;byk=byk+(-wk(i))*dy./r3;bzk=bzk+(-wk(i))*dz./r3;end
    bxm=bxm+I_apdl(k)*C_k(k)*1e-7*bxk;
    bym=bym+I_apdl(k)*C_k(k)*1e-7*byk;
    bzm=bzm+I_apdl(k)*C_k(k)*1e-7*bzk;
end
sc=500;

% (b) Error: 6 coils +-100um
Ns=coil_data(1).N;
pxs=coil_data(1).px;pys=coil_data(1).py;pzs=coil_data(1).pz;
all_rel=[];
for k=1:6
    bf=coil_data(k).b_fem; wk=coil_data(k).KI_w;
    bxk=zeros(Ns,1);byk=bxk;bzk=bxk;
    for i=1:6,dx=pxs-pos(1,i);dy=pys-pos(2,i);dz=pzs-pos(3,i);r3=(dx.^2+dy.^2+dz.^2).^(3/2);
    bxk=bxk+(-wk(i))*dx./r3;byk=byk+(-wk(i))*dy./r3;bzk=bzk+(-wk(i))*dz./r3;end
    bmod=[C_k(k)*1e-7*bxk;C_k(k)*1e-7*byk;C_k(k)*1e-7*bzk];
    res=bmod-bf;
    em=sqrt(res(1:Ns).^2+res(Ns+1:2*Ns).^2+res(2*Ns+1:3*Ns).^2);
    fm=sqrt(bf(1:Ns).^2+bf(Ns+1:2*Ns).^2+bf(2*Ns+1:3*Ns).^2);
    all_rel=[all_rel;100*em./fm];
end
me=mean(all_rel);Ntotal=length(all_rel);

fig=figure('Position',[50 50 1800 750]);

% (a) Both FEM and Model at SAME random points
subplot('Position',[0.05 0.08 0.42 0.82]);
h1=quiver3(px*1e6,py*1e6,pz*1e6,bxf*sc,byf*sc,bzf*sc,0,'b-','LineWidth',1.0);hold on;
h2=quiver3(px*1e6,py*1e6,pz*1e6,bxm*sc,bym*sc,bzm*sc,0,'r-','LineWidth',1.0);
xlabel('x (\mum)','FontSize',14);ylabel('y (\mum)','FontSize',14);zlabel('z (\mum)','FontSize',14);
xlim([-50 50]);ylim([-50 50]);zlim([-50 50]);
legend([h1 h2],'FEM','Model','Location','northeast','FontSize',13);
title('(a)','FontSize',20,'FontWeight','bold');
view(135,25);grid on;set(gca,'FontSize',13);box on;

% (b)
ax2=subplot('Position',[0.55 0.08 0.35 0.78]);
scatter(1:Ntotal,all_rel,3,'b','filled','MarkerFaceAlpha',0.4);hold on;
yline(me,'r--','LineWidth',2);
xlabel('Index of Points (6 coils)','FontSize',14);ylabel('Fitting Error (%)','FontSize',14);
title('(b)','FontSize',20,'FontWeight','bold');
grid on;set(gca,'FontSize',13);ylim([0 min(max(all_rel)*1.05,3)]);

annotation('textbox',[0.905 0.12 0.09 0.04],...
    'String',['mean = ' num2str(me,'%.2f') '%'],...
    'Color','r','FontSize',12,'FontWeight','bold',...
    'EdgeColor','none','BackgroundColor','none',...
    'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','none');
annotation('textbox',[0.55 0.88 0.42 0.06],...
    'String','e_i = |{\bfB}_{model} - {\bfB}_{FEM}|_i  /  |{\bfB}_{FEM}|_i',...
    'FontSize',12,'EdgeColor','k','BackgroundColor','w',...
    'HorizontalAlignment','center','FitBoxToText','on','Interpreter','tex','Margin',4);

set(fig,'PaperPositionMode','auto');
print(fig,fullfile('..','..','figures','analytic','fitting_J_idealKI_quiver_error.png'),'-dpng','-r200');
fprintf('Saved. mean=%.2f%%\n',me);

function [c,Ck]=cjoint(x,cd,km)
    p=reshape(x,3,6);c=0;Ck=zeros(6,1);
    for k=1:6, pw=[cd(k).px,cd(k).py,cd(k).pz];
    bu=evcf(p,cd(k).KI_w,km,pw);
    Ck(k)=(bu'*cd(k).b_fem)/(bu'*bu);
    r=Ck(k)*bu-cd(k).b_fem;c=c+sum(r.^2);end
end
function b=evcf(p,w,km,pw)
    N=size(pw,1);bx=zeros(N,1);by=bx;bz=bx;
    for i=1:6,dx=pw(:,1)-p(1,i);dy=pw(:,2)-p(2,i);dz=pw(:,3)-p(3,i);
    r3=(dx.^2+dy.^2+dz.^2).^(3/2);
    bx=bx+(-w(i))*dx./r3;by=by+(-w(i))*dy./r3;bz=bz+(-w(i))*dz./r3;end
    b=km*[bx;by;bz];
end
