%% Print per-pole positions from [J] ideal K_I
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'core'));
cnst = mt_constants();
K_I = eye(6) - ones(6)/6;
cube_half = 100e-6;
coil_data = struct();
for k=1:6
    dk=import_ansys_data(fullfile('..','..','results',sprintf('coil%d',k),'filleted'),'wp',sprintf('coil%d',k));
    mask = abs(dk.x)<cube_half & abs(dk.y)<cube_half & abs(dk.z)<cube_half;
    N_max=5000;
    if sum(mask)>N_max, idx_all=find(mask); rng(0); idx_sub=idx_all(randperm(length(idx_all),N_max)); mask=false(size(mask)); mask(idx_sub)=true; end
    coil_data(k).px=dk.x(mask);coil_data(k).py=dk.y(mask);coil_data(k).pz=dk.z(mask);
    coil_data(k).bmag=sqrt(dk.bx(mask).^2+dk.by(mask).^2+dk.bz(mask).^2);
    coil_data(k).b_fem=[dk.bx(mask);dk.by(mask);dk.bz(mask)];
    coil_data(k).N=sum(mask);
    Iv=zeros(6,1);Iv(k)=1;coil_data(k).KI_w=K_I*Iv;
end
KI_data=load('../../data/KI_fit.mat');
opts=optimset('TolX',1e-10,'TolFun',1e-22,'MaxIter',50000,'MaxFunEvals',50000,'Display','off');
[xo,~]=fminsearch(@(x) cjoint(x,coil_data,1e-7),KI_data.pos(:),opts);
pos=reshape(xo,3,6);

fprintf('%-4s %-6s %8s %8s %8s %10s %12s\n','Pole','Type','x[um]','y[um]','z[um]','|l|[um]','Tip->Charge');
fprintf('%s\n',repmat('-',1,60));
for i=1:6
    if cnst.pole_is_lower(i),ts='Lower';else,ts='Upper';end
    ell_i=norm(pos(:,i))*1e6;
    fprintf('P%d   %-6s %8.1f %8.1f %8.1f %10.1f %12.1f\n',i,ts,pos(1,i)*1e6,pos(2,i)*1e6,pos(3,i)*1e6,ell_i,ell_i-500);
end
lower=logical(cnst.pole_is_lower);
fprintf('%s\n',repmat('-',1,60));
fprintf('Lower avg: |l|=%.1f, Tip->Charge=%.1f\n',mean(vecnorm(pos(:,lower)))*1e6, mean(vecnorm(pos(:,lower)))*1e6-500);
fprintf('Upper avg: |l|=%.1f, Tip->Charge=%.1f\n',mean(vecnorm(pos(:,~lower)))*1e6, mean(vecnorm(pos(:,~lower)))*1e6-500);

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
