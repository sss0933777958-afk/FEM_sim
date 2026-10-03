%% rho_dominance.m -- is rho the variable that decides whether a smooth fit exists?
%
%  The question this answers: in the two-stage search (smoothness gate, then test
%  NMAE), can lambda be treated as THE search direction with rho pinned, or does rho
%  decide the outcome?
%
%  Method.  Two radii, chosen because the diagnostic showed their winners are set by
%  DIFFERENT mechanisms:
%      R = 150   the NMAE minimum sits inside a wide smooth window -- overfitting
%                decides the answer, the gate never binds
%      R = 450   NMAE falls monotonically past the gate -- the GATE decides
%  For each, sweep rho over four decades.  At every rho, walk lambda from 1e1 down to
%  1e-12 in tenth-decade steps (131 of them) and record
%      n_pass      how many of the 131 steps pass the d3 sign test
%      window      the lambda range that passes, and how many disconnected pieces
%      NMAE_smooth the best test NMAE among the passing steps      <- the answer
%      NMAE_any    the best test NMAE ignoring the gate            <- the ceiling
%      edf         effective dof  sum_i d_i/(d_i+lam)  at both of those points,
%                  the candidate one-dimensional complexity axis
%
%  If n_pass and NMAE_smooth swing wildly with rho, rho is the dominant variable and
%  cannot be pinned.  If they barely move, pinning rho and walking lambda is sound.
%
%  [ADDED 2026-09-16]

clear;  clc;

MODEL='long2016_hexapole_halfcut';  GEOM='tip40um';  VARIANT='maxwell';
RHOL = [40 80 160 300 400 566 600 800 1200 2000 4000 8000];
LAMS = 10.^(1:-0.1:-12);
RCASE= [150 450];
POLE=[1 3 6];  SG=[1 1 -1];  NQ=301;  NTRMAX=4000;  FTRAIN=0.8;  SEED=0;

here=fileparts(mfilename('fullpath'));  FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX)));  FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT=fullfile(FMX,'utils','data');

cfg=model_config(MODEL,GEOM);  raw=extract_maxwell_data(cfg,'all',VARIANT);
ad=build_actuator_data(raw,cfg);  Pa=ad.Pa*1e6;  Ba=ad.Ba;  rr=vecnorm(Pa,2,2);

OUT=struct('R',{},'rho',{},'npass',{},'nseg',{},'lo',{},'hi',{},'nm_s',{},'lam_s',{}, ...
           'edf_s',{},'nm_a',{},'lam_a',{},'edf_a',{});
for R0 = RCASE
    in=rr<=R0; P0=Pa(in,:); B0=Ba(in,:,:); N0=nnz(in);
    rng(SEED); ix=randperm(N0); n1=round(FTRAIN*N0); i1=ix(1:n1); i2=ix(n1+1:end);
    if n1>NTRMAX, i1=i1(1:NTRMAX); end
    Ptr=P0(i1,:); Btr=B0(i1,:,:); ntr=numel(i1);
    Pte=P0(i2,:); Bte=B0(i2,:,:); nte=numel(i2);
    nB0=sum(vecnorm(reshape(permute(Bte,[1 3 2]),[],3),2,2));
    sA=linspace(-R0,R0,NQ).';  AXQ=cell(1,3);
    for a=1:3, Pq=zeros(NQ,3); Pq(:,a)=sA; AXQ{a}=Pq; end
    d2c=sum(Ptr.^2,2).';  D2t=sum(Ptr.^2,2)+d2c-2*(Ptr*Ptr.');  D2t(D2t<0)=0;
    d2e=sum(Pte.^2,2)+d2c-2*(Pte*Ptr.');  Y=reshape(Btr,ntr,18);

    fprintf('\n================ R = %d um   (ntr %d, nte %d) ================\n', R0, ntr, nte);
    fprintf('%6s %7s %5s %10s %10s %11s %9s %11s %9s\n', ...
            'rho','npass','seg','lam lo','lam hi','NMAE smooth','edf','NMAE any','edf');
    for rho = RHOL
        Phi=exp(-D2t/rho^2);  [V,d]=eig((Phi+Phi.')/2,'vector');  VtY=V.'*Y;
        Kte=exp(-d2e/rho^2);
        PH=cell(1,3); UU=cell(1,3);
        for a=1:3
            d2q=sum(AXQ{a}.^2,2)+d2c-2*(AXQ{a}*Ptr.');
            PH{a}=exp(-d2q/rho^2);  UU{a}=AXQ{a}(:,a)-Ptr(:,a).';
        end
        nm=nan(size(LAMS));  ok=false(size(LAMS));  ed=nan(size(LAMS));
        for q=1:numel(LAMS)
            lam=LAMS(q);  if min(d)+lam<=0, continue, end
            W=V*(VtY./(d+lam));  ed(q)=sum(d./(d+lam));
            good=true;
            for a=1:3
                k=POLE(a); Wk=W(:,(k-1)*3+(1:3)); ph=PH{a}; U=UU{a};
                b0=ph*Wk;  b1=((-2/rho^2)*(U.*ph))*Wk;
                b2=((-2/rho^2)*ph+(4/rho^4)*(U.^2.*ph))*Wk;
                b3=((12/rho^4)*(U.*ph)-(8/rho^6)*(U.^3.*ph))*Wk;
                v=SG(a)*2*sum(3*b1.*b2+b0.*b3,2); g=sign(v); nz=find(g~=0);
                if sum(diff(g(nz))~=0)>0, good=false; break, end
            end
            ok(q)=good;
            Bp=reshape(Kte*W,nte,3,6);
            nm(q)=sum(vecnorm(reshape(permute(Bp-Bte,[1 3 2]),[],3),2,2))/nB0*100;
        end
        ip=find(ok);  nseg=0;  lo=NaN; hi=NaN; nms=NaN; lams=NaN; edfs=NaN;
        if ~isempty(ip)
            nseg=1+sum(diff(ip)>1);  lo=min(LAMS(ip));  hi=max(LAMS(ip));
            [nms,jb]=min(nm(ip));  lams=LAMS(ip(jb));  edfs=ed(ip(jb));
        end
        [nma,ja]=min(nm);  lama=LAMS(ja);  edfa=ed(ja);
        OUT(end+1)=struct('R',R0,'rho',rho,'npass',numel(ip),'nseg',nseg,'lo',lo,'hi',hi, ...
                          'nm_s',nms,'lam_s',lams,'edf_s',edfs,'nm_a',nma,'lam_a',lama,'edf_a',edfa); %#ok<SAGROW>
        fprintf('%6g %4d/%3d %5d %10.1e %10.1e %11.4f %9.1f %11.4f %9.1f\n', ...
                rho, numel(ip), numel(LAMS), nseg, lo, hi, nms, edfs, nma, edfa);
        clear Phi V d
    end
    s = OUT([OUT.R]==R0);
    v = [s.nm_s];  v = v(isfinite(v));
    fprintf('--> across rho: npass %d..%d of %d | best-smooth NMAE %.4f..%.4f %% (spread %.1fx)\n', ...
            min([s.npass]), max([s.npass]), numel(LAMS), min(v), max(v), max(v)/min(v));
end
save(fullfile(DAT,'rho_dominance.mat'),'OUT','RHOL','LAMS','RCASE');
fprintf('\nsaved rho_dominance.mat\n');
