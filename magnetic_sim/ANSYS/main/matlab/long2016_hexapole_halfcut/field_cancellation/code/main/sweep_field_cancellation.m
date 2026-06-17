%% sweep_field_cancellation.m
%  Sweep the 64 source/sink sign combinations of the 6 poles and find which one
%  cancels the magnetic field over the R=50um workspace SPHERE (the actual
%  manipulation volume about the WP), and which is strongest; plot the xy slice.
%
%  Model: each pole drives +/-1 (equal magnitude, sign = current direction =
%  source/sink). Field = linear superposition of the 6 single-coil raw FEM
%  solutions:  B(p) = sum_j s_j * B_FEM,j(p).  NO all-source flip (physical test),
%  NO interpolation (real FEM nodes only).
%
%  SCORE per combo = mean|B| over the R=50um sphere about the WP centre
%  (347 real nodes). mean (not min) so a single near-zero node can't masquerade
%  as cancellation; min/max also recorded. s and -s give identical |B| (B -> -B):
%  64 combos = 32 distinct |B| patterns.
%
%  Plotting still uses the 500um xy slice (visual context); the title score is the
%  R=50um-sphere mean|B|.
%
%  Reads ANSYS_data/.../coil1..6/standard ('all'); writes only figures (no results/).

clear; clc;
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\common');        % ansys_path
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');                  % mt_constants/import_ansys_data
addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\matlab\long2016_hexapole_halfcut\field_cancellation\code\plot');  % plot_field_xy

model  = 'long2016_hexapole_halfcut';
cnst   = mt_constants();
R_ws   = cnst.R_norm;           % workspace radius [m] = 500 um (xy-slice plot range)
R_avg  = 50e-6;                 % SCORE region: R=50 um sphere about WP centre
zslab  = 25e-6;                 % xy-slice half-thickness [m]
N_I    = 6;

%% ---- load 6 coils; build score-sphere (R=50um) + plot-slice (500um) node sets ----
fprintf('loading 6 coils (standard ''all'') ...\n');
d1   = import_ansys_data(ansys_path(model,'coil1','standard'),'all','coil1');
zwp  = d1.z - cnst.SPH_OFST;                                  % WP-frame z
mSph = (d1.x.^2 + d1.y.^2 + zwp.^2) <= R_avg^2;              % R=50um sphere (score)
mDsc = (abs(zwp) < zslab) & (d1.x.^2 + d1.y.^2 <= R_ws^2);   % 500um xy slice (plot)
getf = @(d,m) [d.bx(m), d.by(m), d.bz(m)];
Fs = cell(1,N_I); Fd = cell(1,N_I);
Fs{1}=getf(d1,mSph); Fd{1}=getf(d1,mDsc);
for j = 2:N_I
    cn = sprintf('coil%d',j);
    dj = import_ansys_data(ansys_path(model,cn,'standard'),'all',cn);
    Fs{j}=getf(dj,mSph); Fd{j}=getf(dj,mDsc);
end
BXs=cell2mat(cellfun(@(f)f(:,1),Fs,'uni',0)); BYs=cell2mat(cellfun(@(f)f(:,2),Fs,'uni',0)); BZs=cell2mat(cellfun(@(f)f(:,3),Fs,'uni',0));
BXd=cell2mat(cellfun(@(f)f(:,1),Fd,'uni',0)); BYd=cell2mat(cellfun(@(f)f(:,2),Fd,'uni',0)); BZd=cell2mat(cellfun(@(f)f(:,3),Fd,'uni',0));
Xd=d1.x(mDsc); Yd=d1.y(mDsc);
fprintf('nodes: R=50um sphere=%d (score), 500um xy slice=%d (plot)\n', nnz(mSph), nnz(mDsc));

%% ---- sweep 64 sign combos, SCORE = mean|B| over the R=50um sphere ----
S = zeros(64,6);
for k = 0:63, S(k+1,:) = 2*bitget(k,1:6) - 1; end           % 0/1 -> -1/+1, P1..P6
meanB = zeros(64,1); minB = zeros(64,1); maxB = zeros(64,1);
for k = 1:64
    s = S(k,:);
    bmag = sqrt((BXs*s.').^2 + (BYs*s.').^2 + (BZs*s.').^2);
    meanB(k)=mean(bmag); minB(k)=min(bmag); maxB(k)=max(bmag);
end

%% ---- rank by mean|B| over the R=50um sphere (print one rep per +/- pair: 32) ----
[~,ord] = sort(meanB,'ascend');
isrep = S(:,1)==1;
fprintf('\n rank | sign(P1 P2 P3 P4 P5 P6) | mean|B|   min|B|    max|B|  [mT]  (over R=50um sphere)\n');
fprintf('------+------------------------+-------------------------------------\n');
cnt = 0;
for ii = 1:64
    k = ord(ii); if ~isrep(k), continue; end
    cnt = cnt+1;
    fprintf(' %2d   | %+d %+d %+d %+d %+d %+d        | %7.3f  %7.3f  %7.3f\n', ...
            cnt, S(k,:), meanB(k)*1e3, minB(k)*1e3, maxB(k)*1e3);
end

%% ---- best cancel / strongest / all-source reference ----
[~,kmin] = min(meanB);  [~,kmax] = max(meanB);
s_as = [-1 +1 -1 +1 +1 -1];  kas = find(all(S==s_as,2),1);  % all-source (flip lower P1/P3/P6)
fprintf('\n[min cancel ] s=[%+d %+d %+d %+d %+d %+d]  mean|B|=%.4f  min=%.4f  max=%.4f mT\n', S(kmin,:), meanB(kmin)*1e3, minB(kmin)*1e3, maxB(kmin)*1e3);
fprintf('[max field  ] s=[%+d %+d %+d %+d %+d %+d]  mean|B|=%.4f  min=%.4f  max=%.4f mT\n', S(kmax,:), meanB(kmax)*1e3, minB(kmax)*1e3, maxB(kmax)*1e3);
fprintf('[all-source ] s=[%+d %+d %+d %+d %+d %+d]  mean|B|=%.4f  min=%.4f  max=%.4f mT\n', S(kas,:),  meanB(kas)*1e3,  minB(kas)*1e3,  maxB(kas)*1e3);

%% ---- plot xy slice (500um) for min cancel, max, all-source; score = R=50um-sphere mean|B| ----
plot_field_xy(S(kmin,:), 'mincancel', Xd, Yd, BXd, BYd, BZd, meanB(kmin));
plot_field_xy(S(kmax,:), 'max',       Xd, Yd, BXd, BYd, BZd, meanB(kmax));
plot_field_xy(S(kas,:),  'allsource', Xd, Yd, BXd, BYd, BZd, meanB(kas));
fprintf('done.\n');
