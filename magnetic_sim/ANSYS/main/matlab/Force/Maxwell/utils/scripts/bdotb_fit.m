%% bdotb_fit.m -- b.b along the three actuator axes, RBF vs solid harmonics
%
%  Both models are fitted on THE SAME data: every FEM node with |p| <= R (1771 of
%  them at R = 150).  rbf_field loads them, and its info.P / info.B6 are handed
%  straight to sph_field, so there is no chance of the two seeing different sets.
%
%  Parameters are the winners of the R = 150 row of range_scan.m -- the lowest test
%  NMAE among those passing the smoothness gate:
%      RBF   rho = 460 um, lam = 1e-7   -> 0.1123 %
%      sph   L = 7  (K = 63)            -> 0.1107 %
%  (Those were selected on the 80 % training split; the curves here use all nodes,
%  which is the convention grad3.m already follows.)
%
%  Axes:  x_a <- P1 excited,  y_a <- P3,  z_a <- P6.  Each curve is drawn against the
%  axis as it is, NOT reoriented: P6 sits on -z_a, so the z_a curve rises towards
%  NEGATIVE s while the other two rise towards positive s.  That asymmetry is the
%  geometry and should be visible.  (grad3.m / sph_d3.m flip that axis because there
%  the plotted quantity is a vector COMPONENT and the flip keeps its sign meaningful;
%  b.b is a scalar, so there is nothing to keep consistent and flipping would only
%  hide which side the excited pole is on.)  FLIP_ZA is kept as a switch, set false.
%
%  Writes data/long2016_hexapole_halfcut/.mat/bdotb.mat for plot/.../bdotb.m
%  [ADDED 2026-09-16]

clear;  clc;

MODEL='long2016_hexapole_halfcut';  GEOM='tip40um';  VARIANT='maxwell';
R0      = 150;                 % data range [um]
RHO     = 460;   LAM = 1e-7;   % RBF winner at R = 150
LSPH    = 7;                   % harmonic degree winner at R = 150
NS      = 601;                 % samples along each axis
POLE    = [1 3 6];             % x_a <- P1, y_a <- P3, z_a <- P6
FLIP_ZA = false;               % draw each axis as it is; P6 lies on -z_a
AXN     = {'x_a','y_a','z_a'};

here=fileparts(mfilename('fullpath'));  FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX)));  FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT=fullfile(FMX,'utils','data');

% ---- both fits, one data set ------------------------------------------------
[bR,gR,iR] = rbf_field(R0, struct('rho',RHO,'lam',LAM));
fprintf(['RBF : rho %g um, lam %.0e, %d nodes' newline], RHO, LAM, size(iR.P,1));
[bS,gS,iS] = sph_field(R0, struct('P',iR.P,'B',iR.B6), struct('L',LSPH,'Rn',R0));
fprintf(['sph : L %d (K %d), same %d nodes, train rms %.4g mT' newline], LSPH, iS.K, size(iR.P,1), iS.rms);

% ---- b.b and its gradient along each axis ------------------------------------
%   dbb is the component of grad(b.b) ALONG that same axis, in mT^2/um -- the
%   quantity the force is proportional to.  It is NOT sign-flipped anywhere: on z_a
%   b.b rises towards -s, so its derivative there is negative, and that is the real
%   direction the force points (towards P6).
s = linspace(-R0, R0, NS).';
AX = struct('name',{},'s',{},'bbR',{},'bbS',{},'dbbR',{},'dbbS',{},'pole',{},'flipped',{});
for a = 1:3
    Pq = zeros(NS,3);  Pq(:,a) = s;
    I  = zeros(6,1);   I(POLE(a)) = 1;
    bbR  = sum(bR(Pq,I).^2, 2);                 % [mT^2]
    bbS  = sum(bS(Pq,I).^2, 2);
    % 'du' is passed EXPLICITLY to both: rbf_field defaults to it but sph_field
    % defaults to 'jac' (the Jacobian of B), so relying on the default silently
    % plots a different quantity for one of the two models.
    GR   = gR(Pq,I,'du');   dbbR = GR(:,a);     % [mT^2/um]
    GS   = gS(Pq,I,'du');   dbbS = GS(:,a);
    sd   = s;
    if a == 3 && FLIP_ZA
        sd = -s;  [sd,ord] = sort(sd);
        bbR = bbR(ord);  bbS = bbS(ord);  dbbR = -dbbR(ord);  dbbS = -dbbS(ord);
    end
    AX(a) = struct('name',AXN{a},'s',sd,'bbR',bbR,'bbS',bbS,'dbbR',dbbR,'dbbS',dbbS, ...
                   'pole',POLE(a),'flipped',(a==3 && FLIP_ZA)); %#ok<SAGROW>
    fprintf('%s (P%d): b.b  RBF %7.1f..%7.1f | sph %7.1f..%7.1f mT^2 (gap %.2f %%)\n', ...
            AXN{a}, POLE(a), min(bbR), max(bbR), min(bbS), max(bbS), ...
            max(abs(bbS-bbR)./max(bbR,eps))*100);
    fprintf('%s (P%d): dbb  RBF %+7.3f..%+7.3f | sph %+7.3f..%+7.3f mT^2/um (max abs gap %.4f)\n', ...
            AXN{a}, POLE(a), min(dbbR), max(dbbR), min(dbbS), max(dbbS), max(abs(dbbS-dbbR)));
end

save(fullfile(DAT,'bdotb.mat'), 'AX','R0','RHO','LAM','LSPH','NS','POLE','FLIP_ZA');
fprintf(['wrote bdotb.mat' newline]);
