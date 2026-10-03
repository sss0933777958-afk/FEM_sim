%% xa_grad.m -- d(b.b)/dx_a on the x_a axis: one field model at four data ranges,
%%              plus the 18-parameter charge model from the flux calibration
%
%  METHOD selects the field model; everything else is identical between the two, so
%  the two runs are directly comparable.
%
%    METHOD = 'rbf'   Gaussian RBF, kernel (rho, lam) per range
%    METHOD = 'sph'   solid harmonics, degree L per range
%
%  FIVE CURVES, one quantity, one unit (mT^2/um), P1 excited at 1 A.
%
%  THE FOUR FIELD-MODEL CURVES.  Fitted on ALL nodes inside each radius
%  (1771 / 8225 / 22455 / 47693) with that range's winner from range_scan.m -- the
%  lowest test NMAE among the settings whose gradient passes the smoothness gate
%  (d3(b.b) keeps one sign on all three actuator axes over |s| <= R):
%
%      R     rbf (rho, lam)      sph (L, K)     sph test NMAE
%      150   460, 1e-7           7, 63          0.1107 %
%      250   280, 1e-5           6, 48          0.2339 %
%      350   500, 1e-3           6, 48          1.3765 %
%      450   610, 1e-2           2,  8         27.7294 %   <- gate-limited, see below
%
%  Those settings were SELECTED on the scan's 80/20 split and are re-fitted here on
%  every node, so the curves show what each data range gives at its best rather than
%  what it scores on held-out points.  They are therefore not guaranteed optimal for
%  the all-node fit -- the relative strength of lam (and the conditioning behind the
%  choice of L) shifts when the node count changes.
%
%  WATCH R = 450 FOR sph.  The smoothness gate gets very hard to satisfy at large R,
%  so the surviving degree collapses to L = 2 (8 coefficients) with a 27.7 % test
%  NMAE.  That curve measures the gate, not the harmonic basis.
%
%  THE CHARGE MODEL (18 parameters, flux calibration):
%      grad(b.b)_k = (gB^2 / l_hat) * (Kbar_I * I)' * L_k * (Kbar_I * I)
%  L comes from build_L, which works in pbar = P/l_hat, so the kernel is
%  dimensionless and gB is its gain [mT/A].  The same expression written with a
%  real-distance kernel carries l_hat to the FIFTH power; the two are identical
%  because gB(dimensionless) = gB(real)/l_hat^2.  Checked numerically: power 1 lands
%  on the RBF curves, power 5 is 5.8e11 too small.
%
%  This is the GRADIENT, not the force -- the force is 0.5*mgB*UF = 22.55 times it,
%  deliberately not applied so that all five curves share one axis.
%
%  Writes data/long2016_hexapole_halfcut/.mat/xa_grad_<METHOD>.mat  [ADDED 2026-09-16]

clear;  clc;

METHOD  = 'rbf';                               % 'rbf' | 'sph'
MODEL='long2016_hexapole_halfcut';  GEOM='tip40um';  VARIANT='maxwell';
RS      = [150 250 350 450 480];               % data ranges [um]
RHO     = [460 280 500 610 690];               % rbf: winners from range_scan.mat;
LAM     = [1e-7 1e-5 1e-3 1e-2 1e-2];          %   R=480 interpolated (450->610, 500->740)
LSPH    = [7 6 6 2 2];                         % sph: winners from range_scan.mat
REVAL   = 150;                                 % NMAE is always scored on |p| <= REVAL
CALNAME = 'calib_current_maxwell_R150_eighteen.mat';   % flux 18-parameter source
SPAN    = 150;   NS = 601;                     % plotted segment |s| <= SPAN
POLE    = 1;                                   % x_a <- P1 excited

METHOD = validatestring(METHOD, {'rbf','sph'}, mfilename, 'METHOD');
here=fileparts(mfilename('fullpath'));  FMX=fileparts(fileparts(here));
MAIN=fileparts(fileparts(fileparts(FMX)));  FLX=fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT=fullfile(FMX,'utils','data');

cfg = model_config(MODEL, GEOM);
s   = linspace(-SPAN, SPAN, NS).';
Pq  = [s zeros(NS,2)];                         % the x_a axis
I   = zeros(6,1);  I(POLE) = 1;                % 1 A on P1

% the nodes are needed for the fixed evaluation set (and, for sph, for the fit too)
raw = extract_maxwell_data(cfg, 'all', VARIANT);
ad  = build_actuator_data(raw, cfg);
Pa  = ad.Pa*1e6;  Ba = ad.Ba;  rr = vecnorm(Pa,2,2);
ie  = rr <= REVAL;   Pe = Pa(ie,:);   Be = Ba(ie,:,:);   Ne = nnz(ie);
nBe = sum(vecnorm(reshape(permute(Be,[1 3 2]),[],3),2,2));
fprintf(['evaluation set: %d nodes inside R <= %d um (in-sample for every range)' newline], Ne, REVAL);

% ---- the field model at each range, all nodes --------------------------------
CV = struct('name',{},'R',{},'par',{},'ntr',{},'g',{},'nmae',{});
for q = 1:numel(RS)
    R0 = RS(q);   t0 = tic;
    switch METHOD
        case 'rbf'
            [bq, gq, iq] = rbf_field(R0, struct('rho',RHO(q),'lam',LAM(q)));
            nn  = size(iq.P,1);
            par = sprintf('rho %g, lam %.0e', RHO(q), LAM(q));
            nm  = sprintf('RBF R = %d \\mum', R0);
        case 'sph'
            in  = rr <= R0;
            [bq, gq, iq] = sph_field(R0, struct('P',Pa(in,:),'B',Ba(in,:,:)), ...
                                    struct('L',LSPH(q),'Rn',R0,'quiet',true));
            nn  = nnz(in);
            par = sprintf('L %d, K %d', LSPH(q), iq.K);
            nm  = sprintf('Harmonics R = %d \\mum', R0);
    end
    G = gq(Pq, I, 'du');        % 'du' explicit: the two models default differently
    g = G(:,1);                 % the x_a component
    % NMAE of this range's model on the FIXED R <= REVAL node set, all six excitations
    Bp = zeros(Ne,3,6);
    for k = 1:6, Ik = zeros(6,1); Ik(k) = 1;  Bp(:,:,k) = bq(Pe, Ik);  end
    nmae = sum(vecnorm(reshape(permute(Bp-Be,[1 3 2]),[],3),2,2)) / nBe * 100;
    CV(end+1) = struct('name',nm, 'R',R0, 'par',par, 'ntr',nn, 'g',g, 'nmae',nmae); %#ok<SAGROW>
    fprintf(['%s R=%3d (%s, %6d nodes, %5.0f s): grad %.4f .. %.4f | NMAE(R<=%d) = %.4f %%' newline], ...
            upper(METHOD), R0, par, nn, toc(t0), min(g), max(g), REVAL, nmae);
    clear gq iq
end

% ---- the 18-parameter charge model ------------------------------------------
C   = load(fullfile(FLX,'data',MODEL,'.mat',CALNAME));
lum = C.l_hat*1e6;                             % the .mat stores metres
assert(C.USE_BIAS == 1, 'xa_grad:notEighteen', 'calibration is not the 18-parameter one');
L   = build_L(Pq, lum, C.e, cfg.Pc_base);
Gc  = base_model(L, C.KI_bar, C.Fmap(:,POLE), C.gI_hat^2 / lum);
gc  = squeeze(Gc(1,:,1)).';
fprintf(['charge model (l_hat %.2f um, gI_hat %.4f): %.4f .. %.4f mT^2/um' newline], ...
        lum, C.gI_hat, min(gc), max(gc));
Bp = zeros(Ne,3,6);
for k = 1:6
    Sk = build_S_(Pe, lum, C.e, cfg.Pc_base);      % dimensionless Coulomb kernel
    Ik = zeros(6,1); Ik(k) = 1;
    Bp(:,:,k) = reshape(Sk * (C.gI_hat * C.KI_bar * C.Fmap(:,k)), 3, Ne).';
end
nmae_c = sum(vecnorm(reshape(permute(Bp-Be,[1 3 2]),[],3),2,2)) / nBe * 100;
fprintf(['eighteen parameters: NMAE(R<=%d) = %.4f %%' newline], REVAL, nmae_c);
CV(end+1) = struct('name','Eighteen parameters', 'R',C.R_select*1e6, ...
                   'par','18 par.', 'ntr',NaN, 'g',gc, 'nmae',nmae_c);

save(fullfile(DAT, sprintf('xa_grad_%s.mat',METHOD)), ...
     'CV','s','SPAN','POLE','RS','RHO','LAM','LSPH','METHOD','CALNAME','lum','REVAL','Ne');
fprintf(['wrote xa_grad_%s.mat' newline], METHOD);

% ---- local helper -------------------------------------------------------------
function S = build_S_(P, l_hat, e17, Pc_base)
%BUILD_S_  the flux-side dimensionless Coulomb kernel, 3Np x 6 (mirrors fitting.m)
    E = zeros(3,6);
    E(:,1)=e17(1:3);  E(:,2)=e17(4:6);  E(:,3)=e17(7:9);
    E(:,4)=e17(10:12); E(:,5)=e17(13:15);
    E(1,6)=e17(16);   E(2,6)=e17(17);
    E(3,6)=e17(1)-e17(4)+e17(8)-e17(11)+e17(15);
    Pc = Pc_base + E;
    Np = size(P,1);   pbar = P / l_hat;   S = zeros(3*Np,6);
    for k = 1:6
        d  = pbar - Pc(:,k).';
        r3 = sum(d.^2,2).^1.5;
        S(:,k) = reshape((d ./ r3).', 3*Np, 1);
    end
end
