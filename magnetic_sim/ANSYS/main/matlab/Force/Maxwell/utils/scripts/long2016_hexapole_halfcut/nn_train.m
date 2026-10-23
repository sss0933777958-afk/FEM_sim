%% nn_train.m -- train the NN inverse model (nn_net.m) on the labels inside the 90 um cube.
%  Samples: every direction j (ls_fit.m order, N_j = 2522) x every 30 um grid node of |x|, |y|, |z| <= A
%  (A = 45 um, all octants, N_p = 27)  ->  N_j N_p samples.
%    input   x = [r_hat(phi_j, theta_j); p / R]   (R = 30 um = outermost node: p / R in {-1, 0, 1}; actuator frame)
%    output  Ihat_opt                             (unit force 1 pN, pole order P1..P6)
%  LAB = 'grb'   : global-minimum labels of nn_label_sgn.mat
%        'chain' : position-continued labels of chain.mat
%  SGN = 'dir'   : re-align the label signs (whole 6-vectors, force unchanged): the centre labels across neighbouring
%                  directions (max-|cos| spanning tree from (0, 0)), then every position to its direction's centre
%        'none'  : labels as stored
%  Training: trainnet, Adam, mean-squared error, mini-batch MB (last partial batch used), shuffled every epoch.
%  RESUME = true: continue from the net saved in the output file (EP more epochs at LR; loss appended;
%  ep_run, LR_run record the epochs and learning rate of every run).
%  CONV = true (LOSS = 'label'): stop before EP once the epoch-mean training loss changes by less than TOLE
%  (relative, 1e-4 = 0.01 %) between consecutive epochs for KE = 10 epochs in a row.
%  Force check: e_ij = || fhat_NN - fhat_d || / || fhat_d || * 100 %,  fhat_d = r_hat,  fhat_NN_k = I_NN' Q_k I_NN
%  (Q_k as nn_label.m).
%  Output: utils/data/long2016_hexapole_halfcut/nn_train[_chain].mat
%    net             trained dlnetwork
%    loss            training loss per iteration
%    e_ij            N_j x N_p   [%]   relative force error of the network output
%    e_I             N_j x N_p   [A]   || I_NN - Ihat_opt ||
%    phi_j, theta_j, pbar, l_hat, R, A, MB, EP, LR, LAB, NH
%
%  LOSS = 'phys': no labels; the force model is in the loss,
%    L = || fhat(p_i, Ihat_NN) - r_hat ||^2 + LAM || Ihat_NN ||^2      (mean over the mini-batch)
%    positions: GRID = 'origin': multiples of DP inside [-A, A] on each axis (A = 45, DP = 30 um: -30, 0, 30 -> 27 points)
%               GRID = 'edge'  : -A:DP:A (DP = 10 um: 10^3 = 1000 points); actuator frame;
%    directions as ls_fit.m (phi, theta every 5 deg, N_j = 2522); input x = [r_hat; p / R] (R = 90 um);
%    custom loop (dlfeval, adamupdate), shuffled every epoch.  RESUME continues the saved run.
%    INIT = 'chain': start from the net trained on the chain labels (nn_train_chain.mat, same R) -> nn_phys_chain.mat.
%    LRD < 1: exponential learning-rate decay, LR of epoch e = LR * LRD^(e-1) -> nn_phys_exp.mat.
%    ep_conv: first epoch at which the epoch-mean loss changed by < TOLE (relative) for KE epochs in a row (NaN: never).
%    Output nn_phys.mat: net, loss (per epoch: total, force term, LAM term), e_ij (N_j x N_p [%]),
%    nI (N_j x N_p [A], || I_NN ||), phi_j, theta_j, P, R, A, DP, LAM, MB, EP, LR, NH, ttrain, ep_run, LR_run.
if ~exist('LOSS','var'), LOSS = 'label'; end
if ~exist('LAB','var'), LAB = 'grb'; end
if ~exist('R','var'), R = 30;  if strcmp(LOSS, 'phys'), R = 90; end, end
if ~exist('SGN','var'), SGN = 'dir'; end
if ~exist('A','var'), A = 45; end
if ~exist('DP','var'), DP = 30; end
if ~exist('LAM','var'), LAM = 0.01; end
if ~exist('GRID','var'), GRID = 'origin'; end
if ~exist('MB','var'), MB = 256; end
if ~exist('EP','var'), EP = 500; end
if ~exist('LR','var'), LR = 1e-3; end
if ~exist('NH','var'), NH = [64 64]; end
if ~exist('RESUME','var'), RESUME = false; end
if ~exist('CONV','var'), CONV = false; end
if ~exist('TOLE','var'), TOLE = 1e-4; end
if ~exist('KE','var'), KE = 10; end
if ~exist('INIT','var'), INIT = ''; end
if ~exist('LRD','var'), LRD = 1; end
clearvars -except LOSS LAB R A DP LAM GRID MB EP LR NH RESUME SGN CONV TOLE KE INIT LRD;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
addpath(fullfile(CALROOT, 'function'), here);                  % here: nn_net.m
MODEL = 'long2016_hexapole_halfcut';  GEOM = 'tip40um';
DAT = fullfile(CALROOT, 'utils', 'data', MODEL);
if strcmp(LOSS, 'phys')
    S = load(fullfile(DAT, 'nn_label_sgn.mat'), 'ph', 'th', 'MATNAME');
    % ---- samples: 1000 grid positions x 2522 directions ----
    g = DP * (ceil(-A/DP):floor(A/DP));                         % GRID = 'origin'
    if strcmp(GRID, 'edge'), g = -A:DP:A; end
    [GX, GY, GZ] = ndgrid(g, g, g);
    P = [GX(:), GY(:), GZ(:)];  N_p = size(P, 1);
    [PH, TH] = ndgrid(S.ph, S.th);
    U = TH < 360 & (abs(PH) < 90 | TH == 0);
    phi_j = PH(U);  theta_j = TH(U);  N_j = numel(phi_j);
    r_hat = [cosd(phi_j).*cosd(theta_j), cosd(phi_j).*sind(theta_j), sind(phi_j)];         % N_j x 3
    [J, Pi] = ndgrid(1:N_j, 1:N_p);  J = J(:);  Pi = Pi(:);       % sample s = (j, p)
    X = single([r_hat(J,:), P(Pi,:) / R]);                         % N_s x 6
    N_s = size(X, 1);
    % ---- force model: fhat_k = I' Q_k I  (Q_k as nn_label.m) ----
    mgB = 0.0451;  UF = 1e3;                                   % M g_B [A um^2/mT]; unit factor -> pN
    cal = load_flux_calib(MODEL, GEOM, S.MATNAME, 'current');
    gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;                 % F g_I [pN/A^2]
    FH  = sqrt(gF) * cal.Mbar_m;                               % F H_I [sqrt(pN)/A]
    L   = build_L(P, cal.l_m, cal.e_m, cal.Pc_base);           % N_p x 6 x 6 x 3
    Q   = zeros(6, 6, 3, N_p);
    for p = 1:N_p
        for k = 1:3
            Q_k = FH.' * squeeze(L(p,:,:,k)) * FH;  Q(:,:,k,p) = (Q_k + Q_k.')/2;
        end
    end
    Q = single(Q);  RH = single(r_hat);
    % ---- training ----
    OUT = 'nn_phys.mat';  if strcmp(INIT, 'chain'), OUT = 'nn_phys_chain.mat'; end
    if LRD ~= 1, OUT = 'nn_phys_exp.mat'; end
    lossE = zeros(0, 3);  ep_run = [];  LR_run = [];  t_prev = 0;
    if RESUME
        W = load(fullfile(DAT, OUT));
        assert(W.R == R && W.A == A && W.DP == DP && strcmp(W.GRID, GRID) && W.LAM == LAM && isequal(W.NH, NH), 'nn_train:resume', 'saved run differs');
        net = W.net;  lossE = W.loss;  ep_run = W.ep_run;  LR_run = W.LR_run;  t_prev = W.ttrain;
        avg = W.avg;  avgSq = W.avgSq;  it = W.it;
    else
        net = nn_net(NH);  avg = [];  avgSq = [];  it = 0;
        if strcmp(INIT, 'chain')                                % pretrained on the chain labels
            W0 = load(fullfile(DAT, 'nn_train_chain.mat'), 'net', 'R', 'A', 'NH');
            assert(W0.R == R && W0.A == A && isequal(W0.NH, NH), 'nn_train:init', 'nn_train_chain.mat has another R / A / NH');
            net = W0.net;
        end
    end
    rng(numel(ep_run));
    NB = ceil(N_s / MB);
    t0 = tic;
    for ep = 1:EP
        idx = randperm(N_s);  acc = zeros(1, 3);
        for b = 1:NB
            s  = idx((b-1)*MB+1 : min(b*MB, N_s));
            Xb = dlarray(X(s,:).', 'CB');
            it = it + 1;
            [Lb, gr, parts] = dlfeval(@phys_loss, net, Xb, Q(:,:,:,Pi(s)), RH(J(s),:).', LAM);
            [net, avg, avgSq] = adamupdate(net, gr, avg, avgSq, it, LR * LRD^(ep-1));
            acc = acc + numel(s) * double(extractdata([Lb, parts]));
        end
        lossE(end+1, :) = acc / N_s;                                %#ok<SAGROW>
        if ep == 1 || mod(ep, 10) == 0 || ep == EP
            fprintf('epoch %d/%d: loss %.4g (force %.4g, LAM |I|^2 %.4g), %.0f s\n', ep, EP, lossE(end,:), toc(t0));
        end
    end
    ttrain = t_prev + toc(t0);
    ep_run = [ep_run, EP];  LR_run = [LR_run, LR];
    % ---- force check ----
    e_ij = zeros(N_s, 1);  nI = zeros(N_s, 1);
    for c = 1:ceil(N_s / 1e5)
        s = (c-1)*1e5+1 : min(c*1e5, N_s);
        I = double(extractdata(predict(net, dlarray(X(s,:).', 'CB'))));   % 6 x n
        f = zeros(3, numel(s));
        for k = 1:3, f(k,:) = sum(I .* reshape(pagemtimes(double(squeeze(Q(:,:,k,Pi(s)))), reshape(I, 6, 1, [])), 6, []), 1); end
        e_ij(s) = vecnorm(f - r_hat(J(s),:).', 2, 1).' * 100;     % || fhat_d || = 1
        nI(s)   = vecnorm(I, 2, 1).';
    end
    e_ij = reshape(e_ij, N_j, N_p);  nI = reshape(nI, N_j, N_p);
    loss = lossE;
    rel = abs(diff(loss(:,1))) ./ loss(1:end-1,1);  run_ = 0;  ep_conv = NaN;
    for e = 1:numel(rel)
        if rel(e) < TOLE, run_ = run_ + 1; else, run_ = 0; end
        if run_ >= KE, ep_conv = e + 1;  break; end
    end
    fprintf('convergence (%d epochs in a row < %g %%): epoch %g\n', KE, 100*TOLE, ep_conv);
    save(fullfile(DAT, OUT), 'net', 'loss', 'e_ij', 'nI', 'phi_j', 'theta_j', 'P', 'R', 'A', 'DP', 'GRID', 'LAM', 'MB', 'EP', 'LR', 'NH', ...
         'ttrain', 'ep_run', 'LR_run', 'avg', 'avgSq', 'it', 'INIT', 'LRD', 'ep_conv', 'TOLE', 'KE');
    fprintf('LOSS phys: %d samples, MB %d, epochs %s, LR %s, %.0f s total | final loss %.4g\n', N_s, MB, mat2str(ep_run), ...
            mat2str(LR_run), ttrain, loss(end, 1));
    fprintf('e_ij [%%]: median %.3g, 95th pct %.3g, max %.3g | ||I|| [A]: median %.3g, max %.3g -> %s\n', median(e_ij(:)), ...
            prctile(e_ij(:), 95), max(e_ij(:)), median(nI(:)), max(nI(:)), fullfile(DAT, OUT));
    return
end
S   = load(fullfile(DAT, 'nn_label_sgn.mat'), 'P_act', 'ph', 'th', 'IOPT', 'l_hat', 'MATNAME');
l_hat = S.l_hat;

% ---- samples ----
n = find(all(abs(S.P_act) <= A + 1e-9, 2));  N_p = numel(n);
P = S.P_act(n,:);  pbar = P / l_hat;
[PH, TH] = ndgrid(S.ph, S.th);
U = TH < 360 & (abs(PH) < 90 | TH == 0);
phi_j = PH(U);  theta_j = TH(U);  N_j = numel(phi_j);
Ihat_opt = reshape(S.IOPT(:,:,:,n), 6, [], N_p);  Ihat_opt = Ihat_opt(:, U(:), :);        % 6 x N_j x N_p
if strcmp(LAB, 'chain')
    C = load(fullfile(DAT, 'chain.mat'), 'I_chain', 'pbar', 'phi_j', 'theta_j', 'A');
    [hit, loc] = ismembertol(pbar, C.pbar, 1e-12, 'ByRows', true);
    assert(all(hit) && C.A == A && isequal(C.phi_j, phi_j), 'nn_train:chain', 'chain.mat does not match');
    Ihat_opt = C.I_chain(:, :, loc);
end
neg_before = NaN;  neg_after = NaN;
if strcmp(SGN, 'dir')
    % (1) centre labels: one sign per direction so that neighbouring directions agree (max-|cos| spanning tree
    %     over the direction grid from (0, 0), as nn_sign.m); Ihat and -Ihat give the same force
    p0 = find(all(P == 0, 2));
    V  = Ihat_opt(:, :, p0);  V = V ./ vecnorm(V);             % 6 x N_j
    NPH = numel(S.ph);  NTH = numel(S.th);
    G = zeros(NPH, NTH);  G(U) = 1:N_j;
    for a = find(abs(S.ph(:)) == 90).', G(a, :) = G(a, 1); end  % pole: one direction
    G(:, end) = G(:, 1);                                       % theta = 360 is theta = 0
    s = [reshape(G(1:end-1, 1:end-1), [], 1); reshape(G(:, 1:end-1), [], 1)];   % along phi; along theta (wraps)
    t = [reshape(G(2:end,   1:end-1), [], 1); reshape(G(:, 2:end),   [], 1)];
    E = unique(sort([s t], 2), 'rows');  E = E(E(:,1) ~= E(:,2), :);
    c = sum(V(:, E(:,1)) .* V(:, E(:,2)), 1).';
    j0 = find(phi_j == 0 & theta_j == 0);
    [~, pred] = minspantree(graph(E(:,1), E(:,2), 1 - abs(c) + eps, N_j), 'Root', j0);
    ord = bfsearch(graph(find(pred > 0), pred(pred > 0), ones(nnz(pred > 0), 1), N_j), j0);
    sg = zeros(N_j, 1);  sg(j0) = 1;
    for k = ord(2:end).', sg(k) = sg(pred(k)) * sign(V(:, k).' * V(:, pred(k))); end
    Ihat_opt = Ihat_opt .* sg.';
    neg_before = nnz(c < 0);  neg_after = nnz(c .* sg(E(:,1)) .* sg(E(:,2)) < 0);
    % (2) every position of a direction takes the sign of that direction's centre label
    d0 = sum(Ihat_opt .* Ihat_opt(:, :, p0), 1);
    Ihat_opt = Ihat_opt .* (sign(d0) + (d0 == 0));
    fprintf('sign alignment: neighbouring directions with cos < 0 at the centre %d -> %d (of %d)\n', neg_before, neg_after, numel(c));
end
r_hat = [cosd(phi_j).*cosd(theta_j), cosd(phi_j).*sind(theta_j), sind(phi_j)];             % N_j x 3
[J, Pi] = ndgrid(1:N_j, 1:N_p);                                % sample s = (j, p)
X = [r_hat(J(:),:), P(Pi(:),:) / R];                           % N_s x 6
Y = reshape(permute(Ihat_opt, [2 3 1]), [], 6);                % N_s x 6, same (j, p) order
N_s = size(X, 1);

% ---- training ----
OUT = 'nn_train.mat';  if strcmp(LAB, 'chain'), OUT = 'nn_train_chain.mat'; end
rng(0);
loss0 = [];  ep_run = [];  LR_run = [];  t_prev = 0;
if RESUME                                                      % continue the saved run
    W = load(fullfile(DAT, OUT), 'net', 'loss', 'R', 'A', 'MB', 'NH', 'EP', 'LR', 'ttrain');
    assert(W.R == R && W.A == A && isequal(W.NH, NH), 'nn_train:resume', 'saved run has another R / A / NH');
    net = W.net;  loss0 = W.loss(:);  t_prev = W.ttrain;
    V = load(fullfile(DAT, OUT));
    if isfield(V, 'ep_run'), ep_run = V.ep_run;  LR_run = V.LR_run;  else, ep_run = W.EP;  LR_run = W.LR; end
    rng(numel(ep_run));                                        % new shuffle order for the new run
else
    net = nn_net(NH);
end
opts = trainingOptions('adam', 'MiniBatchSize',MB, 'MaxEpochs',EP, 'InitialLearnRate',LR, ...
                       'Shuffle','every-epoch', 'Verbose',true, 'VerboseFrequency',50*ceil(N_s/MB), ...
                       'Plots','none', 'ExecutionEnvironment','cpu');
if CONV, opts.OutputFcn = @(in) conv_stop(in, TOLE, KE); end
t0 = tic;
[net, info] = trainnet(X, Y, net, 'mse', opts);
ttrain = t_prev + toc(t0);
loss = [loss0; info.TrainingHistory.Loss(:)];
ep_done = max(info.TrainingHistory.Epoch);                     % < EP when CONV stopped the run
ep_run = [ep_run, ep_done];  LR_run = [LR_run, LR];
if CONV, fprintf('CONV: %d of %d epochs (stopped: %d)\n', ep_done, EP, ep_done < EP); end

% ---- force check (force model as nn_label.m) ----
I_NN = double(predict(net, X));                                % N_s x 6
mgB = 0.0451;  UF = 1e3;                                       % M g_B [A um^2/mT]; unit factor -> pN
cal = load_flux_calib(MODEL, GEOM, S.MATNAME, 'current');
gF  = UF * mgB / (2*cal.l_m) * cal.gB_m^2;                     % F g_I [pN/A^2]
FH  = sqrt(gF) * cal.Mbar_m;                                   % F H_I [sqrt(pN)/A]
L   = build_L(P, cal.l_m, cal.e_m, cal.Pc_base);               % N_p x 6 x 6 x 3
fhat_NN = zeros(N_s, 3);
for p = 1:N_p
    s = find(Pi(:) == p);
    for k = 1:3
        Q_k = FH.' * squeeze(L(p,:,:,k)) * FH;  Q_k = (Q_k + Q_k.')/2;
        fhat_NN(s,k) = sum((I_NN(s,:)*Q_k) .* I_NN(s,:), 2);
    end
end
fhat_d = r_hat(J(:),:);
e_ij = reshape(vecnorm(fhat_NN - fhat_d, 2, 2) ./ vecnorm(fhat_d, 2, 2) * 100, N_j, N_p);
e_I  = reshape(vecnorm(I_NN - Y, 2, 2), N_j, N_p);

save(fullfile(DAT, OUT), 'net', 'loss', 'e_ij', 'e_I', 'phi_j', 'theta_j', 'pbar', 'l_hat', 'R', 'A', 'MB', 'EP', 'LR', 'LAB', 'NH', ...
     'ttrain', 'ep_run', 'LR_run', 'SGN', 'neg_before', 'neg_after');
fprintf('LAB %s: %d samples, MB %d, epochs %s, LR %s, %.0f s total | final loss %.3g\n', LAB, N_s, MB, mat2str(ep_run), ...
        mat2str(LR_run), ttrain, loss(end));
fprintf('e_ij [%%]: median %.3g, 95th pct %.3g, max %.3g -> %s\n', median(e_ij(:)), prctile(e_ij(:), 95), max(e_ij(:)), fullfile(DAT, OUT));

function [Lb, gr, parts] = phys_loss(net, Xb, Qb, Rb, LAM)
% LOSS = 'phys':  mean_b || fhat(p_b, I_b) - r_hat_b ||^2 + LAM mean_b || I_b ||^2,  fhat_k = I' Q_k I.
I  = stripdims(forward(net, Xb));                         % 6 x B
B  = size(I, 2);
O  = reshape(I, 6, 1, B) .* reshape(I, 1, 6, B);          % I I^T, 6 x 6 x B
f  = [sum(O .* reshape(Qb(:,:,1,:), 6, 6, B), [1 2]); ...
      sum(O .* reshape(Qb(:,:,2,:), 6, 6, B), [1 2]); ...
      sum(O .* reshape(Qb(:,:,3,:), 6, 6, B), [1 2])];    % 3 x 1 x B
f  = reshape(f, 3, B);
lf = mean(sum((f - Rb).^2, 1));
lr = LAM * mean(sum(I.^2, 1));
Lb = lf + lr;
gr = dlgradient(Lb, net.Learnables);
parts = [lf, lr];
end

function stop = conv_stop(in, TOLE, KE)
% CONV: epoch-mean training loss; stop after KE consecutive epochs with |L_e - L_(e-1)| / L_(e-1) < TOLE
persistent ep acc n Lprev cnt
stop = false;
if strcmp(in.State, 'start') || isempty(ep), ep = 1;  acc = 0;  n = 0;  Lprev = NaN;  cnt = 0; end
if ~strcmp(in.State, 'iteration') || isempty(in.TrainingLoss), return, end
if in.Epoch ~= ep && n > 0                                     % epoch ep finished
    Le = acc / n;
    if abs(Le - Lprev) / Lprev < TOLE, cnt = cnt + 1;  else, cnt = 0; end
    Lprev = Le;  ep = in.Epoch;  acc = 0;  n = 0;
    stop = cnt >= KE;
end
acc = acc + double(in.TrainingLoss);  n = n + 1;
end
