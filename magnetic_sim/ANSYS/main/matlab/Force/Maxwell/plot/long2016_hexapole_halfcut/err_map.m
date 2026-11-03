% err_map.m -- relative force error of the fitted LS inverse model on one mid-plane, inside the fitting nodes (+-30 um),
% one direction (PH, TH) [deg] (continuous map, interpolated shading).  The origin is the black dot.
% Reads utils/data/long2016_hexapole_halfcut/err_map_p<PH>_t<TH>_<PL>.mat (utils/scripts/.../err_map.m).  Load -> plot only.
% Output: figures/long2016_hexapole_halfcut/hexapole/err_map_p<PH>_t<TH>_<PL>.png
% MODE = 'dir': 3D surface of e(theta, phi) at one position (err_map_dir.mat) -> err_map_dir.png
% MODE = 'pts': map e_f(theta, phi) every 2.5 deg at the centre and the +x_a, +y_a, +z_a edge points, NN and LS flux5,
%               calibration CLP (nn_train.m / ls_fit.m PX, set S25) -> emap_<nn|ls>_<CLP>_<o|x|y|z>.png
% MODE = 'loss': NN training loss per epoch (file NNF, default nn_phys.mat) up to epoch NE ([] = all, 'conv' = the
%               saved convergence epoch ep_conv: the horizontal axis then ends there) -> nn_loss.png
% MODE = 'R'  : at each of the 27 fitting nodes the largest e over the 2522 directions vs R = |p|, LS (ls_fit_chain.mat,
%               chain-label fit) and NN (file NNF, physics loss) overlaid -> err_R.png
% MODE = 'hist': histogram of the force error e_f of LS and NN overlaid (one shared bin width) on one data set,
%               SET = 'T12' (T1 + T2, grid directions) | 'T34' (T3 + T4, offset directions); LOSS = 'lab' | 'phys'
%               (test file of each in TF below) -> ef_hist_<LOSS>_<SET>[_<RNG>].png.  Each method in % of its own samples;
%               RNG = '' : frame ends at the larger of the two PCT-th percentiles (the share beyond is printed), NB bins
%                          per tick step
%                     'zoom': frame [0 ZM] exactly, bin ZM/100
%                     'full': frame covers the largest e_f of both, NB bins per tick step
if ~exist('MODE','var'), MODE = 'plane'; end
if ~exist('SET','var'), SET = 'T12'; end
if ~exist('LOSS','var'), LOSS = 'lab'; end
if ~exist('PCT','var'), PCT = 99; end
if ~exist('NB','var'), NB = 25; end
if ~exist('RNG','var'), RNG = ''; end
if ~exist('ZM','var'), ZM = 1; end                            % RNG = 'zoom': frame end [%]
if ~exist('GRP','var'), GRP = 'mth'; end                      % 'mth': LS vs NN on SET | 'set': LS only, T1 + T2 vs T3 + T4
                                                               %        -> ef_hist_<LOSS>_ls[_<RNG>].png
if ~exist('TL','var'), TL = []; end                           % GRP = 'set': only e_f > TL [%], each set in % of ALL its
                                                               %   samples; frame from the PCT-th pct of the kept values
                                                               %   -> ef_hist_<LOSS>_ls_gt<TL>.png
                                                               % GRP = 'inp': NN physics loss, input [r_hat; p/45] (seed 3)
                                                               %   vs input [p/45; phi; theta] (seed SA), all of T1-T4
                                                               %   -> ef_hist_phys_inp[_<RNG>].png
if ~exist('SA','var'), SA = 1; end                            % GRP = 'inp': seed of the angle-input net (best: 1)
                                                               % GRP = 'fwd': forward force-model NN (fwd_train.m), e_f on
                                                               %   its training samples, one series (no legend)
                                                               %   -> ef_hist_fwd<FWT>.png
if ~exist('FWT','var'), FWT = '_sph_r45_M10_e5000'; end       % GRP = 'fwd': tag of fwd_train<FWT>.mat
if ~exist('FWF','var'), FWF = 'fwd_train'; end                % GRP = 'fwd': 'fwd_test' -> fwd_test<FWT>.mat (test set,
                                                               %   e.g. FWT = '_sph_r45_mid_M10_e5000'; res stored there)
if ~exist('RES','var'), RES = false; end                      % GRP = 'fwd': true -> residual || fhat - f || [pN]
                                                               %   = e_f / 100 * || f || (f from fwd_data_sph_r45.mat)
                                                               %   instead of e_f -> res_hist_fwd<FWT>.png
if ~exist('VIEW','var'), VIEW = [-37.5 30]; end
if ~exist('PH','var'), PH = 0; end
if ~exist('TH','var'), TH = 0; end
if ~exist('PL','var'), PL = 'yz'; end
if ~exist('NE','var'), NE = []; end
if ~exist('NNF','var'), NNF = 'nn_phys.mat'; end
if ~exist('IVS','var'), IVS = 'train'; end                    % GRP = 'inv': inverse NN 6-630-6 on the sphere, residual
                                                               %   || fhat - r_hat || [10^-3 pN]; 'train' | 'test' (T2 points
                                                               %   x directions offset 2.5 deg) -> res_hist_inv_<IVS>.png
if ~exist('PDO','var'), PDO = 0; end                           % MODE = 'pts': direction offset [deg]; PDA = 5, PDO = 2.5: the offset test directions -> _off
if ~exist('CLIM','var'), CLIM = []; end                        % MODE = 'pts': [lo hi n] fixed colour range and n segments ([] = from the data)
if ~exist('X0M','var'), X0M = false; end                       % MODE = 'hist': horizontal axis starts at the data minimum
if ~exist('ETF','var'), ETF = ''; end                          % GRP = 'eta': '' training data | 'o' the centre -> ..._o_<ETS>
if ~exist('ETS','var'), ETS = 'lt'; end                        % GRP = 'eta': 'lt' eta < 1 (frame ends at 1) | 'gt' eta > 1 (starts at 1)
if ~exist('PDA','var'), PDA = 2.5; end                         % MODE = 'pts': direction step drawn [deg]; 5 = training directions -> _d5
if ~exist('PTM','var'), PTM = {'nn', 'ls'}; end                % MODE = 'pts': methods; 'm1' / 'm2' = Sec. 5.1 method 1 / 2 (ls_fit.m MTD, PX)
if ~exist('CLP','var'), CLP = 'single'; end                   % GRP = 'lsnn': LS vs 6-210-6 NN (flux5) on the test set, 300 points x
                                                               %   (5 deg + offset 2.5 deg directions), CLP 'single' | 'eighteen'
                                                               %   -> ef_hist_flux5_test_<CLP>[_<RNG>].png
clearvars -except PH TH PL MODE VIEW NE NNF SET LOSS PCT NB RNG ZM GRP TL SA FWT RES FWF IVS CLP PTM PDA ETS ETF X0M CLIM PDO;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(here));   % .../Force/Maxwell
FIGD = fullfile(CALROOT, 'figures', 'long2016_hexapole_halfcut', 'hexapole');

% ---- style (figure-style.md) ---------------------------------------------------------
FS = 60;  FSLAB = FS;  LWBOX = 5.0;  CANV = 14.5;  MSO = 22;

if strcmp(MODE, 'loss')                        % ---- NN (physics loss) training loss per epoch, log scale ----
    NN  = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', NNF));
    Le  = NN.loss(:,1);  ne = numel(Le);
    CB  = [0.05 0.10 0.95];
    nice = sort(reshape([1 1.5 2 2.5 3 4 5].' .* 10.^(0:4), [], 1));
    sx  = min(nice(nice >= ne/4));  XL = [0 4*sx];  XT = sx*(1:3);            % three ticks, end gap = step
    if strcmp(NE, 'conv')                                                      % axis ends at the convergence epoch
        ne = NN.ep_conv;  Le = Le(1:ne);  XL = [0 ne];  XT = round(ne*(1:3)/4);
    elseif ~isempty(NE)                                                        % NE: last epoch drawn
        ne = min(NE, ne);  Le = Le(1:ne);  sx = min(nice(nice >= ne/4));  XL = [0 4*sx];  XT = sx*(1:3);
    end
    YL  = [1e-4 1e0];  YT = [1e-3 1e-2 1e-1];                                  % decades, end gap = one decade
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig, 'Units','normalized', 'Position',[0.215 0.175 0.700 0.725]);  hold(ax, 'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    plot(ax, 1:ne, Le, '-', 'LineWidth',5, 'Color',CB);
    set(ax, 'YScale','log', 'XLim',XL, 'YLim',YL, 'XTick',XT, 'YTick',YT, 'FontSize',FS, 'FontWeight','bold', ...
            'LineWidth',LWBOX, 'Box','on', 'TickDir','in', 'Layer','top', 'YMinorTick','off');
    xlabel(ax, '$\mathbf{Epoch}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{Loss}$', 'Interpreter','latex', 'FontSize',FSLAB);
    for x = XL                                                                 % horizontal-axis end numbers
        text(ax, x, YL(1) / 10^(0.022*log10(YL(2)/YL(1))), num2str(x), 'HorizontalAlignment','center', ...
             'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold');
    end
    OUT = fullfile(FIGD, 'nn_loss.png');
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, OUT, '-dpng', '-r200');
    close(fig);
    fprintf('wrote %s | %d epochs, loss %.4g -> %.4g\n', OUT, ne, Le(1), Le(end));
    return
end

if strcmp(MODE, 'R')                           % ---- worst direction at every fitting node vs distance R ----
    DD  = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut');
    LS  = load(fullfile(DD, 'ls_fit_chain.mat'), 'e_ij', 'pbar', 'l_hat');      % LS, chain labels
    NN  = load(fullfile(DD, NNF), 'e_ij', 'P');                                % NN, physics loss
    RL  = vecnorm(LS.pbar, 2, 2) * LS.l_hat;   EL = max(LS.e_ij, [], 1).';      % R [um], e of the worst direction [%]
    RN  = vecnorm(NN.P, 2, 2);                 EN = max(NN.e_ij, [], 1).';
    CB  = [0.05 0.10 0.95];  CR = [0.85 0.10 0.10];  FSLEG = 45;  MSZ = 30;
    XT  = 15:15:45;  XL = [0 60];                                              % three ticks, end gap = step
    nice = sort(reshape([1 1.5 2 2.5 3 4 5].' .* 10.^(-2:3), [], 1));
    sy  = min(nice(nice >= 1.02*max([EL; EN])/4));  YT = sy*(1:3);  YL = [0 4*sy];
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig, 'Units','normalized', 'Position',[0.215 0.175 0.700 0.725]);  hold(ax, 'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    h(1) = plot(ax, RL, EL, 'o', 'MarkerSize',MSZ, 'LineWidth',4, 'Color',CB);
    h(2) = plot(ax, RN, EN, 's', 'MarkerSize',MSZ, 'LineWidth',4, 'Color',CR);
    set(ax, 'XLim',XL, 'YLim',YL, 'XTick',XT, 'YTick',YT, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
            'Box','on', 'TickDir','in', 'Layer','top');
    xlabel(ax, '$\mathbf{R\;(\mu m)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{e\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    for x = XL                                                                 % horizontal-axis end numbers
        text(ax, x, YL(1) - 0.022*diff(YL), num2str(x), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
             'FontSize',FS, 'FontWeight','bold');
    end
    lg = legend(ax, h, {'LS', 'NN'}, 'FontSize',FSLEG, 'FontWeight','bold', 'Location','northeast');
    lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];
    OUT = fullfile(FIGD, 'err_R.png');
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, OUT, '-dpng', '-r200');
    close(fig);
    fprintf('wrote %s | LS: max %.3g %%, mean %.3g %% | NN: max %.3g %%, mean %.3g %% (27 nodes, worst direction)\n', ...
            OUT, max(EL), mean(EL), max(EN), mean(EN));
    return
end

if strcmp(MODE, 'hist')                        % ---- e_f histogram, LS vs NN on one data set ----
    TF  = struct('lab', {{'ls_test_grb_d15.mat', 'nn_test_d15_e1000_s2.mat'}}, ...
                 'phys', {{'ls_test_phys_grb_d15.mat', 'nn_test_phys_exp_d15_s3_e1000.mat'}});
    TS  = {'T1', 'T2'};  if strcmp(SET, 'T34'), TS = {'T3', 'T4'}; end
    DD  = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut');
    if strcmp(GRP, 'fwd')                                                      % forward NN, one series
        L   = load(fullfile(DD, [FWF FWT '.mat']));
        e0  = L.e_f(:);  e1 = [];  LEG = {'Forward NN', ''};  TAG = FWT;  LOSS = 'fwd';
        if RES && isfield(L, 'res'), e0 = L.res(:) * 100;                     % [10^-2 pN]: integer ticks
        elseif RES                                                             % residual [pN], same sample order
            Dt = load(fullfile(DD, 'fwd_data_sph_r45.mat'), 'F');
            e0 = e0 / 100 .* vecnorm(Dt.F, 2, 2) * 100;
        end
        if strcmp(FWF, 'fwd_test'), TAG = ['_test' FWT]; end
        FXF = contains(FWT, '_flux_');                                         % flux net: unit force, no pN
        if FXF                                                                 % NMAE = sum ||fhat_NN - fhat|| / sum ||fhat||
            tk  = regexp(FWT, '_flux_([a-z]+)_x\d+(_sph_r45(_mid)?)_M', 'tokens', 'once');
            Dt  = load(fullfile(DD, ['fwd_data_flux_' tk{1} tk{2} '.mat']), 'F');
            fprintf('NMAE %s: %.4g %%\n', FWT, 100 * sum(L.res) / sum(vecnorm(Dt.F, 2, 2)));
            TAG = FWT;
        end
    else
        L   = load(fullfile(DD, TF.(LOSS){1}), TS{:});  N = load(fullfile(DD, TF.(LOSS){2}), TS{:});
        e0  = [L.(TS{1}).e(:); L.(TS{2}).e(:)];  e1 = [N.(TS{1}).e(:); N.(TS{2}).e(:)];   % LS, NN [%]
        LEG = {'LS', 'NN'};  TAG = SET;
    end
    if strcmp(GRP, 'set')                                                      % LS only, the two data sets
        L   = load(fullfile(DD, TF.(LOSS){1}), 'T1', 'T2', 'T3', 'T4');
        e0  = [L.T1.e(:); L.T2.e(:)];  e1 = [L.T3.e(:); L.T4.e(:)];
        LEG = {'T1 + T2', 'T3 + T4'};  TAG = 'ls';
    elseif strcmp(GRP, 'inp')                                                  % NN, r_hat input vs angle input, T1-T4
        TA  = {'T1', 'T2', 'T3', 'T4'};
        L   = load(fullfile(DD, 'nn_test_phys_exp_d15_s3_e1000.mat'), TA{:});
        N   = load(fullfile(DD, sprintf('nn_test_phys_exp_d15_s%d_ang_e1000.mat', SA)), TA{:});
        e0  = cell2mat(cellfun(@(t) L.(t).e(:), TA(:), 'UniformOutput', false));
        e1  = cell2mat(cellfun(@(t) N.(t).e(:), TA(:), 'UniformOutput', false));
        LEG = {'Unit vector', 'Angles'};  TAG = 'inp';  LOSS = 'phys';
    elseif strcmp(GRP, 'inv')                                                  % inverse NN 6-630-6 (sphere), one series
        if strcmp(IVS, 'test')                                                 % T2 points x directions offset 2.5 deg
            L = load(fullfile(DD, 'nn_test_phys_exp_sph_s0_h630_e1000.mat'), 'S2off');  e0 = L.S2off.e(:);
        else                                                                   % training points x 5 deg directions
            L = load(fullfile(DD, 'nn_phys_exp_sph_s0_h630_e1000.mat'), 'e_ij');  e0 = L.e_ij(:);
        end
        e0  = e0 * 10;  e1 = [];  LEG = {'Inverse NN', ''};  TAG = IVS;  LOSS = 'inv';   % e_f [%] of 1 pN -> residual [10^-3 pN]
    elseif strcmp(GRP, 'lsnn')                                                 % LS vs NN, flux5, test set
        L   = load(fullfile(DD, ['ls_test_flux5_' CLP '_sph_r45.mat']), 'S2', 'S2off');
        N   = load(fullfile(DD, ['nn_test_phys_exp_sph_s0_flux5_' CLP '_h210_e5000.mat']), 'S2', 'S2off');
        e0  = [L.S2.e(:); L.S2off.e(:)];  e1 = [N.S2.e(:); N.S2off.e(:)];
        LEG = {'LS', 'NN'};  TAG = ['test_' CLP];  LOSS = 'flux5';
    elseif strcmp(GRP, 'etaef')                                                % e_f of method 1 at the training data, eta < 1 or > 1
        if strcmp(ETF, 'o')                                                    % the centre: e_f of the grid directions of S25
            L  = load(fullfile(DD, ['ls_test_m1_' CLP '_sph_r45_pts.mat']), 'eta', 'S25');
            M_ = load(fullfile(DD, ['ls_m1_' CLP '_sph_r45.mat']), 'phi_j', 'theta_j');
            [~, k_] = ismember([M_.phi_j M_.theta_j], [L.S25.phi L.S25.theta], 'rows');
            e_ = L.S25.e(k_, 1);  L.eta = L.eta(:, 1);
        else
            L  = load(fullfile(DD, ['ls_m1_' CLP '_sph_r45.mat']), 'eta', 'e_q');
            e_ = L.e_q(1:size(L.eta, 1), :);
        end
        if strcmp(ETS, 'gt'), e0 = e_(L.eta > 1); else, e0 = e_(L.eta < 1); end
        e1 = [];  LEG = {'Method 1', ''};  TAG = ['m1_' CLP '_' ETS];  LOSS = 'etaef';
        if strcmp(ETF, 'o'), TAG = ['m1_' CLP '_o_' ETS]; end
        fprintf('eta %s 1: N %d, e_f mean %.4f, std %.4f, median %.4f, max %.3f, >1%%%% %.2f %%%%\n', ETS, numel(e0), mean(e0), std(e0), median(e0), max(e0), 100*mean(e0 > 1));
    elseif strcmp(GRP, 'eta')                                                  % eta = ||Phi_opt|| / ||Phi_LS||, method 1, centre,
        if strcmp(ETF, 'o')                                                    % the centre (ls_fit.m OPT), 2522 directions
            L = load(fullfile(DD, ['ls_test_m1_' CLP '_sph_r45_pts.mat']), 'eta');  L.eta = L.eta(:,1);
        else
            L = load(fullfile(DD, ['ls_m1_' CLP '_sph_r45.mat']), 'eta');      % 300 training points x 2522 directions
            L.eta = L.eta(:);
        end
        if strcmp(ETS, 'gt'), e0 = L.eta(L.eta > 1); else, e0 = L.eta(L.eta < 1); end
        e1 = [];  LEG = {'Method 1', ''};  TAG = ['m1_' CLP '_' ETS];  LOSS = 'eta';
        if strcmp(ETF, 'o'), TAG = ['m1_' CLP '_o_' ETS]; end
        fprintf('eta < 1: %d, eta > 1: %d, of %d samples\n', nnz(L.eta < 1), nnz(L.eta > 1), numel(L.eta));
    elseif strcmp(GRP, 'lsp')                                                  % LS flux5 (ls_fit.m PX): centre vs +x_a edge,
        L   = load(fullfile(DD, 'ls_test_flux5_single_sph_r45_pts.mat'));      % 5 deg + offset 2.5 deg directions
        e0  = [L.S2.e(:,1); L.S2off.e(:,1)];  e1 = [L.S2.e(:,2); L.S2off.e(:,2)];
        LEG = {'Center', '+x_a edge'};  TAG = 'pts_single';  LOSS = 'flux5';
    end
    if ~isempty(RNG), TAG = [TAG '_' RNG]; end
    n0  = numel(e0);  n1 = numel(e1);                                          % shares are of ALL samples of a set
    if strcmp(GRP, 'set') && ~isempty(TL)
        e0 = e0(e0 > TL);  e1 = e1(e1 > TL);  TAG = sprintf('%s_gt%g', TAG, TL);
    end
    CB  = [0.05 0.10 0.95];  CR = [0.85 0.10 0.10];  FSLEG = 45;  NTK = 3;
    nice = sort(reshape([1 1.5 2 2.5 3 4 5].' .* 10.^(-3:3), [], 1));
    sxc = nice(abs(nice*10 - round(nice*10)) < 1e-9);                         % tick steps: multiples of 0.1
    if strcmp(RNG, 'zoom')                                                     % frame exactly [0 ZM]: 4..6 segments
        for ns = 4:6, sx = ZM/ns;  if abs(sx*10 - round(sx*10)) < 1e-9, break; end, end
        XL  = [0 ZM];  XT = round((1:ns-1)*sx*10)/10;  BINW = ZM/100;
    else
        hi  = max(prctile(e0, PCT), prctile(e1, PCT));                        % neither method cut beyond 1 - PCT
        if strcmp(RNG, 'full'), hi = max([e0; e1]); end
        x0  = 0;  if strcmp(GRP, 'set') && ~isempty(TL), x0 = TL; end          % e_f > TL only: the axis starts at TL
        if X0M, x0 = floor(min([e0; e1])*10)/10; end                           % start at the data minimum (one decimal);
                                                                               % ticks then integer multiples of the step
        sx  = min(sxc(sxc >= (hi - x0)/(NTK+1)));
        XL  = x0 + [0 (NTK+1)*sx];  XT = round((x0 + (1:NTK)*sx)*10)/10;      % three ticks, end gap = step
        BINW = sx/NB;
        if X0M                                                                 % start at the minimum, ticks at multiples of the
            XT = (ceil(x0/sx):floor(XL(2)/sx)-1) * sx;                         % step from 0 (not equally spaced from the start),
            XL = [x0, max(XL(2) - x0 + x0, ceil(max([e0; e1])))];              % end at the step grid or the data maximum
            XL(2) = max(ceil(XT(end)/sx + 1) * sx, ceil(max([e0; e1])));
            if XL(2) - XT(end) > sx, XT = [XT, XT(end) + sx]; end
            XT = XT(XT > x0 & XT < XL(2) - 0.5*sx);
        end
        if strcmp(GRP, 'eta')                                                  % eta < 1: frame ends at 1; eta > 1: starts at 1
            if strcmp(ETS, 'gt'), sx = min(sxc(sxc >= (max(e0) - 1)/(NTK+1)));  XL = [1, 1 + (NTK+1)*sx];
            else,                 sx = min(sxc(sxc >= (1 - min(e0))/(NTK+1)));  XL = [1 - (NTK+1)*sx, 1]; end
            XT = round((XL(1) + (1:NTK)*sx)*10)/10;  BINW = sx/NB;
        end
    end
    edg = linspace(XL(1), XL(2), round(diff(XL)/BINW) + 1);                     % one bin width for both methods; the last
    ctr = (edg(1:end-1) + edg(2:end))/2;                                       % edge is the frame end
    p0  = histcounts(e0, edg) / n0 * 100;  p1 = histcounts(e1, edg) / max(n1, 1) * 100;
    smin = 1.02*max([p0 p1])/(NTK+1);  sy = min(nice(nice >= smin));
    if sy > 1.2*smin, dg = 10^floor(log10(smin))/10;  sy = ceil(smin/dg - 1e-9)*dg; end   % two significant digits
    YL  = [0 (NTK+1)*sy];  YT = (1:NTK)*sy;
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig, 'Units','normalized', 'Position',[0.215 0.175 0.700 0.725]);  hold(ax, 'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    C0  = CB;  if strcmp(GRP, 'fwd') && contains(FWT, '_mid'), C0 = CR; end      % forward NN: test set (T2) in red
    if strcmp(GRP, 'inv') && strcmp(IVS, 'test'), C0 = CR; end                  % inverse NN: test set in red
    h(1) = bar(ax, ctr, p0, 1, 'FaceColor',C0, 'FaceAlpha',0.60, 'EdgeColor','none');
    if n1 > 0, h(2) = bar(ax, ctr, p1, 1, 'FaceColor',CR, 'FaceAlpha',0.60, 'EdgeColor','none'); end
    set(ax, 'XLim',XL, 'YLim',YL, 'XTick',XT, 'YTick',YT, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
            'Box','on', 'TickDir','in', 'Layer','top');
    if strcmp(GRP, 'eta')
        xlabel(ax, '$\mathbf{\eta}$', 'Interpreter','latex', 'FontSize',FSLAB);
    elseif strcmp(GRP, 'lsnn')                                                  % test set: modeling error
        xlabel(ax, '$\mathbf{Modeling\;error\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    elseif strcmp(GRP, 'inv') && strcmp(IVS, 'test')                            % test set: modeling error
        xlabel(ax, '$\mathbf{Modeling\;error\;(10^{-3}\;pN)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    elseif strcmp(GRP, 'inv')
        xlabel(ax, '$\mathbf{Residual\;(10^{-3}\;pN)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    elseif strcmp(GRP, 'fwd') && RES && FXF && contains(FWT, '_mid')            % flux net, test set: unit force
        xlabel(ax, '$\mathbf{Modeling\;error\;(10^{-2})}$', 'Interpreter','latex', 'FontSize',FSLAB);
    elseif strcmp(GRP, 'fwd') && RES && FXF
        xlabel(ax, '$\mathbf{Residual\;(10^{-2})}$', 'Interpreter','latex', 'FontSize',FSLAB);
    elseif strcmp(GRP, 'fwd') && RES
        xlabel(ax, '$\mathbf{Residual\;(10^{-2}\;pN)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    else
        xlabel(ax, '$\mathbf{e_f\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    end
    ylabel(ax, '$\mathbf{Percentage\;count\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    for x = XL                                                                 % horizontal-axis end numbers
        text(ax, x, YL(1) - 0.022*diff(YL), sprintf('%g', x), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
             'FontSize',FS, 'FontWeight','bold');
    end
    if n1 > 0                                                                  % single series: no legend
        lg = legend(ax, h, LEG, 'FontSize',FSLEG, 'FontWeight','bold', 'Location','northeast');
        lg.LineWidth = LWBOX;  lg.ItemTokenSize = [55 25];
    end
    OUT = fullfile(FIGD, ['ef_hist_' LOSS '_' TAG '.png']);
    if strcmp(GRP, 'fwd'), OUT = fullfile(FIGD, ['ef_hist_fwd' TAG '.png']); end
    if strcmp(GRP, 'fwd') && RES, OUT = fullfile(FIGD, ['res_hist_fwd' TAG '.png']); end
    if strcmp(GRP, 'inv'), OUT = fullfile(FIGD, ['res_hist_inv_' IVS '.png']); end
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, OUT, '-dpng', '-r200');
    close(fig);
    fprintf('wrote %s | bin %g %% | frame [%g %g] %% | beyond (share of the set): %s %.2f %%, %s %.2f %%\n', OUT, BINW, ...
            XL(1), XL(2), LEG{1}, 100*nnz(e0 >= XL(2))/n0, LEG{2}, 100*nnz(e1 >= XL(2))/n1);
    fprintf('%s: N %d of %d, mean %.4g %%, std %.4g %% | %s: N %d of %d, mean %.4g %%, std %.4g %%\n', ...
            LEG{1}, numel(e0), n0, mean(e0), std(e0), LEG{2}, numel(e1), n1, mean(e1), std(e1));
    return
end

if strcmp(MODE, 'pts')                         % ---- e(theta, phi) at 4 positions, NN and LS (flux5), every 2.5 deg ----
    DD  = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut');
    F   = struct('nn', ['nn_test_phys_exp_sph_s0_flux5_' CLP '_h210_e5000_pts.mat'], 'ls', ['ls_test_flux5_' CLP '_sph_r45_pts.mat'], ...
                 'm1', ['ls_test_m1_' CLP '_sph_r45_pts.mat'], 'm2', ['ls_test_m2_' CLP '_sph_r45_pts.mat'], ...
                 'm2s1', ['ls_test_m2s1_' CLP '_sph_r45_pts.mat']);
    PN  = {'o', 'x', 'y', 'z'};  ph = -90+PDO:PDA:90-PDO;  th = PDO:PDA:360-PDA+PDO;  AXP = [0.20 0.28 0.56 0.56];
    for m = string(PTM)
        R = load(fullfile(DD, F.(m)), 'S25');  R = R.S25;
        k_ = mod(R.phi - PDO, PDA) == 0 & mod(R.theta - PDO, PDA) == 0;   % PDA = 5: the grid directions; PDO = 2.5: offset
        R.phi = R.phi(k_);  R.theta = R.theta(k_);  R.e = R.e(k_, :);
        [~, a] = ismember(R.phi, ph);  [~, b] = ismember(R.theta, th);
        for q = 1:size(R.e, 2)
            E = nan(numel(ph), numel(th));  E(sub2ind(size(E), a, b)) = R.e(:,q);
            if PDO == 0, E(1,:) = E(1,1);  E(end,:) = E(end,1); end            % phi = +-90: one direction
            lo = floor(min(E(:)));  hi = ceil(max(E(:)));  n = 4;
            if mod(hi - lo, 4) ~= 0 && mod(hi - lo, 3) == 0, n = 3; end
            if ~isempty(CLIM), lo = CLIM(1);  hi = CLIM(2);  n = CLIM(3); end      % fixed colour range, e.g. to match another map
            fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
            ax  = axes(fig, 'Units','normalized', 'Position', AXP);  hold(ax, 'on');
            try, ax.Toolbar = []; end
            try, ax.Interactions = []; end
            imagesc(ax, th, ph, E);                                            % raw values, one cell per direction
            colormap(ax, jet(256));  clim(ax, [lo hi]);
            cb = colorbar(ax);  cb.Limits = [lo hi];  cb.Ticks = linspace(lo, hi, n+1);  cb.LineWidth = LWBOX;
            cb.FontSize = FS;  cb.FontWeight = 'bold';
            ylabel(cb, '$\mathbf{e_f\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
            xlim(ax, [0 360]);  ylim(ax, [-90 90]);
            set(ax, 'XTick',90:90:270, 'YTick',-45:45:45, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'Box','on', ...
                    'Layer','top', 'YDir','normal');
            ax.Position = AXP;
            xlabel(ax, '$\mathbf{\theta\;(deg)}$', 'Interpreter','latex', 'FontSize',FSLAB);
            ylabel(ax, '$\mathbf{\varphi\;(deg)}$', 'Interpreter','latex', 'FontSize',FSLAB);
            for x = [0 360]                                                    % horizontal-axis end numbers
                text(ax, x, -90 - 0.012*180, num2str(x), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
                     'FontSize',FS, 'FontWeight','bold');
            end
            OUT = fullfile(FIGD, sprintf('emap_%s_%s_%s.png', m, CLP, PN{q}));
            if PDA ~= 2.5, OUT = strrep(OUT, '.png', sprintf('_d%g.png', PDA)); end
            if PDO ~= 0, OUT = strrep(OUT, '.png', '_off.png'); end
            set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
            print(fig, OUT, '-dpng', '-r200');
            close(fig);
            fprintf('%s %s %s: p = [%g %g %g] um, e mean %.3f, std %.3f, max %.3f %%, colour [%g %g]\n', m, CLP, PN{q}, R.P(q,:), mean(R.e(:,q)), std(R.e(:,q)), max(R.e(:,q)), lo, hi);
        end
    end
    return
end

if strcmp(MODE, 'dir')                         % ---- e over every direction at one position (3D) ----
    D   = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', 'err_map_dir.mat'));
    [TT, PP] = meshgrid(D.th, D.ph);
    SZ  = 1;  while 4*SZ < max(D.ED(:)), SZ = SZ * 2; end     % e ticks SZ, 2SZ, 3SZ; frame [0 4SZ]
    LIM = [0 360; -90 90; 0 4*SZ];
    fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
    ax  = axes(fig, 'Position', [0.17 0.16 0.72 0.74]);  hold(ax, 'on');
    try, ax.Toolbar = []; end
    try, ax.Interactions = []; end
    surf(ax, TT, PP, D.ED, D.ED, 'EdgeColor','none', 'FaceColor','interp');
    colormap(ax, jet(256));  clim(ax, LIM(3,:));
    grid(ax, 'off');  box(ax, 'off');  pbaspect(ax, [1 1 0.8]);
    xlim(ax, LIM(1,:));  ylim(ax, LIM(2,:));  zlim(ax, LIM(3,:));
    view(ax, VIEW(1), VIEW(2));
    set(ax, 'XTick',90:90:270, 'YTick',-45:45:45, 'ZTick',SZ*(1:3), 'FontSize',FS, 'FontWeight','bold', ...
            'LineWidth',LWBOX, 'TickLength',[0.025 0.025], 'XTickLabelRotation',0, 'YTickLabelRotation',0, ...
            'ZTickLabelRotation',0);
    draw_box_(ax, LIM, LWBOX);
    ax.Clipping = 'off';
    xlabel(ax, '$\mathbf{\theta\;(deg)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{\varphi\;(deg)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    zlabel(ax, '$\mathbf{e\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
    OUT = fullfile(FIGD, 'err_map_dir.png');
    set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
    print(fig, OUT, '-dpng', '-r200');
    close(fig);
    fprintf('wrote %s | p = [%g %g %g] um | e %.3g..%.3g %%\n', OUT, D.P0, min(D.ED(:)), max(D.ED(:)));
    return
end

D   = load(fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut', sprintf('err_map_p%g_t%g_%s.mat', PH, TH, PL)));
AXP = [0.20 0.28 0.56 0.56];                  % square plot box; room for the colour bar on the right
LIM = [-1 1] * D.A;  TK = (-1:1) * D.A/2;      % frame = outermost fitting nodes, three ticks, end gap = tick step

NAME = {'xy', 'xz', 'yz'};  AX = [1 2; 1 3; 2 3];  AXN = {'x_a', 'y_a', 'z_a'};
[lo, hi, sc] = cbar_ticks_([min(D.E(:)) max(D.E(:))]);
q = find(strcmp(NAME, D.PL));
fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig, 'Units','normalized', 'Position', AXP);  hold(ax, 'on');
try, ax.Toolbar = []; end
try, ax.Interactions = []; end
[UU, VV] = meshgrid(D.u, D.u);                  % interpolated shading: continuous map, no cell blocks
surface(ax, UU, VV, zeros(size(UU)) - 1, D.E, 'EdgeColor','none', 'FaceColor','interp');
plot(ax, 0, 0, 'o', 'MarkerSize',MSO, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
colormap(ax, jet(256));  clim(ax, [lo hi]);
cb = colorbar(ax);  cb.Limits = [lo hi];  cb.Ticks = lo:sc:hi;  cb.LineWidth = LWBOX;
cb.FontSize = FS;  cb.FontWeight = 'bold';
ylabel(cb, '$\mathbf{e\;(\%)}$', 'Interpreter','latex', 'FontSize',FSLAB);
xlim(ax, LIM);  ylim(ax, LIM);
set(ax, 'XTick',TK, 'YTick',TK, 'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, 'Box','on', 'Layer','top');
ax.Position = AXP;
xlabel(ax, ['$\mathbf{' AXN{AX(q,1)} '\;(\mu m)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
ylabel(ax, ['$\mathbf{' AXN{AX(q,2)} '\;(\mu m)}$'], 'Interpreter','latex', 'FontSize',FSLAB);
yl = ylim(ax);
for x = LIM                                     % horizontal-axis end numbers (not in XTick)
    text(ax, x, yl(1) - 0.012*diff(yl), num2str(x), 'HorizontalAlignment','center', 'VerticalAlignment','top', ...
         'FontSize',FS, 'FontWeight','bold');
end
OUT = fullfile(FIGD, sprintf('err_map_p%g_t%g_%s.png', PH, TH, PL));
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
print(fig, OUT, '-dpng', '-r200');
close(fig);
fprintf('wrote %s | plane %s | e %.3g..%.3g %%\n', OUT, D.PL, min(D.E(:)), max(D.E(:)));

% ---------------------------------------------------------------------------------------
function [lo, hi, s] = cbar_ticks_(r)
% colour-bar range and step: 3..4 segments, narrowest range; on a tie the larger step (fewer segments)
    best = Inf;
    for s0 = reshape([1 2 3 4 5].' * 10.^(-3:2), 1, [])
        s0 = round(s0, 10);
        l = round(floor(r(1)/s0)*s0, 10);  h = round(ceil(r(2)/s0)*s0, 10);  m = round((h - l)/s0);
        if m >= 3 && m <= 4 && h - l <= best + 1e-12, best = h - l;  lo = l;  hi = h;  s = s0; end
    end
end

function draw_box_(ax, lim, lw)
% box off + the 9 cube edges drawn by hand (the 3 edges meeting the corner nearest the camera are dropped).
    [Xc,Yc,Zc] = ndgrid(lim(1,:), lim(2,:), lim(3,:));  C = [Xc(:) Yc(:) Zc(:)];
    E = [];
    for i = 1:8
        for j = i+1:8
            if nnz(abs(C(i,:)-C(j,:)) > 0) == 1, E = [E; i j];  end %#ok<AGROW>
        end
    end
    cp = campos(ax);  [~, near] = min(sum((C - cp).^2, 2));
    E = E(~(E(:,1)==near | E(:,2)==near), :);
    for k = 1:size(E,1)
        p1 = C(E(k,1),:);  p2 = C(E(k,2),:);
        plot3(ax, [p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k-', 'LineWidth', lw);
    end
end
