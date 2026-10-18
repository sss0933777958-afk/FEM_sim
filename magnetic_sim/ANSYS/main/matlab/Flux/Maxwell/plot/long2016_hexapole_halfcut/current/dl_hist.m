%% dl_hist.m -- histogram of the noise-induced shift delta_l_hat (linearized, -J~'/J~'').
%  Reads utils/data/long2016_hexapole_halfcut/calib_noise[_bias].mat and utils/data/NTU_hexapole/calib_noise_<CASE>.mat
%  (written by utils/scripts/{long2016_hexapole_halfcut,NTU_hexapole}/calib_noise.m).
%  No fitting, no field reading: load -> plot (plot-scripts-pure).
%  Style: figure rules. Single series -> no legend; two series -> overlaid, blue / red, legend.
%  Horizontal axis symmetric about its center, three equally spaced ticks with the end gaps equal to the spacing,
%  endpoints numbered by text; vertical axis endpoints unnumbered.  Bin width = tick step / 10 (shared).
%  Params: CASE = 'long2016' | 'single_pole' | 'pair_pole' (NTU) | 'ntu_overlay' (single_pole + pair_pole);
%          BIAS (long2016 only: false = single, e = 0 -> dl_hist.png; true = eighteen, e_hat fixed -> dl_hist_bias.png);
%          PDF (single case only): true = refitted l_hat = l_hat_0 + dl_fit as percentage per bin, overlaid
%          with the closed-form Gaussian N(l_hat_0, SIG_L^2) scaled to percentage per bin (pdf * BINW * 100) -> lhat_pdf[_bias].png (long2016) / lhat_pdf.png (NTU);
%          the horizontal axis is centered at l_hat_0 rounded to 1 um; no legend (user request).
%  Output: figures/long2016_hexapole_halfcut/hexapole/dl_hist[_bias].png, figures/NTU_hexapole/<CASE>/dl_hist.png,
%          figures/NTU_hexapole/{single_pole,pair_pole}/dl_hist_sp_pp.png (overlay, same figure in both)
if ~exist('CASE','var'), CASE = 'long2016'; end
if ~exist('BIAS','var'), BIAS = false; end
if ~exist('PDF','var'), PDF = false; end
clc;  clearvars -except CASE BIAS PDF;  close all;
NTK  = 3;              % ticks per axis
here = fileparts(mfilename('fullpath'));
FLX  = fileparts(fileparts(fileparts(here)));                       % .../Flux/Maxwell
DAT  = fullfile(FLX, 'utils', 'data', 'NTU_hexapole');
LEG  = {};
switch CASE
    case 'long2016'
        TAG = '';  if BIAS, TAG = '_bias'; end
        DAT   = fullfile(FLX, 'utils', 'data', 'long2016_hexapole_halfcut');
        DFILE = {['calib_noise' TAG '.mat']};
        OUT   = fullfile(FLX, 'figures', 'long2016_hexapole_halfcut', 'hexapole');
        FNAME = ['dl_hist' TAG '.png'];
        if PDF, FNAME = ['lhat_pdf' TAG '.png']; end
    case 'ntu_overlay'
        DFILE = {'calib_noise_single_pole.mat', 'calib_noise_pair_pole.mat'};
        LEG   = {'Single pole', 'Pair pole'};
        OUT   = {fullfile(FLX, 'figures', 'NTU_hexapole', 'single_pole'), ...
                 fullfile(FLX, 'figures', 'NTU_hexapole', 'pair_pole')};
        FNAME = 'dl_hist_sp_pp.png';
    otherwise
        DFILE = {['calib_noise_' CASE '.mat']};
        OUT   = fullfile(FLX, 'figures', 'NTU_hexapole', CASE);
        FNAME = 'dl_hist.png';
        if PDF, FNAME = 'lhat_pdf.png'; end
end
PDF = PDF && numel(DFILE) == 1;

FS = 60;  FSLEG = 45;  FSLAB = 44;  LWBOX = 5.0;  CANV = 14.5;
LWCRV = 6.0;           % PDF curve
CB = [0.05 0.10 0.95];  CR = [0.85 0.10 0.10];
COL = {CB, CR};

NS = numel(DFILE);  DL = cell(1, NS);  XC = 0;
for k = 1:NS
    if PDF
        q = load(fullfile(DAT, DFILE{k}), 'dl_fit', 'SIG2', 'SIG_L', 'NMC', 'l0');
        XC = round(q.l0);                                          % axis center [um]
        DL{k} = q.l0 + q.dl_fit - XC;                              % refitted l_hat, relative to XC
        MU = q.l0 - XC;  SL = q.SIG_L;
        fprintf('%s: l_hat_0 %.4f um | sigma %.5g mT, NMC %d: refit mean %.4f um, std %.4g um | closed-form sigma_l_hat %.4g um\n', ...
                DFILE{k}, q.l0, sqrt(q.SIG2), q.NMC, mean(DL{k}) + XC, std(DL{k}), SL);
    else
        q = load(fullfile(DAT, DFILE{k}), 'dl_lin', 'SIG2', 'NMC', 'l0');
        DL{k} = q.dl_lin;
        fprintf('%s: l_hat_0 %.4f um | sigma %.5g mT, NMC %d: mean %.4g um, std %.4g um, var %.4g um^2\n', ...
                DFILE{k}, q.l0, sqrt(q.SIG2), q.NMC, mean(DL{k}), std(DL{k}), var(DL{k}));
    end
end
dall = vertcat(DL{:});

% ---- horizontal axis: XC + [-(n+1)/2 s, (n+1)/2 s], ticks XC - s, XC, XC + s; s a multiple of 0.1, a nice step if it wastes <= 20 %
hm = (NTK+1)/2;
nice = sort(reshape(([1 1.5 2 2.5 3 4 5].') .* 10.^(-1:3), [], 1));
if PDF, nice = nice(nice >= 1 & mod(nice, 1) == 0); end            % PDF: integer step -> integer tick labels
smx  = max(abs(dall))/hm;
sx   = min(nice(nice >= smx));
if sx > 1.2*smx && ~PDF, sx = ceil(smx/0.1 - 1e-9) * 0.1; end   % PDF: always the nice step (round tick labels)
if PDF && strcmp(CASE, 'long2016') && ~BIAS, sx = 3; end          % long2016 single: 865-877, ticks 868/871/874 (user choice)
xr = [-hm*sx, hm*sx];   xt = round(((1:NTK) - hm)*sx*10)/10;
BINW = sx / 10;                                                    % bin width [um]
if PDF, BINW = 0.25; end                                           % PDF: 0.25 um bins (peak ~10 %, round vertical ticks)
edg = (floor(xr(1)/BINW):ceil(xr(2)/BINW)) * BINW;
ctr = (edg(1:end-1) + edg(2:end))/2;
p = zeros(NS, numel(ctr));
for k = 1:NS, p(k,:) = histcounts(DL{k}, edg) / numel(DL{k}) * 100; end   % percentage per bin
pmax = max(p(:));
if PDF
    uu = linspace(xr(1), xr(2), 1001);
    pg = exp(-(uu - MU).^2 / (2*SL^2)) / (SL*sqrt(2*pi)) * BINW * 100;   % N(l_hat_0, SIG_L^2) -> % per bin
    pmax = max(pmax, max(pg));
end

% ---- vertical axis: a nice step if it wastes <= 20 %
smin = 1.02*pmax/(NTK+1);
nice = sort(reshape(([1 1.5 2 2.5 3 4 5].') .* 10.^(-3:3), [], 1));
sy   = min(nice(nice >= smin));
if sy > 1.2*smin && ~PDF                                          % PDF: keep the nice step (round tick labels)
    d1 = 10^floor(log10(smin));
    sy = ceil(smin/d1 - 1e-9)*d1;
    if sy > 1.2*smin, dg = d1/10;  sy = ceil(smin/dg - 1e-9)*dg; end
end
yr = [0, (NTK+1)*sy];   yt = (1:NTK)*sy;

fig = figure('Color','w','Units','inches','Position',[0.5 0.5 CANV CANV]);
ax  = axes(fig, 'Units','normalized', 'Position',[0.215 0.175 0.700 0.725]);
try, ax.Toolbar = []; end                                            %#ok<TRYNC>
try, ax.Interactions = []; end                                       %#ok<TRYNC>
hold(ax,'on');  box(ax,'on');
hb = gobjects(1, NS);
for k = 1:NS
    if NS == 1
        hb(k) = bar(ax, ctr, p(k,:), 1, 'FaceColor',COL{k}, 'EdgeColor','none');
    else                                                          % overlaid: translucent so both stay visible
        hb(k) = bar(ax, ctr, p(k,:), 1, 'FaceColor',COL{k}, 'FaceAlpha',0.60, 'EdgeColor','none');
    end
end
if PDF
    hb(2) = plot(ax, uu, pg, '-', 'Color',CR, 'LineWidth',LWCRV);
end
xtl = arrayfun(@(v) sprintf('%g', round((XC + v)*10)/10), xt, 'UniformOutput', false);
set(ax, 'XLim',xr, 'YLim',yr, 'XTick',xt, 'XTickLabel',xtl, 'YTick',yt, ...
        'FontSize',FS, 'FontWeight','bold', 'LineWidth',LWBOX, ...
        'TickDir','in', 'TickLength',[0.018 0.018], 'Layer','top');
if PDF
    xlabel(ax, '$\mathbf{\hat{\ell}\;(\mu m)}$',                      'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$',       'Interpreter','latex', 'FontSize',FSLAB);
else
    xlabel(ax, '$\mathbf{\delta\hat{\ell}\;(\mu m)}$',      'Interpreter','latex', 'FontSize',FSLAB);
    ylabel(ax, '$\mathbf{Percentage\;(\%)}$',       'Interpreter','latex', 'FontSize',FSLAB);
end

% rule 4: horizontal endpoints numbered by hand
yoff = yr(1) - 0.022*diff(yr);
for xv = xr
    text(ax, xv, yoff, sprintf('%g', round((XC + xv)*10)/10), 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontSize',FS, 'FontWeight','bold', 'Clipping','off');
end

if numel(LEG) > 1
    lg = legend(ax, hb, LEG, 'Interpreter','tex', 'FontSize',FSLEG, 'FontWeight','bold', ...
                'Location','northeast', 'NumColumns',1);
    lg.Box = 'on';  lg.EdgeColor = 'k';  lg.LineWidth = LWBOX;  lg.Color = 'w';
    lg.ItemTokenSize = [55 25];
end

set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 CANV CANV], 'PaperSize',[CANV CANV]);
OUT = cellstr(OUT);                                                % overlay goes to both case folders
for k = 1:numel(OUT)
    name = fullfile(OUT{k}, FNAME);
    print(fig, name, '-dpng', '-r200');
    fprintf('saved %s\n', name);
end
close(fig);
