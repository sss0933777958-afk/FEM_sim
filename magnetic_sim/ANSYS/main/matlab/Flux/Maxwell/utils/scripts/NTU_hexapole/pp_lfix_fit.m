% pp_lfix_fit.m -- NTU pair_pole (P1 excitation) or single_pole, centre samples (c500, 350..650 um from the tip):
% 1-D model H(x) = a*I/(x - l)^2 with l held fixed at L_LIST; a re-fitted in closed form for each l
%   a(l) = sum(y.*u) / (I*sum(u.^2)),  u = 1./(x - l).^2,  cost J(l) = sum((y - a I u).^2)  [mT^2].
% l = Inf: a/(x - l)^2 -> constant c, fitted as c = mean(y) (the l -> Inf limit of the same model).
% CASE = 'pair_pole' (default) | 'single_pole'.  L_LIST = fixed l values [um]; default per CASE
% (pair [500 1000 1560.5 3000], single [300 400 500 600]).
% SAMP = 'c500' (default) | 'cyl' (single_pole only: cylinder model, 17 points 100 .. 1700 um from the tip,
%        sp_axis_cyl.mat -> sp_lfix_fit_cyl.mat)
% Output: utils/data/NTU_hexapole/{pp,sp}_lfix_fit[_cyl].mat
%   L_LIST [um], a [mT*um^2/A] (NaN for Inf), c [mT] (value of the constant for Inf), J [mT^2], RMS [mT],
%   d, y = sample distance from the tip [um] and |b_x| [mT];  dc, Hc = model curves on a fine grid [um], [mT]
%   R = residual y - H at the samples, one column per l [mT]
if ~exist('CASE','var'), CASE = 'pair_pole'; end
if ~exist('L_LIST','var'), L_LIST = []; end
if ~exist('SAMP','var'), SAMP = 'c500'; end
clearvars -except CASE L_LIST SAMP;  close all;
SFX = '';  if ~strcmp(SAMP, 'c500'), SFX = ['_' SAMP]; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT  = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
PFX  = struct('pair_pole','pp', 'single_pole','sp');  PFX = PFX.(CASE);
S    = load(fullfile(DAT, [PFX '_axis_' SAMP '.mat']));
if isempty(L_LIST)                                    % fixed l_hat [um] (Inf also supported: constant fit)
    if strcmp(CASE, 'single_pole'), L_LIST = [300 400 500 600]; else, L_LIST = [500 1000 1560.5 3000]; end
end

x  = S.x(:)*1e6;  y = abs(S.b(:,1));  I = S.I;  d = -x;
xc = linspace(min(x), max(x), 301).';  dc = -xc;
n  = numel(L_LIST);  a = nan(1,n);  c = nan(1,n);  J = nan(1,n);  Hc = nan(numel(xc), n);  R = nan(numel(y), n);
for k = 1:n
    l = L_LIST(k);
    if isinf(l)
        c(k) = mean(y);  H = c(k)*ones(size(y));  Hc(:,k) = c(k);
    else
        u = 1./(x - l).^2;  a(k) = (y.'*u)/(I*(u.'*u));
        H = a(k)*I*u;  Hc(:,k) = a(k)*I./(xc - l).^2;
    end
    R(:,k) = y - H;  J(k) = sum(R(:,k).^2);
end
RMS = sqrt(J/numel(y));
save(fullfile(DAT, [PFX '_lfix_fit' SFX '.mat']), 'CASE', 'L_LIST', 'a', 'c', 'J', 'RMS', 'd', 'y', 'dc', 'Hc', 'R', 'I');
for k = 1:n
    fprintf('l = %7.1f um | a = %12.6g mT*um^2/A | J = %8.4f mT^2 | RMS = %.4f mT\n', L_LIST(k), a(k), J(k), RMS(k));
end
if any(isinf(L_LIST)), fprintf('l = Inf -> constant c = %.4f mT\n', c(isinf(L_LIST))); end
