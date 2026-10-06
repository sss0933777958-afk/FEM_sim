% sp_noise_mc.m -- Monte Carlo: add N(0,1) [mT] to each of the 17 |b_x| values, refit (a, b), NT times.
% Model H(x) = a/(x - b)^2; a closed form for each b, b by 1-D minimisation (coarse scan + fminbnd).
if ~exist('CASE','var'), CASE = 'single_pole'; end   % 'single_pole' | 'pair_pole' (P1 excitation)
PRE = struct('single_pole','sp', 'pair_pole','pp');  PRE = PRE.(CASE);
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
S  = load(fullfile(DAT, [PRE '_axis.mat']));
xi = S.x * 1e6;  y0 = abs(S.b(:,1));  I = S.I;        % tip frame [um], |b_x| [mT]
if ~exist('NT','var'), NT = 1000; end   % number of trials
SIG = 1;  rng(0);
bgrid = (max(xi)+0.5 : 0.5 : 3000).';
A = zeros(NT,1);  B = zeros(NT,1);  Jm = zeros(NT,1);
for t = 1:NT
    y  = y0 + SIG*randn(size(y0));
    af = @(b) sum(y./(xi-b).^2) / (I*sum(1./(xi-b).^4));
    Jf = @(b) sum((y - af(b)*I./(xi-b).^2).^2);
    [~, im] = min(arrayfun(Jf, bgrid));
    B(t) = fminbnd(Jf, bgrid(max(im-1,1)), bgrid(min(im+1,end)), optimset('TolX',1e-9));
    A(t) = af(B(t));  Jm(t) = Jf(B(t));
    assert(im > 1 && im < numel(bgrid), 'trial %d: minimum on the scan edge', t);
end
C = corrcoef(A, B);
fprintf('noise-free fit : b = %.4f um  a = %.6g mT*um^2/A\n', S.xc, S.g);
fprintf('MC (N = %d, sigma = %g mT):\n', NT, SIG);
fprintf('  b: mean %.4f  std %.4f um   (linearised std 2.304)\n', mean(B), std(B));
fprintf('  a: mean %.6g  std %.6g mT*um^2/A   (linearised std 4.181e4)\n', mean(A), std(A));
fprintf('  corr(a,b) = %.4f   (linearised 0.9738)\n', C(1,2));
fprintf('  b range %.3f .. %.3f um | a range %.6g .. %.6g\n', min(B), max(B), min(A), max(A));
save(fullfile(DAT, sprintf('%s_noise_mc_N%d.mat', PRE, NT)), 'A', 'B', 'Jm', 'NT', 'SIG', 'xi', 'y0');
