% sp_cost_b.m -- profile cost J(b): for every b, a is re-fitted in closed form, b = 0..500 um.
% H(x) = a*I/(x - b)^2 (b = charge position, a = gain); data y = |b_x| [mT] on the 17 axis samples.
% a*(b) = sum(y.*u)/(I*sum(u.^2)), u = 1./(x - b).^2 ; J(b) = sum((y - a*(b) I u).^2) [mT^2].
if ~exist('CASE','var'), CASE = 'single_pole'; end   % 'single_pole' | 'pair_pole'
if ~exist('SAMP','var'), SAMP = 'tip'; end   % 'tip' | 'c500' (see sp_axis_calc.m)
clearvars -except CASE SAMP;  close all;
SFX = '';  if ~strcmp(SAMP, 'tip'), SFX = ['_' SAMP]; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
MAT = struct('single_pole','sp_axis.mat', 'pair_pole','pp_axis.mat');  MAT = MAT.(CASE);
MAT = fullfile(DAT, strrep(MAT, '.mat', [SFX '.mat']));
S  = load(MAT);
xi = S.x * 1e6;  y = abs(S.b(:,1));  I = S.I;
af = @(bb) sum(y ./ (xi-bb).^2) / (I * sum(1 ./ (xi-bb).^4));   % closed-form a for a given l_hat
Jf = @(bb) sum((y - af(bb) * I ./ (xi-bb).^2).^2);
BR = 500;  if strcmp(SAMP, 'c500'), BR = 2000; end   % sweep 0..BR um (c500: the pair minimum lies near 1560 um)
bs = linspace(0, BR, 50001).';
C  = arrayfun(Jf, bs);  A = arrayfun(af, bs);
[~, im] = min(C);
bmin = fminbnd(Jf, bs(max(im-1,1)), bs(min(im+1,end)), optimset('TolX',1e-9));
Cmin = Jf(bmin);  amin = af(bmin);
fprintf('profile: b_min = %.4f um (joint fit %.4f) | a = %.6g mT*um^2/A | J_min = %.6g\n', bmin, S.xc, amin, Cmin);
fprintf('J(0) = %.4g, J(500) = %.4g, max = %.4g, min/max ratio %.3g\n', C(1), C(end), max(C), Cmin/max(C));
save(strrep(MAT, '_axis', '_cost_b'), 'bs', 'C', 'A', 'bmin', 'Cmin', 'amin');
