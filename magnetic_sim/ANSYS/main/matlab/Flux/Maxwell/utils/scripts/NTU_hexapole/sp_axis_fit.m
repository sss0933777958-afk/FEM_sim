% sp_axis_fit.m -- 1-D charge fit on the 17 axis samples: H(x) = g*I/(x - x_c)^2, data = |b_x| [mT] as exported.
% cost(x_c, g) = sum_i (|b_x(x_i)| - H(x_i))^2 ; g is linear -> closed form for each x_c, 1-D search on x_c.
if ~exist('CASE','var'), CASE = 'single_pole'; end   % 'single_pole' | 'pair_pole'
if ~exist('SAMP','var'), SAMP = 'tip'; end   % 'tip' | 'c500' (see sp_axis_calc.m)
clearvars -except CASE SAMP;  close all;
SFX = '';  if ~strcmp(SAMP, 'tip'), SFX = ['_' SAMP]; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
MAT = struct('single_pole','sp_axis.mat', 'pair_pole','pp_axis.mat');  MAT = MAT.(CASE);
MAT = fullfile(DAT, strrep(MAT, '.mat', [SFX '.mat']));
S  = load(MAT);
xi = S.x * 1e6;                         % tip frame [um] (negative: in front of the tip)
y  = abs(S.b(:,1));  I = S.I;           % |b_x| [mT], drive current [A]
gof  = @(xc) sum(y ./ (xi-xc).^2) / (I * sum(1 ./ (xi-xc).^4));     % closed-form g for a given x_c
cost = @(xc) sum((y - gof(xc) * I ./ (xi-xc).^2).^2);
% scan both admissible regions (x_c cannot lie between samples)
scanA = linspace(max(xi)+0.01, 3000, 300001);      % right of all samples: x_c > -20 um
scanB = linspace(-3000, min(xi)-0.01, 300001);     % left of all samples (unphysical, checked only)
[cA, iA] = min(arrayfun(cost, scanA));  [cB, iB] = min(arrayfun(cost, scanB));
fprintf('scan  right: x_c = %.3f um  cost = %.6g | left: x_c = %.3f um  cost = %.6g\n', scanA(iA), cA, scanB(iB), cB);
dA = scanA(2)-scanA(1);
xc = fminbnd(cost, scanA(max(iA-1,1)), scanA(min(iA+1,end)), optimset('TolX',1e-9));
g  = gof(xc);  J = cost(xc);
r  = y - g * I ./ (xi-xc).^2;
fprintf('x_c = %.4f um   g = %.6g mT*um^2/A (= %.6g mT*mm^2/A)   cost = %.6g mT^2   RMS = %.4g mT\n', ...
        xc, g, g*1e-6, J, sqrt(J/numel(y)));
fprintf('%8s %10s %10s %10s\n', 'd[um]', 'data', 'H', 'err');
fprintf('%8.0f %10.4f %10.4f %10.4f\n', [-xi, y, g*I./(xi-xc).^2, r].');
dc = linspace(0, 400, 801).';  Hc = g * I ./ (-dc - xc).^2;          % curve vs distance from the tip
save(MAT, 'xc', 'g', 'J', 'r', 'dc', 'Hc', '-append');
