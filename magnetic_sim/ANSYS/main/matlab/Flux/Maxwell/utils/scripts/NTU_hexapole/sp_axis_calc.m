% sp_axis_calc.m -- field on the 17 axis samples of the NTU single pole (tip frame), saved for plotting.
if ~exist('SAMP','var'), SAMP = 'tip'; end
% SAMP = 'tip'  : x = -20 .. -340 um, 20 um apart (17 points from the tip outward)
%        'c500' : 17 points on +-150 um around a centre 500 um in front of the tip (x = -650 .. -350 um)
%        'cyl'  : cylinder model (Maxwell single_pole_cyl, export P1_cyc.fld): x = -100 .. -1700 um, 100 um apart
%        'cyl2' : cylinder model, two groups: 100 .. 400 um (100 um apart, 4 points) and 1000 .. 1960 um
%                 (80 um apart, 13 points) from the tip; DX = 0 (non-uniform)
SFX = '';  if ~strcmp(SAMP, 'tip'), SFX = ['_' SAMP]; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
addpath(fullfile(CALROOT, 'function', 'single_pole'));
cfg = sp_config('NTU_hexapole','single_pole');
if startsWith(SAMP, 'cyl'), cfg.fld_files = {'P1_cyc.fld'}; end   % same geometry and frame, air mesh cylinder model
L = sp_load(cfg);
DX = 20e-6;  NP = 17;  I = 1;                         % step [m], count, drive current [A] (50 A-turn / N_c 50)
switch SAMP
    case 'tip',  x = -(1:NP).' * DX;                                      % tip frame [m]
    case 'c500', x = -500e-6 + linspace(-150e-6, 150e-6, NP).';  DX = 300e-6/(NP-1);
    case 'cyl',  DX = 100e-6;  x = -(1:NP).' * DX;
    case 'cyl2', DX = 0;  x = -[100:100:400, 1000:80:1960].' * 1e-6;
    otherwise,   error('unknown SAMP ''%s''', SAMP);
end
assert(numel(x) == NP, 'expected %d samples', NP);
P  = [x, zeros(NP,2)];
assert(all(L.air_stencil(P)), 'a sample touches iron');
b  = L.interp(P);                                     % mT
d  = -x;                                              % distance from the tip [m]
save(fullfile(DAT, ['sp_axis' SFX '.mat']), 'x', 'd', 'b', 'I', 'DX', 'NP');
fprintf('%8s %10s %10s %10s\n', 'd[um]', 'bx[mT]', 'by[mT]', 'bz[mT]');
fprintf('%8.0f %10.3f %10.4f %10.4f\n', [d*1e6, b].');
