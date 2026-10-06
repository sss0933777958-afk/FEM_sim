% sp_axis_calc.m -- field on the 17 axis samples of the NTU single pole (tip frame), saved for plotting.
if ~exist('SAMP','var'), SAMP = 'tip'; end
% SAMP = 'tip'  : x = -20 .. -340 um, 20 um apart (17 points from the tip outward)
%        'c500' : 17 points on +-150 um around a centre 500 um in front of the tip (x = -650 .. -350 um)
SFX = '';  if ~strcmp(SAMP, 'tip'), SFX = ['_' SAMP]; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
addpath(fullfile(CALROOT, 'function', 'single_pole'));
cfg = sp_config('NTU_hexapole','single_pole');  L = sp_load(cfg);
DX = 20e-6;  NP = 17;  I = 1;                         % step [m], count, drive current [A] (50 A-turn / N_c 50)
switch SAMP
    case 'tip',  x = -(1:NP).' * DX;                                      % tip frame [m]
    case 'c500', x = -500e-6 + linspace(-150e-6, 150e-6, NP).';  DX = 300e-6/(NP-1);
    otherwise,   error('unknown SAMP ''%s''', SAMP);
end
P  = [x, zeros(NP,2)];
assert(all(L.air_stencil(P)), 'a sample touches iron');
b  = L.interp(P);                                     % mT
d  = -x;                                              % distance from the tip [m]
save(fullfile(DAT, ['sp_axis' SFX '.mat']), 'x', 'd', 'b', 'I', 'DX', 'NP');
fprintf('%8s %10s %10s %10s\n', 'd[um]', 'bx[mT]', 'by[mT]', 'bz[mT]');
fprintf('%8.0f %10.3f %10.4f %10.4f\n', [d*1e6, b].');
