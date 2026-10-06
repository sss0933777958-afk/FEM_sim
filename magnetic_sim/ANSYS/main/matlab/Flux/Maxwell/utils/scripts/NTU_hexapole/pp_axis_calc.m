% pp_axis_calc.m -- NTU pair_pole: field of both excitations on the 17 axis samples, origin = right (P1) tip.
if ~exist('SAMP','var'), SAMP = 'tip'; end
% SAMP = 'tip'  : x = -20 .. -340 um, 20 um apart (17 points from the tip outward)
%        'c500' : 17 points on +-150 um around a centre 500 um in front of the tip (x = -650 .. -350 um)
SFX = '';  if ~strcmp(SAMP, 'tip'), SFX = ['_' SAMP]; end
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
FLD = {'D:\Maxwell_sim\NTU\export\pair_pole\P1.fld', 'D:\Maxwell_sim\NTU\export\pair_pole\P2.fld'};
TIP = [42.42464, 0, 0.125] * 1e-3;                    % right (P1) tip apex, Maxwell global [m] (OCC on pair_pole.STEP)
DX = 20e-6;  NP = 17;  I = 1;                         % step, count, drive current [A]
switch SAMP
    case 'tip',  x = -(1:NP).' * DX;                                      % tip frame [m]
    case 'c500', x = -500e-6 + linspace(-150e-6, 150e-6, NP).';  DX = 300e-6/(NP-1);
    otherwise,   error('unknown SAMP ''%s''', SAMP);
end
b  = zeros(NP, 3, 2);
for j = 1:2
    fid = fopen(FLD{j}); fgetl(fid); fgetl(fid); C = textscan(fid, '%f %f %f %f %f %f', 'CollectOutput', true); fclose(fid);
    D = C{1};  D(:,1:3) = D(:,1:3)*1e-3;
    xs = unique(D(:,1)) - TIP(1); ys = unique(D(:,2)) - TIP(2); zs = unique(D(:,3)) - TIP(3);
    nx = numel(xs); ny = numel(ys); nz = numel(zs);
    assert(nx*ny*nz == size(D,1) && D(2,3) ~= D(1,3));
    for c = 1:3
        G = permute(reshape(D(:,3+c)*1e3, nz, ny, nx), [3 2 1]);       % mT
        F = griddedInterpolant({xs, ys, zs}, G, 'linear', 'none');
        b(:,c,j) = F(x, zeros(NP,1), zeros(NP,1));
    end
end
d = -x;                                               % distance from the tip [m]
save(fullfile(DAT, ['pp_axis' SFX '.mat']), 'x', 'd', 'b', 'I', 'DX', 'NP', 'TIP', 'FLD');
fprintf('%6s | %9s %8s %8s | %9s %8s %8s\n', 'd[um]', 'P1: bx', 'by', 'bz', 'P2: bx', 'by', 'bz');
fprintf('%6.0f | %9.3f %8.3f %8.3f | %9.3f %8.3f %8.3f\n', [-x*1e6, b(:,:,1), b(:,:,2)].');
