% pp_slice_y0.m -- NTU pair_pole, P1 excitation: field on the real FEM nodes of the plane y = 0 (magnetic circuit).
% Frame (same as pp_axis_calc.m / pp_sampling_plot.m): origin = right (P1) tip apex, pole axis +x (into the right
% plate), z from the plate mid-plane; left (P2) tip at x = -1000 um.  No interpolation: the arrows sit on the .fld
% grid nodes (20 um grid) with the exported b_x, b_z and |b| at those nodes.  So that the arrows do not look like a
% regular lattice, the nodes are picked in random order (rng 0) and a node is kept only if no kept node lies
% within DMIN of it (blue-noise subset of the real nodes).
% Output: utils/data/NTU_hexapole/pp_slice_y0.mat
%   X, Z [um], Bx, Bz, Bm [mT] (1 A excitation, field as exported);  XW, ZW window [um], DMIN, FLD, TIP
clearvars;  close all;
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Flux/Maxwell
DAT  = fullfile(CALROOT, 'utils', 'data', 'NTU_hexapole');
FLD  = 'D:\Maxwell_sim\NTU\export\pair_pole\P1.fld';
TIP  = [42.42464, 0, 0.125] * 1e3;                    % right (P1) tip apex, Maxwell global [um]
XW   = [-1100 100];  ZW = [-600 600];                  % window [um] (square, both tips inside)
DMIN = 34;                                            % minimum spacing between kept nodes [um]

fid = fopen(FLD);  fgetl(fid);  fgetl(fid);
C = textscan(fid, '%f %f %f %f %f %f', 'CollectOutput', true);  fclose(fid);
D = C{1};  clear C
x = round(D(:,1)*1e3 - TIP(1), 6);  y = round(D(:,2)*1e3 - TIP(2), 6);  z = round(D(:,3)*1e3 - TIP(3), 6);
h  = min(diff(unique(x)));
in = find(abs(y) < 1e-6 & x >= XW(1) & x <= XW(2) & z >= ZW(1) & z <= ZW(2));     % all y = 0 nodes in the window
rng(0);  ord = in(randperm(numel(in)));
sel = false(size(ord));  px = zeros(0,1);  pz = zeros(0,1);
for k = 1:numel(ord)
    i = ord(k);
    if isempty(px) || min((px - x(i)).^2 + (pz - z(i)).^2) >= DMIN^2
        sel(k) = true;  px(end+1,1) = x(i);  pz(end+1,1) = z(i);                     %#ok<AGROW>
    end
end
keep = ord(sel);
X  = x(keep);  Z = z(keep);
Bx = D(keep,4)*1e3;  Bz = D(keep,6)*1e3;  Bm = vecnorm(D(keep,4:6), 2, 2)*1e3;      % mT
save(fullfile(DAT, 'pp_slice_y0.mat'), 'X', 'Z', 'Bx', 'Bz', 'Bm', 'XW', 'ZW', 'DMIN', 'FLD', 'TIP');
fprintf('y = 0 nodes: %d in window (grid %.0f um), kept %d with spacing >= %d um | |b| %.2f .. %.1f mT\n', ...
        numel(in), h, numel(X), DMIN, min(Bm), max(Bm));
