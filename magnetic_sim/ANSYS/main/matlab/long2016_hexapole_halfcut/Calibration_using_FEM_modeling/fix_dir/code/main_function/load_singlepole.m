function [Pw, BmT] = load_singlepole(shape_dir, cnst, dataset)
%LOAD_SINGLEPOLE  Read one single-pole shape's FEM dataset (coord + B), WP frame.
%   [Pw, BmT] = LOAD_SINGLEPOLE(shape_dir, cnst, dataset)
%   dataset (optional, default 'tip'): 'tip' (R<2mm ball @ WP) or 'all' (full domain).
%   Reads <shape_dir>\coord_<ds>.dat (node x y z [m]) + bfield_<ds>.dat (node Bx By Bz |B| [T]),
%   matches by node id, shifts z into the WP frame (z - SPH_OFST). NO sign flip / no frame
%   rotation here (caller does actuator rotation + all-source flip).
%   Returns:
%     Pw   Nx3 [m]  node positions (WP frame)
%     BmT  Nx3 [mT] raw FEM B at those nodes (ANSYS Tesla ×1e3), no sign flip
    if nargin < 3 || isempty(dataset), dataset = 'tip'; end
    C = read_dat(fullfile(shape_dir, ['coord_'  dataset '.dat']), 4);   % node x y z
    B = read_dat(fullfile(shape_dir, ['bfield_' dataset '.dat']), 5);   % node Bx By Bz |B|
    [~, ia, ib] = intersect(C(:,1), B(:,1));                  % match by node id
    xyz  = C(ia, 2:4);
    bxyz = B(ib, 2:4);
    Pw  = [xyz(:,1), xyz(:,2), xyz(:,3) - cnst.SPH_OFST];     % WP frame [m]
    BmT = bxyz * 1e3;                                         % Tesla -> mT (raw)
end

function M = read_dat(f, ncol)
% Robust ANSYS PRINT reader: regexp-extract numbers per line (handles E-notation
% negatives concatenated without spaces, e.g. '...E-003-0.15...E-003'); keep only
% data rows with EXACTLY ncol numbers + integer positive node id.
    txt   = fileread(f);
    lines = regexp(txt, '\r?\n', 'split');
    pat   = '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?';
    rows  = cell(numel(lines),1);  k = 0;
    for i = 1:numel(lines)
        nums = regexp(lines{i}, pat, 'match');
        if numel(nums) == ncol
            v = str2double(nums);
            if v(1) > 0 && v(1) == round(v(1))       % node id = positive integer
                k = k+1;  rows{k} = v;
            end
        end
    end
    M = cell2mat(rows(1:k));
end
