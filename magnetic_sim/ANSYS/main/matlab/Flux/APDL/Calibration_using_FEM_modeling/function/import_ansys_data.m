function data = import_ansys_data(results_dir, dataset, coil_name)
% IMPORT_ANSYS_DATA  Read ANSYS POST1 exported coordinate + B-field data
%   data = import_ansys_data('results/coil1')                    % reads 'all' dataset
%   data = import_ansys_data('results/coil1', 'wp')              % reads working-point dataset
%   data = import_ansys_data('results/coil2', 'all', 'coil2')   % specify coil name
%   Returns struct with fields: node_id, x, y, z, bx, by, bz, bsum (SI units)

    if nargin < 2, dataset = 'all'; end
    if nargin < 3, coil_name = 'coil1'; end

    coord_file  = fullfile(results_dir, [coil_name '_coord_'  dataset '.dat']);
    bfield_file = fullfile(results_dir, [coil_name '_bfield_' dataset '.dat']);

    % --- Read coordinate file (space-separated, readmatrix handles it) ---
    %   [FIX] On large files readmatrix can mis-detect extra columns (banner line)
    %   -> trailing NaN cols. Check NaN only in the data cols (node,x,y,z) so valid
    %   rows are not wrongly dropped (was: coord -> 0 rows for >1M-node files).
    coord = readmatrix(coord_file, 'FileType', 'text');
    nc = min(4, size(coord,2));
    coord = coord(~any(isnan(coord(:,1:nc)), 2), :);

    % --- Read B-field file (fixed-width format, negatives may concatenate) ---
    % ANSYS PRNSOL uses fixed-width columns; consecutive negative values
    % like "-0.123E-002-0.456E-003" have no space separator.
    % Fix: insert space before any '-' or '+' that immediately follows a digit.
    txt = fileread(bfield_file);
    txt = regexprep(txt, '(\d)([-+])', '$1 $2');
    tmpfile = [tempname '.dat'];
    fid = fopen(tmpfile, 'w'); fprintf(fid, '%s', txt); fclose(fid);
    bfield = readmatrix(tmpfile, 'FileType', 'text');
    delete(tmpfile);
    nb = min(5, size(bfield,2));
    bfield = bfield(~any(isnan(bfield(:,1:nb)), 2), :);

    % --- Merge on node ID (PRNSOL skips SOURC36 coil nodes) ---
    [node_ids, ic, ib] = intersect(coord(:,1), bfield(:,1));
    fprintf('  Matched %d nodes (coord: %d, bfield: %d)\n', ...
            length(node_ids), size(coord,1), size(bfield,1));

    data.node_id = node_ids;
    data.x       = coord(ic, 2);    % [m]
    data.y       = coord(ic, 3);
    data.z       = coord(ic, 4);
    data.bx      = bfield(ib, 2);   % [T]
    data.by      = bfield(ib, 3);
    data.bz      = bfield(ib, 4);
    data.bsum    = bfield(ib, 5);
end
