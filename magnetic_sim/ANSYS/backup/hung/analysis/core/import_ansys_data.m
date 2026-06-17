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

    % --- Read coordinate file (filter to numeric data lines) ---
    txt_c = fileread(coord_file);
    lines_c = strsplit(txt_c, '\n');
    keep_c = cellfun(@(s) ~isempty(regexp(s, '^\s*\d', 'once')), lines_c);
    txt_c_clean = strjoin(lines_c(keep_c), '\n');
    tmpfile_c = [tempname '.dat'];
    fid = fopen(tmpfile_c, 'w'); fprintf(fid, '%s', txt_c_clean); fclose(fid);
    coord = readmatrix(tmpfile_c, 'FileType', 'text');
    delete(tmpfile_c);
    coord = coord(~any(isnan(coord), 2), :);
    if size(coord, 2) >= 4
        coord = coord(:, 1:4);
    end

    % --- Read B-field file (fixed-width format, negatives may concatenate) ---
    % ANSYS PRNSOL uses fixed-width columns; consecutive negative values
    % like "-0.123E-002-0.456E-003" have no space separator.
    % Fix: insert space before any '-' or '+' that immediately follows a digit.
    % Also: Hung post scripts produce MAPDL banner headers — filter to
    %       lines starting with a number (node ID).
    txt = fileread(bfield_file);
    txt = regexprep(txt, '(\d)([-+])', '$1 $2');
    lines = strsplit(txt, '\n');
    % Keep only lines that start with optional spaces then a digit (data rows)
    keep = cellfun(@(s) ~isempty(regexp(s, '^\s*\d', 'once')), lines);
    txt_clean = strjoin(lines(keep), '\n');
    tmpfile = [tempname '.dat'];
    fid = fopen(tmpfile, 'w'); fprintf(fid, '%s', txt_clean); fclose(fid);
    bfield = readmatrix(tmpfile, 'FileType', 'text');
    delete(tmpfile);
    bfield = bfield(~any(isnan(bfield), 2), :);
    % Keep only rows with exactly 5 columns (node, BX, BY, BZ, BSUM)
    if size(bfield, 2) >= 5
        bfield = bfield(:, 1:5);
    end

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
