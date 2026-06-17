function [d, BX, BY, BZ, BSUM] = load_prpath(filepath)
% LOAD_PRPATH  Parse ANSYS POST1 PRPATH output file
%   [d, BX, BY, BZ, BSUM] = load_prpath(filepath)
%   Returns column vectors: distance [m], B components [T], magnitude [T].
%   Handles ANSYS fixed-width format with concatenated negatives.

    txt = fileread(filepath);
    % Fix ANSYS fixed-width: insert space before +/- that follows a digit
    txt = regexprep(txt, '(\d)([-+])', '$1 $2');
    tmpfile = [tempname '.dat'];
    fid = fopen(tmpfile, 'w');
    fprintf(fid, '%s', txt);
    fclose(fid);
    M = readmatrix(tmpfile, 'FileType', 'text');
    delete(tmpfile);
    M = M(~any(isnan(M), 2), :);

    d    = M(:,1);
    BX   = M(:,2);
    BY   = M(:,3);
    BZ   = M(:,4);
    BSUM = M(:,5);
end
