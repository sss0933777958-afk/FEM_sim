function make_scaled_coil1(scale, outname)
% MAKE_SCALED_COIL1  Linearly scale coil1/standard FEM data to a new current.
%   make_scaled_coil1(2,'standard_2A')  -> data/coil1/standard_2A (B x2)
%   make_scaled_coil1(3,'standard_3A')  -> data/coil1/standard_3A (B x3)
%   Linear FEM (mu_r=280 const) -> B is exactly proportional to coil current,
%   so 1A baseline x scale = (scale)A field. NO re-solve needed.
%   coord_*.dat copied verbatim (node positions unchanged); bfield_*.dat has
%   its B columns (Bx,By,Bz,Bsum) multiplied by `scale`, written space-separated
%   so import_ansys_data() reads it back identically.

    rr  = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\data\coil1';
    src = fullfile(rr,'standard');
    dst = fullfile(rr,outname);
    if ~exist(dst,'dir'); mkdir(dst); end

    for ds = {'all','wp'}
        d = ds{1};
        % --- coord: copy verbatim ---
        copyfile(fullfile(src,sprintf('coil1_coord_%s.dat',d)), ...
                 fullfile(dst,sprintf('coil1_coord_%s.dat',d)));
        % --- bfield: parse (same as import_ansys_data), scale B cols, rewrite ---
        bf = fullfile(src,sprintf('coil1_bfield_%s.dat',d));
        txt = fileread(bf);
        txt = regexprep(txt,'(\d)([-+])','$1 $2');     % split concatenated negatives
        tmp = [tempname '.dat']; fid=fopen(tmp,'w'); fprintf(fid,'%s',txt); fclose(fid);
        M = readmatrix(tmp,'FileType','text'); delete(tmp);
        nb = min(5,size(M,2));
        M = M(~any(isnan(M(:,1:nb)),2),:);             % keep valid rows (node,Bx,By,Bz,Bsum)
        M(:,2:5) = M(:,2:5) * scale;                   % linear current scaling
        % write: node(int) + Bx By Bz Bsum (scientific), space separated
        out = fullfile(dst,sprintf('coil1_bfield_%s.dat',d));
        fid = fopen(out,'w');
        fprintf(fid,'%10d %18.9E %18.9E %18.9E %18.9E\n', M(:,1:5).');
        fclose(fid);
        fprintf('  %s: %d rows scaled x%g -> %s\n', d, size(M,1), scale, out);
    end

    % --- verify round-trip ---
    addpath('G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\backup\hexapole-long2016\analysis');
    q = import_ansys_data(dst,'all','coil1');
    bmax = max(sqrt(q.bx.^2+q.by.^2+q.bz.^2));
    fprintf('VERIFY %s: matched %d nodes, |B|max = %.4g T (baseline x%g)\n', ...
            outname, numel(q.x), bmax, scale);
end
