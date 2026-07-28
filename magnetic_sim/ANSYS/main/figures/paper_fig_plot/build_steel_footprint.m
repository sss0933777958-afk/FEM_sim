function build_steel_footprint()
% 一次性:讀 ANSYS NWRITE 鐵件節點(MAT 2)→ 存**鐵件節點 ID**(steel_ids.mat,供精確比對場節點)。
%   (順便存 xy footprint mask 備用。)NWRITE 有負號黏連,先 regexprep 補空白。
    src = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\ANSYS_data\long2016_hexapole_halfcut\db\mesh\_steel_export\steel_nodes.dat';
    here = fileparts(mfilename('fullpath'));
    txt = fileread(src);  txt = regexprep(txt, '(\d)([-+])', '$1 $2');
    tmp = [tempname '.dat'];  fid=fopen(tmp,'w'); fprintf(fid,'%s',txt); fclose(fid);
    M = readmatrix(tmp,'FileType','text');  delete(tmp);
    M = M(~any(isnan(M(:,1:3)),2), :);
    steel_ids = M(:,1);                                   % 鐵件節點 ID
    save(fullfile(here,'steel_ids.mat'), 'steel_ids');
    fprintf('steel node IDs saved: %d\n', numel(steel_ids));
end
