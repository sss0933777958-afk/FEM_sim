function T = emit_tex()
%EMIT_TEX  低階 LaTeX 排版 helper（供 emit_results 等用）。回 handle struct：
%   T.mat(fid,name,M,collab,fmode)   6×6 矩陣（含欄/列標籤）
%   T.e(fid,name,E,collab,fmode)     3×6 偏移矩陣（列標 e_x/e_y/e_z）
%   T.scalar_unit(fid,name,val,unit) 純量 + 單位
%   數值標註慣例（figure-style）：矩陣 |指數|≥2 才抽 ×10ⁿ；純量 |指數|≥2 才抽，否則印原值。
    T.mat         = @emit_mat;
    T.e           = @emit_e;
    T.scalar_unit = @emit_scalar_unit;
end

function emit_mat(fid, name, M, collab, fmode)
    if nargin < 5, fmode = ''; end
    fac = '';  Ms = M;
    if strcmp(fmode,'auto')
        mx = max(abs(M(:)));
        if mx > 0
            ex = floor(log10(mx));
            if abs(ex) >= 2, Ms = M/10^ex; fac = sprintf('10^{%d}\\,', ex); end
        end
    end
    [nr, nc] = size(Ms);                                         % [MODIFIED] 任意尺寸(原寫死 6×6);6×6 呼叫輸出不變
    fprintf(fid, '\\[\n%s = %s\\begin{bmatrix}\n', name, fac);   % 標準矩陣(不用欄位/標籤表格);正值不標 + 號
    for i = 1:nr
        for j = 1:nc
            if j > 1, fprintf(fid, ' & '); end
            fprintf(fid, '%9.4f', Ms(i,j));
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\end{bmatrix}\n\\]\n');
end

function emit_e(fid, name, E, collab, fmode)
    if nargin < 5, fmode = ''; end
    fac = '';  Es = E;
    if strcmp(fmode,'auto')
        mx = max(abs(E(:)));
        if mx > 0
            ex = floor(log10(mx));
            if abs(ex) >= 2, Es = E/10^ex; fac = sprintf('10^{%d}\\,', ex); end
        end
    end
    rowlab = {'e_x','e_y','e_z'};
    fprintf(fid, '\\[\n%s = %s\\begin{array}{c|cccccc}\n', name, fac);
    fprintf(fid, ' & %s & %s & %s & %s & %s & %s \\\\ \\hline\n', collab{:});
    for i = 1:3
        fprintf(fid, '%s', rowlab{i});
        for j = 1:6, fprintf(fid, ' & %9.4f', Es(i,j)); end   % 正值不標 + 號
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\end{array}\n\\]\n');
end

function emit_scalar_unit(fid, name, val, unit)
    ge = floor(log10(abs(val)));
    if abs(ge) <= 1
        fprintf(fid, '\\[ %s = %.4f~\\mathrm{%s} \\]\n', name, val, unit);
    else
        fprintf(fid, '\\[ %s = %.3f\\times10^{%d}~\\mathrm{%s} \\]\n', name, val/10^ge, ge, unit);
    end
end
