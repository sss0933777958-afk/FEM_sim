function T = emit_tex()
%EMIT_TEX  Shared LaTeX table emitters for calibration result PDFs.
%   Returns a struct of function handles so the same helpers are reused by
%   emit_model_results / emit_charge_Q / gen_B_matrix (no per-file duplication).
%     T.mat(fid,name,M,lab,fmode)                       6x6 matrix (auto 10^n if |exp|>=2)
%     T.e(fid,name,E,collab,fmode)                      3x6 bias matrix (rows e_x/e_y/e_z)
%     T.scalar_unit(fid,name,val,unit)                  scalar w/ unit (10^0/10^1 not factored)
%     T.labeled(fid,name,M,rowlab,collab,fmode,caption) labeled 6x6 + caption line
%   Numeric-annotation convention per figure-style.md (only |exp|>=2 factored out).
    T.mat         = @emit_mat;
    T.e           = @emit_e;
    T.scalar_unit = @emit_scalar_unit;
    T.labeled     = @emit_labeled_matrix;
end

function emit_mat(fid, name, M, lab, fmode)
    Ms = M; fac = '';
    if strcmp(fmode,'auto'), mx=max(abs(M(:))); if mx>0, e=floor(log10(mx)); if abs(e)>=2, Ms=M/10^e; fac=sprintf('10^{%d}\\,',e); end, end, end
    fprintf(fid,'\\[\n%s = %s\\begin{array}{c|cccccc}\n ', name, fac);
    for j=1:6, fprintf(fid,'& %s ', lab{j}); end
    fprintf(fid,'\\\\\\hline\n');
    for i=1:6, fprintf(fid,'%s ', lab{i}); for j=1:6, fprintf(fid,'& %+9.4f ', Ms(i,j)); end, fprintf(fid,'\\\\\n'); end
    fprintf(fid,'\\end{array}\n\\]\n');
end

function emit_e(fid, name, E, collab, fmode)
    Es = E; fac = '';
    if strcmp(fmode,'auto'), mx=max(abs(E(:))); if mx>0, e=floor(log10(mx)); if abs(e)>=2, Es=E/10^e; fac=sprintf('10^{%d}\\,',e); end, end, end
    erow = {'e_x','e_y','e_z'};
    fprintf(fid,'\\[\n%s = %s\\begin{array}{c|cccccc}\n ', name, fac);
    for j=1:6, fprintf(fid,'& %s ', collab{j}); end
    fprintf(fid,'\\\\\\hline\n');
    for i=1:3, fprintf(fid,'%s ', erow{i}); for j=1:6, fprintf(fid,'& %+9.4f ', Es(i,j)); end, fprintf(fid,'\\\\\n'); end
    fprintf(fid,'\\end{array}\n\\]\n');
end

function emit_scalar_unit(fid, name, val, unit)
    ge = floor(log10(abs(val)));
    if abs(ge)<=1
        fprintf(fid,'\\[\n %s = %.4f~\\mathrm{%s}\n\\]\n', name, val, unit);
    else
        fprintf(fid,'\\[\n %s = %.3f\\times10^{%d}~\\mathrm{%s}\n\\]\n', name, val/10^ge, ge, unit);
    end
end

function emit_labeled_matrix(fid, name_tex, M, rowlab, collab, factor_mode, caption)
    if strcmp(factor_mode,'auto')
        mx = max(abs(M(:)));
        if mx > 0, e = floor(log10(mx)); else, e = 0; end
        if abs(e) >= 2, Ms = M / 10^e; fac = sprintf('10^{%d}\\,', e); else, Ms = M; fac = ''; end
    else
        Ms = M; fac = '';
    end
    fprintf(fid,'\\[\n%s = %s\\begin{array}{c|cccccc}\n', name_tex, fac);
    fprintf(fid,' ');
    for j = 1:6, fprintf(fid,'& %s ', collab{j}); end
    fprintf(fid,'\\\\\\hline\n');
    for i = 1:6
        fprintf(fid,'%s ', rowlab{i});
        for j = 1:6, fprintf(fid,'& %+9.4f ', Ms(i,j)); end
        fprintf(fid,'\\\\\n');
    end
    fprintf(fid,'\\end{array}\n\\]\n');
    fprintf(fid,'\\noindent\\small %s\n\n', caption);
end
