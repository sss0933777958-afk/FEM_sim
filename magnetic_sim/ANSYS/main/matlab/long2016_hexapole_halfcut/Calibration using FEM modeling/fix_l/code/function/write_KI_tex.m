function write_KI_tex(fname, shape, R_um, I_A, Khat, ell, gB, errpct)
%WRITE_KI_TEX  Write a minimal, results-only LaTeX script for one fit.
%   WRITE_KI_TEX(fname, shape, R_um, I_A, Khat, ell, gB, errpct)
%   Emits ONLY the results: Khat_I^FEM matrix, ell, gB, and the relative RMS error.
%   Applies the all-source presentation flip (upper poles P2,P4,P5 = columns 2,4,5)
%   so the diagonal is all-positive, per the project sign convention. No prose.
    Khat(:, [2 4 5]) = -Khat(:, [2 4 5]);            % all-source presentation flip

    ge = floor(log10(abs(gB)));  gm = gB/10^ge;       % gB -> mantissa x 10^exp

    fid = fopen(fname, 'w');
    fprintf(fid, '%% K_I fit -- shape=%s, R=%d um, I=%g A (all-source convention)\n', shape, R_um, I_A);
    fprintf(fid, '\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n\n');
    fprintf(fid, '\\[\n\\widehat{\\mathbf{K}}_I^{\\mathrm{FEM}} =\n\\begin{bmatrix}\n');
    for i = 1:6
        fprintf(fid, '%.4f', Khat(i,1));
        for j = 2:6, fprintf(fid, ' & %.4f', Khat(i,j)); end
        if i < 6, fprintf(fid, ' \\\\\n'); else, fprintf(fid, '\n'); end
    end
    fprintf(fid, '\\end{bmatrix}\n\\]\n\n');
    fprintf(fid, '\\[\n\\widehat{\\ell} = %.3f~\\mathrm{mm}, \\qquad \\widehat{g}_B = %.3f\\times10^{%d}.\n\\]\n\n', ...
            ell*1e3, gm, ge);
    fprintf(fid, 'Relative RMS field error over the region: $%.2f\\%%$.\n\n', errpct);
    fprintf(fid, '\\end{document}\n');
    fclose(fid);
end
