function write_KI_tex(fname, shape, R_um, I_A, K_bar, ell, ghat_I_B, errpct)
%WRITE_KI_TEX  Write a minimal, results-only LaTeX script for one fit.
%   WRITE_KI_TEX(fname, shape, R_um, I_A, K_bar, ell, ghat_I_B, errpct)
%   Emits ONLY the results: K-bar matrix, ell, ^Bg_I, and the relative RMS error.
%   [MODIFIED] load_coils is now flip-sink (field already all-source) -> K_bar diagonal is
%   already all-positive; NO presentation flip needed (was K_bar(:,[2 4 5])=-...). No prose.
    ge = floor(log10(abs(ghat_I_B)));  gm = ghat_I_B/10^ge;   % ^Bg_I -> mantissa x 10^exp

    fid = fopen(fname, 'w');
    fprintf(fid, '%% K_I fit -- shape=%s, R=%d um, I=%g A (all-source convention)\n', shape, R_um, I_A);
    fprintf(fid, '\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n\n');
    fprintf(fid, '\\[\n\\bar{\\mathbf{K}} =\n\\begin{bmatrix}\n');
    for i = 1:6
        fprintf(fid, '%.4f', K_bar(i,1));
        for j = 2:6, fprintf(fid, ' & %.4f', K_bar(i,j)); end
        if i < 6, fprintf(fid, ' \\\\\n'); else, fprintf(fid, '\n'); end
    end
    fprintf(fid, '\\end{bmatrix}\n\\]\n\n');
    fprintf(fid, '\\[\n\\widehat{\\ell} = %.1f~\\mu\\mathrm{m}, \\qquad {}^{B}\\widehat{g}_{I} = %.3f\\times10^{%d}~\\mathrm{mT/A}.\n\\]\n\n', ...
            ell, gm, ge);
    fprintf(fid, 'Relative RMS field error over the region: $%.2f\\%%$.\n\n', errpct);
    fprintf(fid, '\\end{document}\n');
    fclose(fid);
end
