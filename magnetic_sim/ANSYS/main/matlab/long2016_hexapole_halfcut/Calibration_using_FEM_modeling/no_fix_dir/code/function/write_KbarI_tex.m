function write_KbarI_tex(fname, shape, R_um, I_A, K_bar, ell, ghat_I_B, e_hat, coil_sign, errpct)
%WRITE_KBARI_TEX  Write a minimal, results-only LaTeX script for one bias fit.
%   WRITE_KBARI_TEX(fname, shape, R_um, I_A, K_bar, ell, ghat_I_B, e_hat, coil_sign, errpct)
%   Emits ONLY the results: K-bar matrix, ell, ^Bg_I, the bias e-hat (3x6), and the
%   relative RMS error. Applies the all-source display flip via coil_sign (upper poles
%   P2,P4,P5 -> diagonal all-positive). No prose.
%   e_hat is the 17-vector; e6z is reconstructed for the 3x6 display.
    Kdisp = K_bar .* coil_sign;                          % all-source presentation flip

    E36 = zeros(3,6);                                    % rebuild 3x6 bias (incl constrained e6z)
    E36(:,1) = e_hat(1:3);   E36(:,2) = e_hat(4:6);
    E36(:,3) = e_hat(7:9);   E36(:,4) = e_hat(10:12);
    E36(:,5) = e_hat(13:15);
    E36(1,6) = e_hat(16);    E36(2,6) = e_hat(17);
    E36(3,6) = e_hat(1) - e_hat(4) + e_hat(8) - e_hat(11) + e_hat(15);

    ge = floor(log10(abs(ghat_I_B)));  gm = ghat_I_B/10^ge;   % ^Bg_I -> mantissa x 10^exp

    fid = fopen(fname, 'w');
    fprintf(fid, '%% no_fix_l 18-param bias fit -- shape=%s, R=%d um, I=%g A (all-source convention)\n', shape, R_um, I_A);
    fprintf(fid, '\\documentclass{article}\n\\usepackage{amsmath}\n\\begin{document}\n\n');
    % --- K-bar ---
    fprintf(fid, '\\[\n\\bar{\\mathbf{K}} =\n\\begin{bmatrix}\n');
    for i = 1:6
        fprintf(fid, '%.4f', Kdisp(i,1));
        for j = 2:6, fprintf(fid, ' & %.4f', Kdisp(i,j)); end
        if i < 6, fprintf(fid, ' \\\\\n'); else, fprintf(fid, '\n'); end
    end
    fprintf(fid, '\\end{bmatrix}\n\\]\n\n');
    % --- ell, ^Bg_I ---
    fprintf(fid, '\\[\n\\widehat{\\ell} = %.1f~\\mu\\mathrm{m}, \\qquad {}^{B}\\widehat{g}_{I} = %.3f\\times10^{%d}~\\mathrm{mT/A}.\n\\]\n\n', ell, gm, ge);
    % --- bias e-hat (3x6, columns e1..e6) ---
    fprintf(fid, '\\[\n\\widehat{\\mathbf{e}} =\n\\begin{bmatrix}\n');
    for r = 1:3
        fprintf(fid, '%+.4f', E36(r,1));
        for k = 2:6, fprintf(fid, ' & %+.4f', E36(r,k)); end
        if r < 3, fprintf(fid, ' \\\\\n'); else, fprintf(fid, '\n'); end
    end
    fprintf(fid, '\\end{bmatrix}\n\\]\n\n');
    fprintf(fid, 'Relative RMS field error over the region: $%.2f\\%%$.\n\n', errpct);
    fprintf(fid, '\\end{document}\n');
    fclose(fid);
end
