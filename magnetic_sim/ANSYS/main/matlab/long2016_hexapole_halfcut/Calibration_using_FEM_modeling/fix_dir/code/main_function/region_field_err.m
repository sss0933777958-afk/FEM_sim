function errpct = region_field_err(coil, J)
%REGION_FIELD_ERR  Relative RMS field error over the sampled region [%].
%   errpct = REGION_FIELD_ERR(coil, J)
%   errpct = 100 * sqrt( J / sum||B_FEM||^2 ), where J is the fit cost
%   (sum of squared residuals). Matches the definition used by sweep_KI_radius.m.
    sumB2 = 0;
    for k = 1:numel(coil)
        sumB2 = sumB2 + sum(coil(k).bfem.^2);
    end
    errpct = 100*sqrt(J / sumB2);
end
