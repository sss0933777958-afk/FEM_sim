function errpct = region_field_err(Bstack, J)
%REGION_FIELD_ERR  Relative RMS field error over the sampled region [%].
%   errpct = REGION_FIELD_ERR(Bstack, J)
%   errpct = 100 * sqrt( J / sum||B_FEM||^2 ), J = sum of squared residuals,
%   sum over all stacked simulations. Same definition as the fix_l driver.
    sumB2  = sum(Bstack(:).^2);
    errpct = 100*sqrt(J / sumB2);
end
