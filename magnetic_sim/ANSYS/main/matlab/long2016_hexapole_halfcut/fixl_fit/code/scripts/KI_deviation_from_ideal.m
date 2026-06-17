%% KI_deviation_from_ideal.m
%  Which ball radius gives K_hat_I^FEM closest to the ideal K_hat = eye - ones/6?
%  Uses the source-corrected K_hat (upper-pole cols 2,4,5 flipped), Frobenius norm.
data = 'G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\MATLAB_data\long2016_hexapole_halfcut\charge_fit\fit_KI_ball';
Kid  = eye(6) - ones(6)/6;  nid = norm(Kid,'fro');
R = 50:50:500;  dev = zeros(size(R));
fprintf('  R[um]   ||Khat-Kideal||_F   rel[%%]\n');
for i = 1:numel(R)
    S = load(fullfile(data,sprintf('fit_KI_R%03d.mat',R(i))),'Khat');
    K = S.Khat;  K(:,[2 4 5]) = -K(:,[2 4 5]);
    dev(i) = norm(K-Kid,'fro');
    fprintf('  %4d      %8.4f          %6.2f\n', R(i), dev(i), dev(i)/nid*100);
end
[~,j] = min(dev);
fprintf('\nCLOSEST to ideal: R = %d um  (dev = %.4f, rel = %.2f%%)\n', R(j), dev(j), dev(j)/nid*100);
