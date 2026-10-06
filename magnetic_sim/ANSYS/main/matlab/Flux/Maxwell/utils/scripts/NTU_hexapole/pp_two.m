% pp_two.m -- test only: two-charge 1-D model on pair_pole P1 excitation, 40 axis samples.
% H(x) = a1/(x-b1)^2 + a2/(x-b2)^2 (right charge b1 > -20 um, left charge b2 < -800 um), a1, a2 closed form.
FLD = 'D:\Maxwell_sim\NTU\export\pair_pole\P1.fld';  TIP = [42.42464, 0, 0.125]*1e-3;
fid = fopen(FLD); fgetl(fid); fgetl(fid); C = textscan(fid, '%f %f %f %f %f %f', 'CollectOutput', true); fclose(fid);
D = C{1};  D(:,1:3) = D(:,1:3)*1e-3;
xs = unique(D(:,1))-TIP(1); ys = unique(D(:,2))-TIP(2); zs = unique(D(:,3))-TIP(3);
F = griddedInterpolant({xs,ys,zs}, permute(reshape(D(:,4)*1e3, numel(zs), numel(ys), numel(xs)), [3 2 1]));
xi = -(1:40).'*20;  y = abs(F(xi*1e-6, zeros(40,1), zeros(40,1)));
fprintf('|b_x| [mT] at d = 500, 700, 800 um: %.2f %.2f %.2f (min over samples %.2f at d = %d)\n', y(25), y(35), y(40), min(y), -xi(y==min(y)));
U  = @(p) [1./(xi-p(1)).^2, 1./(xi-p(2)).^2];
Jf = @(p) sum((y - U(p)*(U(p)\y)).^2);
best = inf;
for b1 = 0:20:400, for b2 = -1600:20:-820
    J = Jf([b1 b2]); if J < best, best = J; p0 = [b1 b2]; end
end, end
p = fminsearch(Jf, p0, optimset('TolX',1e-8,'TolFun',1e-10,'MaxFunEvals',1e5));
a = U(p)\y;  r = y - U(p)*a;
fprintf('two charges: b1 = %.2f um, b2 = %.2f um (left tip at -1000), a1 = %.4f, a2 = %.4f mT*mm^2/A | cost %.1f  RMS %.3f mT\n', ...
        p, a*1e-6, sum(r.^2), sqrt(mean(r.^2)));
for k = 1:3:40, id = k:min(k+2,40); rm = sqrt(mean(r(id).^2));
    fprintf('  %4d-%-4d um  RMS %.3f  rel %.2f%%\n', -xi(id(1)), -xi(id(end)), rm, rm/mean(y(id))*100); end
