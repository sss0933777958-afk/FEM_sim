% pp_n40.m -- pair_pole P1 excitation, 1-D model H = a/(x-b)^2 on |b_x| [mT], 17 vs 40 axis samples (20 um step).
FLD = 'D:\Maxwell_sim\NTU\export\pair_pole\P1.fld';  TIP = [42.42464, 0, 0.125]*1e-3;
fid = fopen(FLD); fgetl(fid); fgetl(fid); C = textscan(fid, '%f %f %f %f %f %f', 'CollectOutput', true); fclose(fid);
D = C{1};  D(:,1:3) = D(:,1:3)*1e-3;
xs = unique(D(:,1))-TIP(1); ys = unique(D(:,2))-TIP(2); zs = unique(D(:,3))-TIP(3);
nx = numel(xs); ny = numel(ys); nz = numel(zs);
F = griddedInterpolant({xs,ys,zs}, permute(reshape(D(:,4)*1e3, nz, ny, nx), [3 2 1]));
for NP = [17 40]
    xi = -(1:NP).' * 20;  y = abs(F(xi*1e-6, zeros(NP,1), zeros(NP,1)));     % |b_x| [mT] (I = 1 A)
    af = @(b) sum(y./(xi-b).^2) / sum(1./(xi-b).^4);
    Jf = @(b) sum((y - af(b)./(xi-b).^2).^2);
    bs = linspace(max(xi)+0.01, 3000, 300001);  [~, im] = min(arrayfun(Jf, bs));
    b = fminbnd(Jf, bs(max(im-1,1)), bs(min(im+1,end)), optimset('TolX',1e-9));  a = af(b);
    r = y - a./(xi-b).^2;
    fprintf('\nN = %d (d = 20..%d um): b = %.2f um  a = %.4f mT*mm^2/A  cost = %.1f  RMS = %.3f mT\n', ...
            NP, -xi(end), b, a*1e-6, sum(r.^2), sqrt(mean(r.^2)));
    fprintf('  %7s %10s %10s %10s\n', 'points', 'd[um]', 'RMS', 'rel');
    for k = 1:3:NP
        id = k:min(k+2,NP);  rm = sqrt(mean(r(id).^2));
        fprintf('  %3d-%-3d %4d-%-5d %10.3f %9.2f%%\n', id(1), id(end), -xi(id(1)), -xi(id(end)), rm, rm/mean(y(id))*100);
    end
end
