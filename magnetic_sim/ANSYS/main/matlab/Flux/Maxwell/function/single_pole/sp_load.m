function L = sp_load(cfg)
%SP_LOAD  Read the single-pole .fld into a grid in the tip frame (mT, all-source) with an iron mask.
%   L = SP_LOAD(cfg)
%   Frame: origin = cfg.WP (tip apex), axes = Maxwell global axes (pole axis is already +x).
%   Returns
%     .xs .ys .zs     grid vectors [m] (tip frame)
%     .B              nx x ny x nz x 3 field [mT], multiplied by cfg.s_source
%     .iron           nx x ny x nz logical, true = node inside the pole plate
%     .interp(P)      trilinear field at P (Np x 3 [m]) -> Np x 3 [mT]
%     .air_stencil(P) true where all 8 grid nodes around P are air (trilinear value is pure air field)
%   The iron mask models the plate near the tip only (fillet + tip wedge); it errors out if the
%   .fld box reaches past cfg.TIP_WEDGE_XMAX, where the outline changes.

    f = fullfile(cfg.fld_dir, cfg.fld_files{1});
    fid = fopen(f, 'r');
    if fid < 0, error('sp_load:open', 'cannot open %s', f); end
    fgetl(fid); fgetl(fid);                                   % two header lines
    C = textscan(fid, '%f %f %f %f %f %f', 'CollectOutput', true);
    fclose(fid);
    D = C{1};
    D(:,1:3) = D(:,1:3) * 1e-3;                               % mm -> m (Maxwell global)

    xs = unique(D(:,1));  ys = unique(D(:,2));  zs = unique(D(:,3));
    nx = numel(xs);  ny = numel(ys);  nz = numel(zs);
    assert(nx*ny*nz == size(D,1), 'sp_load: %s is not a full grid', f);
    assert(D(2,3) ~= D(1,3) && D(2,1) == D(1,1), 'sp_load: expected z to vary fastest in %s', f);
    toG = @(v) permute(reshape(v, nz, ny, nx), [3 2 1]);

    % iron mask in the Maxwell global frame
    [X, Y, Z] = ndgrid(xs, ys, zs);
    assert(max(xs) < cfg.TIP_WEDGE_XMAX, ...
           'sp_load: .fld box reaches x = %.3f mm, beyond the tip wedge (%.3f mm)', max(xs)*1e3, cfg.TIP_WEDGE_XMAX*1e3);
    th = cfg.TIP_HALF_ANG * pi/180;
    r  = cfg.POLE_TIP_R;
    xc = cfg.pole_tip(1) + r;                                 % fillet centre on the axis
    xv = xc - r/sin(th);                                      % virtual wedge apex
    xt = xc - r*sin(th);                                      % fillet-to-wedge tangent x
    yc = cfg.pole_tip(2);
    inz   = Z >= cfg.PLATE_Z(1) & Z <= cfg.PLATE_Z(2);
    wedge = X >= xt & abs(Y - yc) <= (X - xv) * tan(th);
    ball  = (X - xc).^2 + (Y - yc).^2 <= r^2;
    iron  = inz & (wedge | ball);

    sgn = cfg.s_source(1);
    B = zeros(nx, ny, nz, 3);
    for c = 1:3, B(:,:,:,c) = sgn * 1e3 * toG(D(:,3+c)); end   % T -> mT, all-source

    % shift to the tip frame
    L.xs = xs - cfg.WP(1);  L.ys = ys - cfg.WP(2);  L.zs = zs - cfg.WP(3);
    L.B = B;  L.iron = iron;  L.fld = f;
    Fi = cell(1,3);
    for c = 1:3, Fi{c} = griddedInterpolant({L.xs, L.ys, L.zs}, B(:,:,:,c), 'linear', 'none'); end
    L.interp = @(P) [Fi{1}(P), Fi{2}(P), Fi{3}(P)];
    L.air_stencil = @(P) stencil_air(P, L);
    fprintf('  sp_load: %d x %d x %d nodes, %d iron  (%s)\n', nx, ny, nz, nnz(iron), f);
end

function ok = stencil_air(P, L)
    ok = false(size(P,1), 1);
    for i = 1:size(P,1)
        ix = find(L.xs <= P(i,1), 1, 'last');  iy = find(L.ys <= P(i,2), 1, 'last');  iz = find(L.zs <= P(i,3), 1, 'last');
        if isempty(ix) || isempty(iy) || isempty(iz) || ix == numel(L.xs) || iy == numel(L.ys) || iz == numel(L.zs)
            continue                                          % outside the grid
        end
        blk = L.iron(ix:ix+1, iy:iy+1, iz:iz+1);
        ok(i) = ~any(blk(:));
    end
end
