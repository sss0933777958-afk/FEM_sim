%% xa_interp.m -- b on the x_a axis, straight from the FEM export by OUR trilinear
%%                 interpolation (not from any fitted model), all six excitations.
%
%  POINTS.  1001 equidistant points on the x_a axis of the ACTUATOR frame,
%  s = -150 .. +150 um, i.e. P_act = (s, 0, 0).  Their actuator-frame coordinates
%  are what the next step needs, so they are what is stored; the measure-frame
%  coordinates are kept alongside only because that is what indexes the lattice.
%
%      P_meas = P_act * R_act  +  [0 0 SPH_OFST]        (row-vector form)
%      B_act  = B_meas * R_act'
%
%  INTERPOLATION.  The project's own regular-grid trilinear -- trilerp and the
%  lattice builder are reproduced here line for line from
%  Flux/Maxwell/function/conv_design_ws.m rather than shared, which is the standing
%  rule for this engine (no common loader; each caller carries its own copy).
%
%  SIGN.  cfg.s_source is all +1 for this model (the Maxwell export is already
%  all-source), so nothing is flipped and this b is the same convention the RBF and
%  spherical-harmonic fits were given.
%
%  Output: data/long2016_hexapole_halfcut/.mat/xa_interp.mat
%
%  [ADDED 2026-09-17]

clearvars;  clc;
HERE = fileparts(mfilename('fullpath'));
FMX  = fileparts(fileparts(HERE));
MAIN = fileparts(fileparts(fileparts(FMX)));
FLX  = fullfile(MAIN,'matlab','Flux','Maxwell');
addpath(fullfile(FMX,'function'), fullfile(FMX,'function'), ...
        fullfile(FLX,'function'), fullfile(FLX,'common_path'), fullfile(FLX,'utils'));
DAT = fullfile(FMX,'utils','data');

MODEL='long2016_hexapole_halfcut';  GEOM='tip40um';  VARIANT='maxwell';
SPAN = 150;      % |s| <= SPAN [um]
NQ   = 1001;     % equidistant points on the axis
AXL  = 1:3;      % all three actuator axes (1 = x_a, 2 = y_a, 3 = z_a)
NODETOL = 1e-6;  % a sample counts as sitting ON a node if every fractional lattice
                 % coordinate is within this many CELLS of an integer (2e-11 m)

ANM  = {'x_a','y_a','z_a'};
cfg  = model_config(MODEL, GEOM);
raw  = extract_maxwell_data(cfg, 'all', VARIANT);        % [T], Maxwell/measure frame
RACT = cfg.R_act;
assert(all(cfg.s_source == 1), 'xa_interp:sign', 'this script assumes an all-source export');

%% ---- the points, defined in the ACTUATOR frame -------------------------------
s = linspace(-SPAN, SPAN, NQ).';

%% ---- the lattice (reproduced from conv_design_ws.m) --------------------------
N   = numel(raw.x);   N_I = size(raw.B,3);
q   = @(vv) unique(round(vv*1e9))/1e9;                   % quantise to nm
xg  = q(raw.x);   yg = q(raw.y);   zg = q(raw.z);
n   = [numel(xg) numel(yg) numel(zg)];
assert(prod(n) == N, 'xa_interp:notGrid', '.fld is not a full regular grid');
o   = [xg(1) yg(1) zg(1)];
h   = [(xg(end)-xg(1))/(n(1)-1), (yg(end)-yg(1))/(n(2)-1), (zg(end)-zg(1))/(n(3)-1)];
ix  = round((raw.x-o(1))/h(1))+1;  iy = round((raw.y-o(2))/h(2))+1;  iz = round((raw.z-o(3))/h(3))+1;
dev = max([max(abs(o(1)+(ix-1)*h(1)-raw.x)), max(abs(o(2)+(iy-1)*h(2)-raw.y)), ...
           max(abs(o(3)+(iz-1)*h(3)-raw.z))]);
assert(dev < 1e-9, 'xa_interp:notUniform', 'grid is not uniform (max dev %.3e m)', dev);
lin = sub2ind(n, ix, iy, iz);
V   = zeros(prod(n), 3*N_I);
for m = 1:N_I
    for cc = 1:3, V(lin,(m-1)*3+cc) = raw.B(:,cc,m);  end
end
lat = struct('n',n, 'o',o, 'h',h, 'N_I',N_I, 'V',V);
fprintf(['  [lattice] %d x %d x %d, pitch %.1f/%.1f/%.1f um' newline], n, h*1e6);

%% ---- one axis at a time -------------------------------------------------------
P_act = cell(1,3);   P_glb = cell(1,3);   B = cell(1,3);   CHK = struct([]);
for a = AXL
    Pa = zeros(NQ,3);   Pa(:,a) = s;                     % [um], actuator frame
    Pg = (Pa*1e-6) * RACT;                               % -> measure [m]
    Pg(:,3) = Pg(:,3) + cfg.SPH_OFST;                    % -> global z
    [Bt, inbox, loc] = trilerp_(lat, Pg(:,1), Pg(:,2), Pg(:,3));
    assert(all(inbox), 'xa_interp:outbox', 'axis %d: %d points outside the box', a, nnz(~inbox));
    Bm = 1e3 * Bt;                                       % T -> mT, measure frame
    Ba = zeros(NQ,3,N_I);
    for m = 1:N_I, Ba(:,:,m) = Bm(:,:,m) * RACT.';  end   % -> actuator frame
    P_act{a} = Pa;   P_glb{a} = Pg;   B{a} = Ba;

    % ---- do any samples sit ON a lattice node? --------------------------------
    % A node is where all three fractional lattice coordinates are integers.  frac
    % is the distance to the nearest integer, in cells; d_node converts it to um.
    F    = [loc.i0 + loc.tx, loc.j0 + loc.ty, loc.k0 + loc.tz];
    frac = abs(F - round(F));
    onnode = all(frac < NODETOL, 2);
    d_node = sqrt(sum((frac .* (h*1e6)).^2, 2));         % [um] to the nearest node
    onface = sum(frac < NODETOL, 2);                     % 1 = on a face, 2 = on an edge
    CHK(a).onnode = onnode;   CHK(a).d_node = d_node;   CHK(a).onface = onface;
    CHK(a).n_on   = nnz(onnode);
    fprintf([newline '%s: %d of %d samples land ON a node%s' newline], ANM{a}, ...
            nnz(onnode), NQ, repmat(' ', 1, 0));
    if any(onnode)
        fprintf(['    at s = %s um' newline], mat2str(s(onnode).', 6));
    end
    fprintf(['    distance to nearest node: min %.4f, median %.4f, max %.4f um' newline], ...
            min(d_node(~onnode)), median(d_node), max(d_node));
    fprintf(['    samples with >=1 integer coordinate (on a cell face): %d;  >=2 (edge): %d' newline], ...
            nnz(onface>=1), nnz(onface>=2));
    fprintf(['    distinct cells %d, distinct nodes read %d' newline], ...
            size(unique([loc.i0 loc.j0 loc.k0],'rows'),1), ...
            size(unique([loc.i0 loc.j0 loc.k0; loc.i0+1 loc.j0 loc.k0; ...
                         loc.i0 loc.j0+1 loc.k0; loc.i0 loc.j0 loc.k0+1; ...
                         loc.i0+1 loc.j0+1 loc.k0; loc.i0+1 loc.j0 loc.k0+1; ...
                         loc.i0 loc.j0+1 loc.k0+1; loc.i0+1 loc.j0+1 loc.k0+1],'rows'),1));
    for m = 1:N_I
        fprintf(['    P%d: |b| %8.4f .. %8.4f mT' newline], m, ...
                min(vecnorm(Ba(:,:,m),2,2)), max(vecnorm(Ba(:,:,m),2,2)));
    end
end

save(fullfile(DAT,'xa_interp.mat'), 's','P_act','P_glb','B','CHK','AXL','SPAN','NQ', ...
     'NODETOL','MODEL','GEOM','VARIANT','RACT','lat');
fprintf([newline 'wrote %s' newline], fullfile(DAT,'xa_interp.mat'));

%% ---- local: trilerp, reproduced from conv_design_ws.m ------------------------
function [Bq, ok, loc] = trilerp_(lat, xq, yq, zq)
    nx = lat.n(1);   ny = lat.n(2);   nz = lat.n(3);
    fx = (xq - lat.o(1)) / lat.h(1);
    fy = (yq - lat.o(2)) / lat.h(2);
    fz = (zq - lat.o(3)) / lat.h(3);
    TOL = 1e-9;
    ok = fx >= -TOL & fx <= nx-1+TOL & fy >= -TOL & fy <= ny-1+TOL ...
       & fz >= -TOL & fz <= nz-1+TOL;
    fx = fx(ok);   fy = fy(ok);   fz = fz(ok);
    i0 = min(max(floor(fx), 0), nx-2);   tx = min(max(fx - i0, 0), 1);
    j0 = min(max(floor(fy), 0), ny-2);   ty = min(max(fy - j0, 0), 1);
    k0 = min(max(floor(fz), 0), nz-2);   tz = min(max(fz - k0, 0), 1);
    a = i0 + 1;   b = j0 + 1;   c = k0 + 1;
    sz = [nx ny nz];
    L  = @(da,db,dc) sub2ind(sz, a+da, b+db, c+dc);
    V  = lat.V;
    V000 = V(L(0,0,0),:);   V100 = V(L(1,0,0),:);
    V010 = V(L(0,1,0),:);   V110 = V(L(1,1,0),:);
    V001 = V(L(0,0,1),:);   V101 = V(L(1,0,1),:);
    V011 = V(L(0,1,1),:);   V111 = V(L(1,1,1),:);
    lp = @(A,Bb,t) A + t .* (Bb - A);
    B00 = lp(V000, V100, tx);   B01 = lp(V001, V101, tx);
    B10 = lp(V010, V110, tx);   B11 = lp(V011, V111, tx);
    C0  = lp(B00,  B10,  ty);   C1  = lp(B01,  B11,  ty);
    Vq  = lp(C0,   C1,   tz);
    Bq  = reshape(Vq, [], 3, lat.N_I);
    loc = struct('i0',i0, 'j0',j0, 'k0',k0, 'tx',tx, 'ty',ty, 'tz',tz);
end
