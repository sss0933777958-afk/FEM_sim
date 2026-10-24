function S = rho_solve(Q, o)
% rho_solve -- minimum-norm unit current over the whole direction grid at ONE position (actuator frame).
% Shared solver of rho_map.m (one position) and cap_map.m (positions along rays); the method is documented in the
% rho_map.m header:
%   stage 1  fmincon (sqp) multistart per direction + Lagrange PSD certificate (CERT)
%   stage 2  (o.GRB) Gurobi on every direction without the certificate, local fmincon polish, keep the smaller
%            ||I_hat||^2 only if the polish succeeded (PLSH); GLOB = CERT | Gurobi OPTIMAL
% Q   6x6x3 force quadratic forms Q_k at the position
% o   struct, all optional:
%     DA     direction grid step [deg] (5)          NS    random starts per direction (20)
%     GRB    run stage 2 (true)                      GTHR  Gurobi Threads, 0 = Gurobi default (0)
%     TLIM   Gurobi TimeLimit per direction [s] (60) SEED  rng seed of the random starts (0)
%     I_INIT first start of the first direction ([0 0 0 0 0 1]')
%     X0G    6 x NPH x NTH extra start per direction, e.g. the solution at a neighbouring position ([])
%     S1     stage-1 struct (fields IOPT FHAT LAM RHO ERR MINEIG NFEAS) -> stage 1 is skipped ([])
%     RETRY  a direction with no feasible start gets further batches of NS random starts (random scale 0.05..1)
%            until one is feasible (false; at most 200 batches, then it stays NaN with a warning)
%     GMAX   stage 2 on at most GMAX of the uncertified directions, a random subset (seed SEED) (Inf = all)
%     DOFF   direction grid offset [deg] (0): ph = -90+DOFF:DA:90-DOFF, th = DOFF:DA:360+DOFF (DOFF = DA/2: the
%            directions halfway between the grid ones, no pole)
% S   struct: ph th, per-direction NPH x NTH arrays (IOPT FHAT LAM RHO ERR MINEIG NFEAS CERT GLOB KKT and the stage-2
%     fields GRBD GSTAT GOBJ GBND GGAP DOBJ PLSH), stage-1 copies (*0), cap, iso, gp, t1 / t2 (stage times [s])
o = fill_(o, struct('DA',5, 'NS',20, 'GRB',true, 'GTHR',0, 'TLIM',60, 'SEED',0, 'I_INIT',[0 0 0 0 0 1].', ...
                    'X0G',[], 'S1',[], 'RETRY',false, 'GMAX',Inf, 'DOFF',0));
TOLC = 1e-6;  TOLE = 1e-8;                                  % constraint / PSD tolerances
opt = optimoptions('fmincon', 'Algorithm','sqp', 'SpecifyObjectiveGradient',true, ...
                   'SpecifyConstraintGradient',true, 'OptimalityTolerance',1e-10, 'ConstraintTolerance',1e-12, ...
                   'StepTolerance',1e-14, 'MaxIterations',1000, 'Display','off');
ph = -90+o.DOFF:o.DA:90-o.DOFF;  th = o.DOFF:o.DA:360+o.DOFF;  NPH = numel(ph);  NTH = numel(th);
t0 = tic;

% ---- stage 1: fmincon multistart + PSD certificate ----
if ~isempty(o.S1)
    IOPT = o.S1.IOPT;  FHAT = o.S1.FHAT;  LAM = o.S1.LAM;  RHO = o.S1.RHO;  ERR = o.S1.ERR;
    MINEIG = o.S1.MINEIG;  NFEAS = o.S1.NFEAS;
else
    rng(o.SEED);
    IOPT = zeros(6, NPH, NTH);  FHAT = zeros(3, NPH, NTH);  LAM = zeros(3, NPH, NTH);  RHO = nan(NPH, NTH);
    ERR  = nan(NPH, NTH);  MINEIG = nan(NPH, NTH);  NFEAS = zeros(NPH, NTH);
    x0 = o.I_INIT;
    for a = 1:NPH
        for b = 1:NTH
            if abs(ph(a)) == 90 && b > 1                     % pole: one direction
                IOPT(:,a,b) = IOPT(:,a,1);  FHAT(:,a,b) = FHAT(:,a,1);  LAM(:,a,b) = LAM(:,a,1);  RHO(a,b) = RHO(a,1);
                ERR(a,b) = ERR(a,1);  MINEIG(a,b) = MINEIG(a,1);  NFEAS(a,b) = NFEAS(a,1);  continue;
            end
            r  = [cosd(ph(a))*cosd(th(b)); cosd(ph(a))*sind(th(b)); sind(ph(a))];
            R0 = randn(6, o.NS);  X0 = [x0, 0.35 * R0 ./ vecnorm(R0)];   % chain start + NS random starts
            if ~isempty(o.X0G) && any(o.X0G(:,a,b)), X0 = [X0, o.X0G(:,a,b)]; end   %#ok<AGROW> neighbour start
            best = inf;  nb = 0;
            while true
                for s = 1:size(X0, 2)
                    [I, ~, ef, ~, lm] = fmincon(@(x) obj_(x), X0(:,s), [],[],[],[],[],[], @(x) con_(x, Q, r), opt);
                    f = [I.'*Q(:,:,1)*I; I.'*Q(:,:,2)*I; I.'*Q(:,:,3)*I];
                    if ef <= 0 || norm(f - r) >= TOLC, continue; end
                    NFEAS(a,b) = NFEAS(a,b) + 1;
                    if I.'*I < best
                        best = I.'*I;  lam = -lm.eqnonlin;
                        IOPT(:,a,b) = I;  FHAT(:,a,b) = f;  LAM(:,a,b) = lam;  RHO(a,b) = 1/best;  ERR(a,b) = norm(f - r);
                        MINEIG(a,b) = min(eig(eye(6) - (lam(1)*Q(:,:,1) + lam(2)*Q(:,:,2) + lam(3)*Q(:,:,3))));
                    end
                end
                if isfinite(best) || ~o.RETRY, break; end
                nb = nb + 1;                                % no feasible start yet: another batch
                if nb > 200
                    warning('rho_solve:infeasible', 'phi %g theta %g: no feasible start after %d batches', ph(a), th(b), nb);
                    break
                end
                R0 = randn(6, o.NS);  X0 = (0.05 + 0.95*rand(1, o.NS)) .* R0 ./ vecnorm(R0);
            end
            if nb > 0 && isfinite(best)
                fprintf('    rho_solve: phi %g theta %g feasible after %d extra batch(es) of %d starts\n', ph(a), th(b), nb, o.NS);
            end
            if isfinite(best), x0 = IOPT(:,a,b); end       % best solution seeds the next grid point
        end
    end
end
t1 = toc(t0);
IOPT0 = IOPT;  FHAT0 = FHAT;  LAM0 = LAM;  RHO0 = RHO;  ERR0 = ERR;  MINEIG0 = MINEIG;   % stage-1 (fmincon) copies

% ---- stage 2: Gurobi on the directions without the PSD certificate ----
GRBD = false(NPH, NTH);  GSTAT = strings(NPH, NTH);  GOBJ = nan(NPH, NTH);  GBND = nan(NPH, NTH);
GGAP = nan(NPH, NTH);  DOBJ = nan(NPH, NTH);  PLSH = false(NPH, NTH);
gp = struct('NonConvex',2, 'MIPGap',1e-9, 'MIPGapAbs',1e-12, 'FeasibilityTol',1e-9, 'OptimalityTol',1e-9, ...
            'NumericFocus',2, 'TimeLimit',o.TLIM, 'Seed',0, 'OutputFlag',0, 'Threads',o.GTHR);
if o.GRB
    V = MINEIG0 < -TOLE;  V(abs(ph) == 90, 2:end) = false;  % theta = 360 column solved too (seam check)
    [ia, ib] = find(V);
    if numel(ia) > o.GMAX                                   % random subset of the uncertified directions
        ks = randperm(RandStream('mt19937ar', 'Seed', o.SEED), numel(ia), o.GMAX);  ia = ia(ks);  ib = ib(ks);
    end
    for n = 1:numel(ia)
        a = ia(n);  b = ib(n);
        r  = [cosd(ph(a))*cosd(th(b)); cosd(ph(a))*sind(th(b)); sind(ph(a))];
        I0 = IOPT0(:,a,b);  B = norm(I0);
        m = struct('Q',speye(6), 'obj',zeros(6,1), 'A',sparse(0,6), 'rhs',zeros(0,1), 'sense','', ...
                   'lb',-B*ones(6,1), 'ub',B*ones(6,1), 'start',I0);
        for k = 1:3
            m.quadcon(k) = struct('Qc',sparse(Q(:,:,k)), 'q',zeros(6,1), 'rhs',r(k), 'sense','=');
        end
        m.quadcon(4) = struct('Qc',speye(6), 'q',zeros(6,1), 'rhs',B^2, 'sense','<');   % ball bound
        res = gurobi(m, gp);
        GRBD(a,b) = true;  GSTAT(a,b) = res.status;
        if ~isfield(res, 'x'), continue; end
        GOBJ(a,b) = res.objval;  GBND(a,b) = res.objbound;  GGAP(a,b) = res.mipgap;  DOBJ(a,b) = I0.'*I0 - res.objval;
        [I, ~, ef, ~, lm] = fmincon(@(x) obj_(x), res.x, [],[],[],[],[],[], @(x) con_(x, Q, r), opt);   % polish
        f = [I.'*Q(:,:,1)*I; I.'*Q(:,:,2)*I; I.'*Q(:,:,3)*I];
        PLSH(a,b) = ef > 0 && norm(f - r) < TOLC;          % polish converged and feasible
        lam = -lm.eqnonlin;
        % keep the smaller ||I_hat||^2; an unpolished Gurobi point is never taken: it violates the constraints by up
        % to FeasibilityTol, which alone can make ||I||^2 look ~1e-9 smaller than the exact fmincon solution
        if PLSH(a,b) && I.'*I < I0.'*I0
            IOPT(:,a,b) = I;  FHAT(:,a,b) = f;  LAM(:,a,b) = lam;  RHO(a,b) = 1/(I.'*I);  ERR(a,b) = norm(f - r);
            MINEIG(a,b) = min(eig(eye(6) - (lam(1)*Q(:,:,1) + lam(2)*Q(:,:,2) + lam(3)*Q(:,:,3))));
        end
    end
    for a = find(abs(ph) == 90)                             % pole: copy the single direction
        IOPT(:,a,:) = repmat(IOPT(:,a,1), 1, 1, NTH);  FHAT(:,a,:) = repmat(FHAT(:,a,1), 1, 1, NTH);
        LAM(:,a,:)  = repmat(LAM(:,a,1), 1, 1, NTH);    RHO(a,:) = RHO(a,1);  ERR(a,:) = ERR(a,1);
        MINEIG(a,:) = MINEIG(a,1);  GRBD(a,:) = GRBD(a,1);  GSTAT(a,:) = GSTAT(a,1);  GOBJ(a,:) = GOBJ(a,1);
        GBND(a,:) = GBND(a,1);  GGAP(a,:) = GGAP(a,1);  DOBJ(a,:) = DOBJ(a,1);  PLSH(a,:) = PLSH(a,1);
    end
end
KKT = nan(NPH, NTH);                                        % stationarity residual ||I - sum lambda_k Q_k I||
for a = 1:NPH
    for b = 1:NTH
        I = IOPT(:,a,b);  lam = LAM(:,a,b);
        KKT(a,b) = norm(I - (lam(1)*Q(:,:,1) + lam(2)*Q(:,:,2) + lam(3)*Q(:,:,3))*I);
    end
end

% capacity (trapezoid in theta and phi, solid-angle weight cos(phi)) and isotropy
CERT = MINEIG >= -TOLE;                                    % PSD certificate (sufficient)
GLOB = CERT | GSTAT == "OPTIMAL";                          % global minimum: PSD or Gurobi (gap within its tolerances)
cap  = (trapz(deg2rad(ph), trapz(deg2rad(th), RHO.^3, 2) .* cosd(ph(:))) / (4*pi))^(1/3);
iso  = min(RHO(:)) / max(RHO(:));
t2 = toc(t0) - t1;

S = struct('ph',ph, 'th',th, 'IOPT',IOPT, 'FHAT',FHAT, 'LAM',LAM, 'RHO',RHO, 'ERR',ERR, 'MINEIG',MINEIG, ...
           'NFEAS',NFEAS, 'CERT',CERT, 'GLOB',GLOB, 'KKT',KKT, 'cap',cap, 'iso',iso, ...
           'IOPT0',IOPT0, 'FHAT0',FHAT0, 'LAM0',LAM0, 'RHO0',RHO0, 'ERR0',ERR0, 'MINEIG0',MINEIG0, ...
           'GRBD',GRBD, 'GSTAT',GSTAT, 'GOBJ',GOBJ, 'GBND',GBND, 'GGAP',GGAP, 'DOBJ',DOBJ, 'PLSH',PLSH, ...
           'gp',gp, 't1',t1, 't2',t2);
end

% ---------------------------------------------------------------------------------------
function o = fill_(o, d)
% defaults d for the fields missing in o
    if isempty(o), o = struct(); end
    for f = fieldnames(d).'
        if ~isfield(o, f{1}), o.(f{1}) = d.(f{1}); end
    end
end

function [v, g] = obj_(x)
% ||I_hat||^2 and its gradient
    v = x.'*x;  g = 2*x;
end

function [c, ceq, gc, gceq] = con_(x, Q, r)
% equality constraints I_hat' Q_k I_hat - r_k = 0 (k = 1..3) and their gradients (6 x 3)
    c = [];  gc = [];
    ceq  = [x.'*Q(:,:,1)*x - r(1); x.'*Q(:,:,2)*x - r(2); x.'*Q(:,:,3)*x - r(3)];
    gceq = 2*[Q(:,:,1)*x, Q(:,:,2)*x, Q(:,:,3)*x];
end
