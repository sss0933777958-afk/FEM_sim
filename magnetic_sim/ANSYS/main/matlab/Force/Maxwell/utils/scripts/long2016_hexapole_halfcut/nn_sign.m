%% nn_sign.m -- sign-align the NN labels: I_opt_hat and -I_opt_hat give the same force, so fmincon returns
%  either one at random; this picks one sign per (direction, position) so the labels vary continuously.
%  The whole 6-vector is multiplied by +-1 (never single components), so RHO, ERR and the force are unchanged.
%  Graph: one node per (phi, theta, position); edges join neighbours along phi, along theta (360 -> 0 wraps),
%  and along the 30 um position lattice (x_a, y_a, z_a).  Edge weight 1 - |cos|, cos = angle between the two
%  6-vectors.  The signs are propagated from the root along the maximum-|cos| spanning tree:
%      s(child) = s(parent) * sign(cos(child, parent))
%  Root: centre position p = 0, phi = 0, theta = 0, sign chosen so that its P1 current is positive.
%  Edges left with cos < 0 after the alignment cannot be fixed by a sign (branch switch of I_opt_hat itself).
%
%  Input:  utils/data/long2016_hexapole_halfcut/nn_label.mat (nn_label.m), nn_grid.mat (grid step H)
%  Output: utils/data/long2016_hexapole_halfcut/nn_label_sgn.mat
%    every variable of nn_label.mat, with IOPT replaced by the sign-aligned labels, plus
%    SGN     NPH x NTH x NP   int8 +-1 applied to the nn_label.mat IOPT
%    CDIR, CPOS             cos over the direction edges / the position edges after the alignment
here = fileparts(mfilename('fullpath'));  CALROOT = fileparts(fileparts(fileparts(here)));   % .../Force/Maxwell
DAT = fullfile(CALROOT, 'utils', 'data', 'long2016_hexapole_halfcut');
S   = load(fullfile(DAT, 'nn_label.mat'));
G   = load(fullfile(DAT, 'nn_grid.mat'), 'H');
[~, NPH, NTH, NP] = size(S.IOPT);
NN  = NPH * NTH * NP;
V   = reshape(S.IOPT, 6, NN);                                  % node n = sub2ind([NPH NTH NP], a, b, i)
V   = V ./ vecnorm(V);
id  = reshape(1:NN, NPH, NTH, NP);

% ---- edges ----
e1 = id(1:end-1, :, :);  f1 = id(2:end, :, :);                % along phi
e2 = id(:, [1:end-1 end], :);  f2 = id(:, [2:end 1], :);       % along theta, theta(end) -> theta(1) wrap
K  = round(S.P_act / G.H);  E = eye(3);  ep = [];  fp = [];     % position lattice neighbours
for k = 1:3
    [hit, j] = ismember(K + E(k,:), K, 'rows');
    ep = [ep; find(hit)];  fp = [fp; j(hit)];                   %#ok<AGROW>
end
e3 = id(:, :, ep);  f3 = id(:, :, fp);
s = [e1(:); e2(:); e3(:)];  t = [f1(:); f2(:); f3(:)];
ND = numel(e1) + numel(e2);                                    % first ND edges: direction edges
c  = sum(V(:, s) .* V(:, t), 1).';                             % cos of every edge (before alignment)

% ---- maximum-|cos| spanning tree, signs from the root ----
i0   = find(all(S.P_act == 0, 2));  a0 = find(S.ph == 0);  b0 = find(S.th == 0);
root = id(a0, b0, i0);
[~, pred] = minspantree(graph(s, t, 1 - abs(c) + eps, NN), 'Root', root);
ord  = bfsearch(graph(find(pred > 0), pred(pred > 0), ones(nnz(pred > 0), 1), NN), root);
sg   = zeros(NN, 1);  sg(root) = sign(S.IOPT(1, a0, b0, i0));
for n = ord(2:end).'
    p = pred(n);
    sg(n) = sg(p) * sign(V(:, n).' * V(:, p));
end
assert(all(sg ~= 0), 'nn_sign:unreached', '%d nodes not reached', nnz(sg == 0));

% ---- apply and check ----
SGN  = int8(reshape(sg, NPH, NTH, NP));
IOPT = S.IOPT .* reshape(double(SGN), 1, NPH, NTH, NP);
ca   = c .* sg(s) .* sg(t);                                    % cos after alignment
CDIR = ca(1:ND);  CPOS = ca(ND+1:end);
S.IOPT = IOPT;  S.SGN = SGN;  S.CDIR = CDIR;  S.CPOS = CPOS;
save(fullfile(DAT, 'nn_label_sgn.mat'), '-struct', 'S', '-v7.3');
fprintf('flipped %.1f %% of %d labels -> %s\n', 100*mean(sg < 0), NN, fullfile(DAT, 'nn_label_sgn.mat'));
fprintf('direction edges %d: cos < 0 before %.1f %%, after %.3f %% | cos < 0.5 after %.3f %%\n', ...
        ND, 100*mean(c(1:ND) < 0), 100*mean(CDIR < 0), 100*mean(CDIR < 0.5));
fprintf('position  edges %d: cos < 0 before %.1f %%, after %.3f %% | cos < 0.5 after %.3f %%\n', ...
        numel(CPOS), 100*mean(c(ND+1:end) < 0), 100*mean(CPOS < 0), 100*mean(CPOS < 0.5));
