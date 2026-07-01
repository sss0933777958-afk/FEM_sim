function C = load_coils(results_root, cnst, apdl_to_paper_idx, variant)
%LOAD_COILS  Load the 6 single-coil FEM B-fields (WP frame, air nodes, all-source sign).
%   C = LOAD_COILS(results_root, cnst, apdl_to_paper_idx, variant)
%   Reads coil1..6 'wp' dataset under <results_root>\coilN\<variant>, removes iron
%   nodes (filter_iron_nodes), shifts z into the WP frame (z - SPH_OFST), and makes the
%   field ALL-SOURCE by literal FLIP-SINK: negate only the lower-pole (sink) excitations
%   P1/P3/P6; upper poles P2/P4/P5 are already source (raw) and are kept. Every excited
%   tip then radiates outward -> positive self-charge (diagonal all positive).
%   variant (optional, default 'standard'): FEM variant subfolder, e.g. 'gap200um_mueq'.
%   Returns struct array C(k):
%     .P  Nx3 [m]  air-node positions (WP frame)
%     .Bn Nx3 [mT] FEM B at those nodes (ANSYS Tesla ×1e3 → 原生 mT), all-source (flip-sink: lower P1/P3/P6 negated)
%     .pj scalar   paper pole index excited by coil k
%   Requires import_ansys_data + filter_iron_nodes on the path (hexapole-long2016\analysis).
    if nargin < 4 || isempty(variant), variant = 'standard'; end     % [MODIFIED] optional variant (back-compat)
    C = struct('P',{},'Bn',{},'pj',{});
    for k = 1:6
        cn  = sprintf('coil%d', k);
        d   = import_ansys_data(fullfile(results_root, cn, variant), 'wp', cn);  % [MODIFIED] variant subfolder
        air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize',false));
        zwp = d.z - cnst.SPH_OFST;
        C(k).P  =  [d.x(air),  d.y(air),  zwp(air)];
        sgn = 1; if ismember(apdl_to_paper_idx(k), [1 3 6]), sgn = -1; end   % [MODIFIED] flip-sink：只翻下極 sink
        C(k).Bn = sgn*1e3*[d.bx(air), d.by(air), d.bz(air)];   % ANSYS Tesla → ×1e3 原生 mT（Unit Sheet）
        C(k).pj =  apdl_to_paper_idx(k);
    end
end
