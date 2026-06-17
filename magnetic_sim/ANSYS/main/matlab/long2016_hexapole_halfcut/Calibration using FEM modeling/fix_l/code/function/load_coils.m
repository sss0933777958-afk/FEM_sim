function C = load_coils(results_root, cnst, apdl_to_paper_idx)
%LOAD_COILS  Load the 6 single-coil FEM B-fields (WP frame, air nodes, all-source sign).
%   C = LOAD_COILS(results_root, cnst, apdl_to_paper_idx)
%   Reads coil1..6 'wp' dataset under <results_root>\coilN\standard, removes iron
%   nodes (filter_iron_nodes), shifts z into the WP frame (z - SPH_OFST), and NEGATES
%   B (source convention: each excited tip radiates outward -> positive charge).
%   Returns struct array C(k):
%     .P  Nx3 [m]  air-node positions (WP frame)
%     .Bn Nx3 [T]  FEM B at those nodes, sign-negated
%     .pj scalar   paper pole index excited by coil k
%   Requires import_ansys_data + filter_iron_nodes on the path (hexapole-long2016\analysis).
    C = struct('P',{},'Bn',{},'pj',{});
    for k = 1:6
        cn  = sprintf('coil%d', k);
        d   = import_ansys_data(fullfile(results_root, cn, 'standard'), 'wp', cn);
        air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize',false));
        zwp = d.z - cnst.SPH_OFST;
        C(k).P  =  [d.x(air),  d.y(air),  zwp(air)];
        C(k).Bn = -[d.bx(air), d.by(air), d.bz(air)];
        C(k).pj =  apdl_to_paper_idx(k);
    end
end
