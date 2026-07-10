function C = load_coils(results_root, cnst, apdl_to_paper_idx, variant)
%LOAD_COILS  Load the 6 single-coil FEM B-fields (ACTUATOR frame, air nodes, all-source sign).
%   C = LOAD_COILS(results_root, cnst, apdl_to_paper_idx, variant)
%   Reads coil1..6 'wp' dataset under <results_root>\coilN\<variant>, removes iron
%   nodes (filter_iron_nodes), shifts z into the WP frame (z - SPH_OFST), makes the
%   field ALL-SOURCE by literal FLIP-SINK (negate only lower-pole sink excitations
%   P1/P3/P6; upper P2/P4/P5 kept), then ROTATES points + fields into the actuator frame
%   R_act = [d̂_P1 d̂_P3 d̂_P5]' — so the fit residual (B_model − B_FEM) is formed in the
%   actuator frame, consistent with the other 3 calibration packages. Frame rotation is
%   norm-preserving → fit results (ell, ^Bg_I, K̄, err) are unchanged vs the old WP frame.
%   In the actuator frame the on-axis charges sit at ell·Pc_base (Pc_base = R_act·d̂).
%   variant (optional, default 'standard'): FEM variant subfolder, e.g. 'gap_200um'.
%   Returns struct array C(k):
%     .P  Nx3 [m]  air-node positions (ACTUATOR frame)
%     .Bn Nx3 [mT] FEM B at those nodes (ANSYS Tesla ×1e3 → mT), all-source, ACTUATOR frame
%     .pj scalar   paper pole index excited by coil k
%   Requires import_ansys_data + filter_iron_nodes on the path (hexapole-long2016\analysis).
    if nargin < 4 || isempty(variant), variant = 'standard'; end     % [MODIFIED] optional variant (back-compat)
    tip   = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];   % [ADDED] 3x6 pole-tip dirs (WP frame)
    dhat  = tip ./ vecnorm(tip);
    R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';                      % [ADDED] measure(WP)→actuator rotation
    C = struct('P',{},'Bn',{},'pj',{});
    for k = 1:6
        cn  = sprintf('coil%d', k);
        d   = import_ansys_data(fullfile(results_root, cn, variant), 'wp', cn);  % [MODIFIED] variant subfolder
        air = filter_iron_nodes(d.x, d.y, d.z, cnst, struct('visualize',false));
        zwp = d.z - cnst.SPH_OFST;
        sgn = 1; if ismember(apdl_to_paper_idx(k), [1 3 6]), sgn = -1; end   % flip-sink：只翻下極 sink
        Pm  = [d.x(air),  d.y(air),  zwp(air)];                       % WP-frame air nodes
        Bm  = sgn*1e3*[d.bx(air), d.by(air), d.bz(air)];              % all-source, mT (ANSYS Tesla ×1e3)
        C(k).P  = ( R_act * Pm.' ).';                                % [MODIFIED] → actuator frame
        C(k).Bn = ( R_act * Bm.' ).';                                % [MODIFIED] → actuator frame (rotation keeps all-source sign)
        C(k).pj =  apdl_to_paper_idx(k);
    end
end
