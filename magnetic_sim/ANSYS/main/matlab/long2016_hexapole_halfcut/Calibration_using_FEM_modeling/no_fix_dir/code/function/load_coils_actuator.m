function D = load_coils_actuator(model, cnst, apdl_to_paper_idx, dataset)
%LOAD_COILS_ACTUATOR  Load the 6-coil FEM field once and rotate into the actuator frame.
%   D = LOAD_COILS_ACTUATOR(model, cnst, apdl_to_paper_idx, dataset)
%   Reads coilN/standard <dataset> via ansys_path, removes iron nodes, shifts z to the
%   WP frame, NEGATES B (all-source convention), and rotates points + fields by the
%   actuator rotation R_act = [uhat vhat what]' (uhat=P1, vhat=P3, what=P5 tip dirs).
%   Returns struct D with:
%     .Pa      Nair x 3       sample points (actuator frame, shared across coils)
%     .r2      Nair x 1       |p|^2 (rotation-invariant) for ball selection
%     .Ba      Nair x 3 x N_I all-source B per simulation (actuator frame)
%     .R_act   3 x 3          measure->actuator rotation
%     .Pc_base 3 x 6          ideal charge lattice [+u -u +v -v +w -w]
%     .F       6 x N_I        current matrix (permutation, rank 6)
%   Requires ansys_path + import_ansys_data + filter_iron_nodes on the path.
    N_I = 6;
    tip   = [cnst.pole_tip_x; cnst.pole_tip_y; cnst.pole_tip_z_wp];   % 3x6 (measure)
    dhat  = tip ./ vecnorm(tip);
    R_act = [dhat(:,1), dhat(:,3), dhat(:,5)].';                      % uhat=P1, vhat=P3, what=P5
    Pc_base = [ 1 -1  0  0  0  0;
                0  0  1 -1  0  0;
                0  0  0  0  1 -1];
    assert(abs(det(R_act)-1) < 1e-9, 'R_act must be a proper rotation');
    assert(max(abs(R_act*dhat - Pc_base), [], 'all') < 1e-9, 'R_act*dhat must equal Pc_base');

    F = zeros(6, N_I);
    for j = 1:N_I, F(apdl_to_paper_idx(j), j) = 1; end

    d1   = import_ansys_data(ansys_path(model,'data','coil1','standard'), dataset, 'coil1');
    air1 = filter_iron_nodes(d1.x,d1.y,d1.z,cnst,struct('visualize',false));
    zwp1 = d1.z - cnst.SPH_OFST;
    P_meas = [d1.x(air1), d1.y(air1), zwp1(air1)];
    Pa  = (R_act * P_meas.').';
    r2  = sum(P_meas.^2, 2);

    Nair = size(P_meas,1);
    Ba = zeros(Nair, 3, N_I);
    for k = 1:N_I
        if k == 1, dk = d1; airk = air1;
        else
            cn = sprintf('coil%d', k);
            dk = import_ansys_data(ansys_path(model,'data',cn,'standard'), dataset, cn);
            airk = filter_iron_nodes(dk.x,dk.y,dk.z,cnst,struct('visualize',false));
        end
        Bk = -[dk.bx(airk), dk.by(airk), dk.bz(airk)];               % all-source
        Ba(:,:,k) = (R_act * Bk.').';
    end

    D = struct('Pa',Pa,'r2',r2,'Ba',Ba,'R_act',R_act,'Pc_base',Pc_base,'F',F);
end
