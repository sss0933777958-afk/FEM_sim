function Pc = make_Pc(e17, Pc_base)
%MAKE_PC  Assemble the 3x6 charge-location pattern Pc = Pc_base + E(e_hat).
%   Pc = MAKE_PC(e17, Pc_base)
%   18-parameter bias model: 17 free offsets e17 fill E, and the 6th-pole z entry is
%   constrained: e6z = e1x - e2x + e3y - e4y + e5z.  (no_fix_l.pdf step 2c.)
    E = zeros(3, 6);
    E(:,1) = e17(1:3);    E(:,2) = e17(4:6);
    E(:,3) = e17(7:9);    E(:,4) = e17(10:12);
    E(:,5) = e17(13:15);
    E(1,6) = e17(16);     E(2,6) = e17(17);
    E(3,6) = e17(1) - e17(4) + e17(8) - e17(11) + e17(15);   % e6z constraint
    Pc = Pc_base + E;
end
