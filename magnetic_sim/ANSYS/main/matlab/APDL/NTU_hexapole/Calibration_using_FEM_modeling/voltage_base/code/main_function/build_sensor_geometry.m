function [sensor_pos, sensor_n] = build_sensor_geometry(cnst)
%BUILD_SENSOR_GEOMETRY  NTU：6 顆 S2 Hall sensor（每根極板上方）中心 + 法線（world 框）。
%   S2 = 板頂上方 0.41mm、pole-local (13.15,−0.99) 旋方位角、n+ = +z（離板朝上）。
%   直接取自 NTU mt_constants（cnst.sensor_pos/sensor_n，已在 mt_constants 依 6 極方位算好）。
%   輸出：sensor_pos 3×6 [m]（world/ANSYS 框，無需再 shift）、sensor_n 3×6（單位法線）。
    sensor_pos = cnst.sensor_pos;      % 3×6 [m]，world
    sensor_n   = cnst.sensor_n;        % 3×6，+z
end
