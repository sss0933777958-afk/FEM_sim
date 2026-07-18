function B = ntu_pole_profile()
%NTU_POLE_PROFILE  NTU 磁極板 2D 輪廓（pole-local，mm）：9 段環離散成 Nx2 polygon。
%   數值來源：apdl/NTU_hexapole/geom/export/MT_Geom_FullAssembly.txt 的 9 段 K/LARC
%   （A~I 九點 + 尖端 r0.005 圓角 ×2 + r10 凹弧 ×2 + 遠端 r5 凸弧拆兩段）。
%   世界座標請自行旋轉方位角：Pw = (R(θ) * B.').'
%
%   [ADDED 2026-07-17] 由 plot/plot_poles_sensors_3d.m 與 plot/plot_p1_sensing_disk.m
%   的重複 local function `build_pole_profile` 抽出共用（原本兩支各一份、無法跨檔呼叫）。
%   行為與原版逐字相同（弧離散點數：尖端 6 / 凹弧 10 / r5 弧 80）。
    A=[0.357 0.005]; Bp=[0.353 0.000]; C=[0.357 -0.005];
    D=[4.509 -0.523]; E=[6.175 -0.877]; F=[19.048 -4.784];
    G=[19.048 4.784]; H=[6.175 0.877]; I=[4.509 0.523];
    Ctip=[0.358 0]; C10lo=[3.271 -10.446]; C5=[20.5 0]; C10up=[3.271 10.446];
    seg = [ arc(Ctip,0.005,A,Bp,6); arc(Ctip,0.005,Bp,C,6); C; D; ...
            arc(C10lo,10,D,E,10); E; F; arc(C5,5,F,G,80); G; H; ...
            arc(C10up,10,H,I,10); I; A ];
    keep=[true; any(abs(diff(seg))>1e-9,2)];  B=seg(keep,:);
end

function P = arc(C, r, P1, P2, n)
    a1=atan2(P1(2)-C(2),P1(1)-C(1)); a2=atan2(P2(2)-C(2),P2(1)-C(1));
    t=linspace(a1,a2,n).';  P=[C(1)+r*cos(t), C(2)+r*sin(t)];
end
