"""從 MT_Sim_P1_gap200.txt 生成 P2-P6 (只差 CURR_ARRAY: coilN=1 餘 0 + 輸出路徑 coilN)。"""
import re,pathlib
d=pathlib.Path(r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\apdl\hung_hexapole\sim\gap200um_mueq")
tmpl=(d/"MT_Sim_P1_gap200.txt").read_text(encoding="utf-8")
for N in range(2,7):
    s=tmpl
    # 輸出/標題/註解 coil1 -> coilN (RESUME 的 mesh_graded 不含 coil1 字樣, 安全全換)
    s=s.replace("coil1","coil%d"%N)
    # R block: 第 N 顆 *1, 其餘 *0
    for k in range(1,7):
        mult=1 if k==N else 0
        s=re.sub(r"R,%d,1,TURNS\*[01],COIL_THK,COIL_H"%k,
                 "R,%d,1,TURNS*%d,COIL_THK,COIL_H"%(k,mult),s)
    (d/("MT_Sim_P%d_gap200.txt"%N)).write_text(s,encoding="utf-8")
    print("wrote MT_Sim_P%d_gap200.txt (coil%d=1)"%(N,N))
