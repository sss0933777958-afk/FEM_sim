# SpaceClaim (IronPython) script: hung_hexapole Full_Assembly.STEP -> Parasolid (.x_t)
# 供 MAPDL ~PARAIN 匯入 (ANSYS 對 Parasolid 遠比 IGES 穩健)
# 執行: SpaceClaim.exe /Headless=True /Splash=False /ExitAfterScript=True /RunScript="<this>"
step = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\CAD_model\hung_hexapole\STEP\Full_Assembly.STEP"
out  = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\apdl\hung_hexapole\geom\import\hung_hexapole_full.x_t"
DocumentOpen.Execute(step)
DocumentSave.Execute(out)
