@echo off
"C:\Program Files\ANSYS2025R2\v252\scdm\SpaceClaim.exe" -Script "C:\Users\pmero\Documents\Lab406\FEM_sim\magnetic_sim\hung\scripts\convert_iges_to_step.py" -Headless -ExitAfterScript
echo Done. Exit code: %ERRORLEVEL%
