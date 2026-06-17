# SpaceClaim IronPython script: batch convert IGES → STEP
# Usage: SpaceClaim.exe -Script convert_iges_to_step.py -Headless -ExitAfterScript
import os
import sys
from SpaceClaim.Api.V261 import *

src_dir = r"C:\Users\pmero\Documents\Lab406\FEM_sim\magnetic_sim\hung\IGES_converted"
dst_dir = r"C:\Users\pmero\Documents\Lab406\magnetic-tweezer-cad\hung\step_for_measure"

files = [f for f in os.listdir(src_dir) if f.lower().endswith('.iges')]
print("Found %d IGES files" % len(files))

for f in files:
    src_path = os.path.join(src_dir, f)
    step_name = f.replace('.iges', '.step').replace('.IGES', '.step')
    dst_path = os.path.join(dst_dir, step_name)

    print("Converting: %s" % f)
    try:
        options = ImportOptions.Create()
        DocumentInsert.Execute(src_path, options)

        export_options = ExportOptions.Create()
        DocumentSave.Execute(dst_path, export_options)

        Window.ActiveWindow.Close()
        print("  -> OK: %s" % step_name)
    except Exception as e:
        print("  -> FAILED: %s" % str(e))

print("Done. %d files processed." % len(files))
