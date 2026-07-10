"""
hung_hexapole air-domain builder (OCP / OpenCascade robust booleans).
MAPDL cannot boolean imported NURBS -> build the air domain here instead:
  read Full_Assembly.STEP -> classify solids by bbox -> keep 49 steel (poles+supports+yoke),
  drop 48 non-magnetic (24 guide posts 4x4 + 24 blocks 10x10x14) -> FUSE steel into one body
  -> air = sphere(R=150mm) CUT steel -> export STEP (steel + air-with-cavities).
Downstream: STEP -> SpaceClaim .x_t -> ac4para .anf -> MAPDL /INPUT -> mesh + NUMMRG,NODE (conformal).
Units: STEP is mm (OCC reads raw mm). Output STEP in mm.
"""
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_StepModelType
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.TopoDS import TopoDS, TopoDS_Compound
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut
from OCP.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCP.gp import gp_Pnt
from OCP.BRep import BRep_Builder
from OCP.TopTools import TopTools_ListOfShape
from OCP.Interface import Interface_Static
Interface_Static.SetCVal_s("write.step.unit", "MM")  # declare MM so downstream (Mechanical) reads mm not m

STEP_IN = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main\CAD_model\hung_hexapole\STEP\Full_Assembly.STEP"
STEP_OUT= r"C:\Users\Kuo\AppData\Local\Temp\claude\G--my-workspace\a41529f7-cb06-482b-bf27-63298351ac0c\scratchpad\hung_air.step"
R_SPHERE = 150.0  # mm

def solids(shape):
    out=[]; e=TopExp_Explorer(shape,TopAbs_SOLID)
    while e.More(): out.append(TopoDS.Solid_s(e.Current())); e.Next()
    return out

def bbox(s):
    # AddOptimal = tight box from actual surfaces (plain Add uses NURBS control pts -> inflates thin posts)
    b=Bnd_Box(); BRepBndLib.AddOptimal_s(s,b,True,False)
    xmn,ymn,zmn,xmx,ymx,zmx=b.Get()
    return xmn,ymn,zmn,xmx,ymx,zmx

def is_nonmagnetic(s):
    xmn,ymn,zmn,xmx,ymx,zmx=bbox(s)
    dx,dy,dz=xmx-xmn,ymx-ymn,zmx-zmn; foot=max(dx,dy); cz=0.5*(zmn+zmx)
    if foot<=6 and dz>=12: return True                         # guide post 4x4 tall/short
    if 8<=foot<=12 and 12<=dz<=16 and (5<cz<10 or 35<cz<41): return True  # 10x10x14 block
    return False

r=STEPControl_Reader(); r.ReadFile(STEP_IN); r.TransferRoots()
allsol=solids(r.OneShape())
print("total solids:",len(allsol))
steel=[s for s in allsol if not is_nonmagnetic(s)]
print("steel (keep):",len(steel)," non-magnetic (drop):",len(allsol)-len(steel))

# fuse all steel into one body (robust OCC boolean)
fused=steel[0]
for i,s in enumerate(steel[1:],2):
    op=BRepAlgoAPI_Fuse(fused,s); op.Build()
    if not op.IsDone(): raise RuntimeError("fuse failed at %d"%i)
    fused=op.Shape()
print("steel fused OK")

# air = sphere - steel
sph=BRepPrimAPI_MakeSphere(gp_Pnt(0,0,0),R_SPHERE).Shape()
cut=BRepAlgoAPI_Cut(sph,fused); cut.Build()
if not cut.IsDone(): raise RuntimeError("air cut failed")
air=cut.Shape()
print("air cut OK; air solids:",len(solids(air)))

# compound: steel + air
comp=TopoDS_Compound(); bb=BRep_Builder(); bb.MakeCompound(comp)
bb.Add(comp,fused); bb.Add(comp,air)
w=STEPControl_Writer(); w.Transfer(comp,STEPControl_StepModelType.STEPControl_AsIs); w.Write(STEP_OUT)
print("wrote",STEP_OUT)
