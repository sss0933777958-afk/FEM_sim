"""
sensor_misplace discrete sphere-grid -> mm STEP for CAD check (deliver-step-for-check.md).

Source = ANSYS IGESOUT from MT_Geom_SensorMisplace.txt (steel MAT_MT + ~129 R0.3 grid spheres),
written to the claude scratchpad. Auto-calibrate scale to mm, write explicit mm STEP.

Grid (per pole): fore-aft line (phi=0, soff=4.572+k*1mm, k=-2..10 =13) + azimuth ring
(soff=4.572, phi=m*30deg; lower m=-3..3 / upper m=0..11, skip m=0). Lower 19, upper 24 -> ~129.
Verifies: sphere count, each sphere center matches P(soff,phi)+(H/2)n_hat, min inter-sphere gap (>=0.6mm).

Run:  python make_sensor_misplace_step.py
"""
import math
from OCP.IGESControl import IGESControl_Reader
from OCP.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCP.Interface import Interface_Static
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.gp import gp_Trsf
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_SOLID
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.TopoDS import TopoDS_Compound
from OCP.BRep import BRep_Builder

ROOT = r"G:\my_workspace\code\FEM_sim\magnetic_sim\ANSYS\main"
SCR  = r"C:\Users\Kuo\AppData\Local\Temp\claude\G--my-workspace\b77e29dd-bd9c-479e-93f7-d3abb6bbd5e1\scratchpad"
IGES = SCR + r"\sensor_misplace.iges"
OUT  = ROOT + r"\model_check\long_fei\sensor_misplace.step"

D2R = math.pi/180.0
PSI=35.2644; INC=36.59; BETA=11.3099; RN=0.5; AIR=0.41; HH=0.10; OFFN=AIR+HH/2; SPHOF=-13.0+RN/math.sqrt(3.0)
SOFF0=4.572; STEP_FA=1.0; STEP_AZ=30
RSHANK=3.047; ABASE_LO=36; ABASE_HI=36
PTH=[0,180,120,300,60,240]; PLO=[1,0,1,0,0,1]

def dirv(e,a): return (math.cos(e*D2R)*math.cos(a*D2R), math.cos(e*D2R)*math.sin(a*D2R), math.sin(e*D2R))
def add(*v): return tuple(sum(c) for c in zip(*v))
def sca(k,v): return (k*v[0],k*v[1],k*v[2])
def cross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def expected_centers():
    cb=math.cos(BETA*D2R); sb=math.sin(BETA*D2R); pts=[]
    for ip in range(6):
        th=PTH[ip]; low=PLO[ip]
        if low: elev=-PSI; ax=dirv(0,th); e2a=-BETA; mlo,mhi=-3,3
        else:   elev=+PSI; ax=dirv(INC,th); e2a=INC+BETA; mlo,mhi=0,11
        tip=(RN*math.cos(elev*D2R)*math.cos(th*D2R), RN*math.cos(elev*D2R)*math.sin(th*D2R), RN*math.sin(elev*D2R)+SPHOF)
        e2=dirv(e2a,th)
        uu=sca(1.0/sb, add(e2, sca(-cb,ax))); vv=cross(ax,uu)
        def cen(soff,phi):
            w=add(sca(math.cos(phi*D2R),uu), sca(math.sin(phi*D2R),vv))
            e=add(sca(cb,ax), sca(sb,w)); n=add(sca(-sb,ax), sca(cb,w))
            return add(tip, sca(soff,e), sca(OFFN,n))
        ablo = ABASE_LO if low else ABASE_HI
        for k in range(-2,11): pts.append((ip+1,"fa",cen(SOFF0+k*STEP_FA,0)))
        for aa in range(15, ablo+1):                       # shank：軸點 + 徑向 û
            pts.append((ip+1,"sh", add(tip, sca(float(aa),ax), sca(RSHANK+OFFN,uu))))
        for m in range(mlo,mhi+1):
            if m!=0: pts.append((ip+1,"az",cen(SOFF0,m*STEP_AZ)))
    return pts

def solids(sh):
    out=[]; e=TopExp_Explorer(sh,TopAbs_SOLID)
    while e.More(): out.append(e.Current()); e.Next()
    return out
def vol(s):
    g=GProp_GProps(); BRepGProp.VolumeProperties_s(s,g); return g.Mass()
def cent(s):
    g=GProp_GProps(); BRepGProp.VolumeProperties_s(s,g); p=g.CentreOfMass(); return (p.X(),p.Y(),p.Z())
def bbox(sh):
    b=Bnd_Box(); BRepBndLib.Add_s(sh,b); return b.Get()

# read + scale
r=IGESControl_Reader(); r.ReadFile(IGES); r.TransferRoots(); sh=r.OneShape()
xm,ym,zm,xM,yM,zM=bbox(sh); span=max(xM-xm,yM-ym,zM-zm)
sc=1.0/25.4 if span>400 else (1000.0 if span<1 else 1.0)
if sc!=1.0:
    t=gp_Trsf(); t.SetScaleFactor(sc); sh=BRepBuilderAPI_Transform(sh,t,True).Shape()
print("read solids=%d span=%.1f scale=%.5f"%(len(solids(sh)),span,sc))
Interface_Static.SetCVal_s("write.step.unit","MM")
w=STEPControl_Writer(); w.Transfer(sh,STEPControl_AsIs); w.Write(OUT)
print("wrote",OUT)

# verify spheres
rd=STEPControl_Reader(); rd.ReadFile(OUT); rd.TransferRoots(); s2=rd.OneShape()
sv=solids(s2)
sph=[s for s in sv if vol(s)<0.5]        # sphere ~0.113 mm^3
exp=expected_centers()
print("spheres found=%d  (expect %d) ; steel solids=%d"%(len(sph),len(exp),len(sv)-len(sph)))
# match each expected to nearest sphere centroid
scen=[cent(s) for s in sph]
maxd=0; worst=None
for ip,tag,ec in exp:
    d,_=min(((( (ec[0]-c[0])**2+(ec[1]-c[1])**2+(ec[2]-c[2])**2)**0.5),c) for c in scen)
    if d>maxd: maxd=d; worst=(ip,tag,ec)
print("max center mismatch (expected->nearest sphere) = %.4f mm  @P%d %s"%(maxd,worst[0],worst[1]))
# min inter-sphere gap
mind=1e9
for i in range(len(scen)):
    for j in range(i+1,len(scen)):
        a,b=scen[i],scen[j]
        d=((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5
        if d<mind: mind=d
print("min inter-sphere center gap = %.4f mm  (need >=0.60 to avoid overlap)"%mind)
vv=[vol(s) for s in sph]
print("sphere vol range = %.4f .. %.4f mm^3 (expect ~0.113 unclipped)"%(min(vv),max(vv)))
# sphere -> steel clearance (catch shank base collision; markers not booleaned so volume won't show it)
steel=[s for s in sv if vol(s)>=0.5]
comp=TopoDS_Compound(); bb=BRep_Builder(); bb.MakeCompound(comp)
for s in steel: bb.Add(comp,s)
WP=(0.0,0.0,SPHOF)
far=[s for s in sph if ((cent(s)[0]-WP[0])**2+(cent(s)[1]-WP[1])**2+(cent(s)[2]-WP[2])**2)**0.5 > 20.0]  # 只查遠場(shank 基座區)球
mind_st=1e9; worst=None
for s in far:
    d=BRepExtrema_DistShapeShape(s,comp).Value()
    if d<mind_st: mind_st=d; worst=cent(s)
print("min (far-shank) sphere->steel = %.4f mm @ (%.2f,%.2f,%.2f) over %d far spheres (>=~0.10 clears base)"%(mind_st,worst[0],worst[1],worst[2],len(far)))
