from PyMpc import *
from PyMpc import MpcOdbVirtualResult as vr
import math

App.clearTerminal()
doc = App.postDocument()
db = doc.getDatabase(1)
R = db.getElementalResult('force', match=MpcOdb.StartsWith)
U = db.getNodalResult('Displacement')
all_stages = db.getStageIDs()
last_stage = all_stages[-1]
all_steps = db.getStepIDs(last_stage)
last_step = all_steps[-1]
opt = MpcOdbVirtualResultEvaluationOptions()
opt.stage = last_stage
opt.step = last_step
field = R.evaluate(opt)
ufield = U.evaluate(opt)
mesh = field.mesh

ele_id = doc.scene.plotSelection[1].info.elements[0]
ele = mesh.getElement(ele_id)
f = []
m = []
u = {}
for node in ele.nodes:
	row = MpcOdbResultField.node(node.id)
	iu = ufield[row]
	u[node.id] = Math.vec3(iu[0], iu[1], 0.0)*1
for i in range(len(ele.nodes)):
	row = MpcOdbResultField.gauss(ele.id, i)
	ir = field[row]
	f.append(Math.vec3(ir[0], ir[1], 0.0))
	m.append(Math.vec3(0.0, 0.0, ir[2]))
F = Math.vec3(0.0,0.0,0.0)
M = Math.vec3(0.0,0.0,0.0)
c = ele.nodes[0].position + u[ele.nodes[0].id]
for iF, iM, node in zip(f, m, ele.nodes):
	F += iF
	M += iM
	d = node.position+u[node.id] - c
	iFM = d.cross(iF)
	M += iFM
print(*F)
print(*M)

print("===========================")
def _pv(x):
	print('['+','.join('{:10.3g}'.format(i if abs(i)>1e-8 else 0) for i in x)+']')
F1 = f[0]
M1 = m[0]
F2 = f[1]
M2 = m[1]
X1 = ele.nodes[0].position + u[ele.nodes[0].id]
X2 = ele.nodes[1].position + u[ele.nodes[1].id]
XC = (X1+X2)/2.0
R1 = Math.vec(6, 0.0)
R2 = Math.vec(6, 0.0)
K = Math.mat(6, 6, 0.0)
# eq X1-XC
print('eq X1-XC')
for i in range(3):
	K[i,i] = -1.0
	K[i+3,i+3] = -1.0
	R1[i] = F1[i]
	R1[i+3] = M1[i]
	R2[i] = F2[i]
	R2[i+3] = M2[i]
D1 = XC-X1
D2 = XC-X2
_pv(R1)
_pv(R2)
K[3,1]= D1.z; K[3,2]=-D1.y
K[4,0]=-D1.z; K[4,2]= D1.x
K[5,0]= D1.y; K[5,1]=-D1.x
FC1 = K.solve(R1)
_pv(FC1)
_pv(K*FC1 - R1)
K[3,1]= D2.z; K[3,2]=-D2.y
K[4,0]=-D2.z; K[4,2]= D2.x
K[5,0]= D2.y; K[5,1]=-D2.x
FC2 = K.solve(R2)
_pv(FC2)
_pv(K*FC2 - R2)
mm1 = -M1[2]
mm2 =  M2[2]
print((mm1+mm2)/2)