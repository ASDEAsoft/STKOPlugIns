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
last_step = all_steps[933]
opt = MpcOdbVirtualResultEvaluationOptions()
opt.stage = last_stage
opt.step = last_step
field = R.evaluate(opt)
ufield = U.evaluate(opt)
mesh = field.mesh

ele_id = doc.scene.plotSelection[1].info.elements[0]
ele = mesh.getElement(ele_id)

