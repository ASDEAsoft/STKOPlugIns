# standard modules
import importlib
import os
import tclimport.tcldoc
importlib.reload(tclimport.tcldoc)
from tclimport.tcldoc import tcldoc
from tclimport.mathutils import *
# PyMpc
from PyMpc import *
# Qt
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
	QApplication,
	QDialog,
	QFileDialog,
	QVBoxLayout,
	QLabel,
	QProgressBar,
	QPlainTextEdit)

# A custom waiting dialog
class WDialog(QDialog):
	def __init__(self, parent=None):
		QDialog.__init__(self, parent)
		self.canClose = False
		# build dialog
		self.setWindowTitle('Importing Tcl...')
		lo = QVBoxLayout(self)
		self.label = QLabel('Process status:')
		lo.addWidget(self.label)
		self.pbar = QProgressBar(self)
		self.pbar.setMaximum(100)
		self.pbar.setTextVisible(True)
		self.pbar.setValue(0)
		lo.addWidget(self.pbar)
		self.log = QPlainTextEdit(self)
		self.log.setReadOnly(True)
		self.log.setOverwriteMode(False)
		self.log.setTextInteractionFlags(Qt.NoTextInteraction)
		self.log.setUndoRedoEnabled(False)
		self.log.setBackgroundVisible(False)
		lo.addWidget(self.log)
		self.setLayout(lo)
	# handle the non-closable stuff
	def closeNow(self):
		self.canClose = True
		self.accept()
		self.canClose = False
	def accept(self):
		if self.canClose:
			QDialog.accept(self)
	def reject(self):
		if self.canClose:
			QDialog.reject(self)
	# print
	def sendLog(self, perc, info):
		self.pbar.setValue(max(0, min(100, int(perc*100.0))))
		self.log.appendPlainText(info)

# Create a new CAE document
App.runCommand("NewPreDocument")

# clear the terminal window
App.clearTerminal()
# get a reference to the current CAE document
doc = App.caeDocument()

# let the user choose the file to import
filename = QFileDialog.getOpenFileName(
	QApplication.activeWindow(), "Open TCL File", ".", 
	"TCL File (*.tcl);;Text files (*.txt);;All files (*.*)")
if filename is not None:
	filename = filename[0]
else:
	raise Exception("No file provided")

with open(filename, 'r') as file:
	filecontents = file.read()

def select(obj):
	if isinstance(obj, MpcInteraction):
		# interaction
		doc.scene.select(obj)
	else:
		# geometry sub
		doc.scene.select(obj[0], obj[1], obj[2])

def unselect():
	doc.scene.unselectAll()

def runCommand_function(a, b=''):
	App.runCommand(a, b)
runCommand = runCommand_function

def add_geom_function(geom):
	doc.addGeometry(geom)
	doc.commitChanges()
	doc.dirty = True
	App.processEvents()
add_geom = add_geom_function

def find_geom(subshape):
	# subshape is a FxOccShape.
	# returns a tuple (MpcGeometry: geom, int: subshape id, MpcSubshapeType: type)
	results = []
	for _, geom in doc.geometries.items():
		ok, type, subshape_id = geom.shape.find(subshape)
		if ok:
			results.append((geom, subshape_id, type))
	if len(results) == 0:
		raise Exception("Cannot find a geometry from the provided subshape")
	return results

def add_inter_function(inter):
	doc.addInteraction(inter)
	doc.commitChanges()
	doc.dirty = True
	App.processEvents()
add_inter = add_inter_function

def add_selection_set_function(sset):
	doc.addSelectionSet(sset)
	doc.commitChanges()
	doc.dirty = True
add_selection_set = add_selection_set_function

def add_locax_function(transf):
	bbox = doc.scene.boundingBox
	p0 = Math.vec3(bbox.minPoint.x - bbox.maxSize*0.35, 0.0, 0.0)
	if transf.vz is not None:
		dz = v3normalize(transf.vz)
		if dz[2] > 0.99:
			p1 = p0 + Math.vec3(1.0, 0.0, 0.0)
			p2 = p0 + Math.vec3(0.0, 1.0, 0.0)
		elif dz[2] < -0.99:
			p1 = p0 + Math.vec3(1.0, 0.0, 0.0)
			p2 = p0 - Math.vec3(0.0, 1.0, 0.0)
		else:
			dx = (0.0, 0.0, 1.0)
			dy = v3normalize(v3cross(dz, dx))
			dx = v3normalize(v3cross(dy, dz))
			p1 = p0 + Math.vec3(dx[0], dx[1], dx[2])
			p2 = p0 + Math.vec3(dy[0], dy[1], dy[2])
		locax = MpcLocalAxes(
			transf.id, 
			'LocalAxes Vz:({:.3g}, {:.3g}, {:.3g})'.format(dz[0], dz[1], dz[2]),
			MpcLocalAxesType.Rectangular,
			p0, p1, p2
			)
	else:
		dx = v3normalize(transf.vx)
		dy = v3normalize(transf.vy)
		dz = v3normalize(v3cross(dx, dy))
		dy = v3normalize(v3cross(dz, dx))
		p1 = p0 + Math.vec3(dx[0], dx[1], dx[2])
		p2 = p0 + Math.vec3(dy[0], dy[1], dy[2])
		locax = MpcLocalAxes(
			transf.id, 
			'LocalAxes Vx:({:.3g}, {:.3g}, {:.3g}) Vy:({:.3g}, {:.3g}, {:.3g})'.format(
			dx[0], dx[1], dx[2], dy[0], dy[1], dy[2]),
			MpcLocalAxesType.Rectangular,
			p0, p1, p2
			)
	doc.addLocalAxes(locax)
	doc.commitChanges()
	doc.dirty = True
	App.processEvents()
add_locax = add_locax_function

def get_meta_phys_prop(name):
	return doc.metaDataPhysicalProperty(name)

def get_phys_prop(id):
	return doc.getPhysicalProperty(id)

def add_phys_prop_function(p):
	doc.addPhysicalProperty(p)
	doc.commitChanges()
	doc.dirty = True
add_phys_prop = add_phys_prop_function

def get_meta_elem_prop(name):
	return doc.metaDataElementProperty(name)

def get_elem_prop(id):
	return doc.getElementProperty(id)

def add_elem_prop_function(p):
	doc.addElementProperty(p)
	doc.commitChanges()
	doc.dirty = True
add_elem_prop = add_elem_prop_function

def get_meta_cond(name):
	return doc.metaDataCondition(name)

def get_cond(id):
	return doc.getCondition(id)

def add_cond_function(p):
	doc.addCondition(p)
	doc.commitChanges()
	doc.dirty = True
add_cond = add_cond_function

def get_meta_step(name):
	return doc.metaDataAnalysisStep(name)

def get_step(id):
	return doc.getAnalysisStep(id)

def add_step_function(p):
	doc.addAnalysisStep(p)
	doc.commitChanges()
	doc.dirty = True
add_step = add_step_function

def get_meta_def(name):
	return doc.metaDataDefinition(name)

def get_def(id):
	return doc.getDefinition(id)

def add_def_function(p):
	doc.addDefinition(p)
	doc.commitChanges()
	doc.dirty = True
add_def = add_def_function

def make_default_mesh():
	mc = doc.meshControls
	for _, gc in mc.geometryControls.items():
		for fc in gc.faceControls:
			fc.topology = MpcMeshAlgoTopology.QuadHexa
			fc.algorithm = MpcMeshAlgo.Structured
	App.runCommand("BuildMesh")

# open the dialog
dialog = WDialog(QApplication.activeWindow())
dialog.show()
def send_log(perc, info):
	dialog.sendLog(perc, info)

# put them all in a callbacks class as required by tcldoc_t
class callbacks:
	def __init__(self):
		self.select = select
		self.unselect = unselect
		self.runCommand = runCommand
		self.addGeometry = add_geom
		self.findGeometry = find_geom
		self.addInteraction = add_inter
		self.addSelectionSet = add_selection_set
		self.addLocalAxes = add_locax
		self.getMetaPhysicalProperty = get_meta_phys_prop
		self.getPhysicalProperty = get_phys_prop
		self.addPhysicalProperty = add_phys_prop
		self.getMetaElementProperty = get_meta_elem_prop
		self.getElementProperty = get_elem_prop
		self.addElementProperty = add_elem_prop
		self.getMetaCondition = get_meta_cond
		self.getCondition = get_cond
		self.addCondition = add_cond
		self.getMetaAnalysisStep = get_meta_step
		self.getAnalysisStep = get_step
		self.addAnalysisStep = add_step
		self.getMetaDefinition = get_meta_def
		self.getDefinition = get_def
		self.addDefinition = add_def
		self.sendLog = send_log

# run
try:
	dialog.sendLog(0.0, 'Parsing TCL file...')
	App.processEvents()
	tcl = tcldoc()
	tcl.parse(filecontents)
	dialog.sendLog(10.0, 'Building STKO model...')
	App.processEvents()
	tcl.build(callbacks())
	doc.randomizeMaterialColors()
	dialog.sendLog(90.0, 'Meshing...')
	App.processEvents()
	make_default_mesh()
	App.runCommand("Regenerate", "2")
	App.updateActiveView()
	dialog.sendLog(100.0, 'Done')
	print(tcl)
	IO.write_clog('File "{}" correctly imported'.format(os.path.basename(filename)))
except:
	import traceback
	IO.write_cerr(traceback.format_exc())
finally:
	dialog.closeNow()
	dialog.deleteLater()