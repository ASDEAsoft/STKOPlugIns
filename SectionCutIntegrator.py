'''
@ todo:
now we avg the stresses.
in case of more gauss points we can extrapolate to nodes ans interpolate to cut integration points?
'''

import os
import math
import time
import itertools
from collections import defaultdict
from PyMpc import *
from PyMpc import MpcOdbVirtualResult as vr
from PySide2.QtCore import (Qt, QObject, QEventLoop, Slot, QTimer)
from PySide2.QtWidgets import (
	QApplication,
	QMessageBox,
	QInputDialog,
	QLineEdit,
	QProgressBar,
	QDialog,
	QVBoxLayout,
	QGridLayout,
	QListWidget,
	QListWidgetItem,
	QLabel,
	QDialogButtonBox,
	QPushButton,
	QCheckBox,
	QWidget,
	QFrame,
	QSpacerItem,
	QSizePolicy,
	)

# some globals vars
_WIN = QApplication.activeWindow()
_TITLE = QApplication.applicationName()
_TIMEOUT = 60*1000

def _asign(x):
	if x == 0.0:
		return (0.0, 0.0)
	elif x > 0.0:
		return (x, 1.0)
	else:
		return (-x, -1.0)
def _acos(x):
	return math.acos(max(-1.0, min(1.0, x)))

class _TP:
	tol = 1e-10
	def __init__(self, p):
		self.p = p
	def __hash__(self):
		return hash((int(self.p.x/_TP.tol), int(self.p.y/_TP.tol), int(self.p.z/_TP.tol)))
	def __eq__(self, b):
		return (int(self.p.x/_TP.tol), int(self.p.y/_TP.tol), int(self.p.z/_TP.tol)) == (int(b.p.x/_TP.tol), int(b.p.y/_TP.tol), int(b.p.z/_TP.tol))
	def __str__(self):
		return str('({},{},{})'.format(self.p.x, self.p.y, self.p.z))
	def __repr__(self):
		return str(self)
class _TPID:
	def __init__(self):
		self.i = 0
	def __str__(self):
		return str(self.i)
	def __repr__(self):
		return str(self)

class _CPlane:
	def __init__(self, mesh, doc = None):
		self.infinite = True
		self.center = Math.vec3(0.0,0.0,0.0) # geometric center
		self.sizey = 0.0 # size in the Y direction
		self.sizez = 0.0 # side in the Z direction
		self.dx = Math.vec3(1.0,0.0,0.0) # local x axis
		self.dy = Math.vec3(0.0,1.0,0.0) # local y axis
		self.dz = Math.vec3(0.0,0.0,1.0) # local z axis
		self.vrep = None
		self.vrep_face = None
		self.doc = doc
		lch = 0.0
		if len(mesh.nodes) > 0:
			bbox = FxBndBox()
			for _, node in mesh.nodes.items():
				bbox.add(node.position)
				self.center += node.position
			lch = (bbox.maxPoint - bbox.minPoint).norm()
			self.center = (bbox.maxPoint + bbox.minPoint)/2.0
		self.setSize(lch, lch) # initial size...
		self.draw()
	def setSize(self, sy, sz):
		self.sizey = sy
		self.sizez = sz
		smean = (sy + sz)/2.0
		self.tolerance = max(1.0e-8, 1.0e-6*smean)
		_TP.tol = self.tolerance*10.0
	def draw(self):
		if self.doc is None:
			return None
		if self.vrep:
			self.doc.removeCustomDrawableEntity(self.vrep)
		if self.vrep_face:
			self.doc.removeCustomDrawableEntity(self.vrep_face)
		self.vrep = FxShape()
		self.vrep_face = FxShape()
		# make triad
		lch = (self.sizey + self.sizez)/2.0
		for dir,col in zip( (self.dx,self.dy,self.dz), 
				(Math.vec3(0.8,0.0,0.0), Math.vec3(0.0,0.8,0.0), Math.vec3(0.0,0.0,0.8)) ):
			edge = FxShapeEdge()
			edge.indices.append(0)
			edge.indices.append(1)
			edge.vertices.append(Math.vertex(self.center, col))
			edge.vertices.append(Math.vertex(self.center + dir * lch/4.0, col))
			self.vrep.edges.append(edge)
		# make rectangle
		edge = FxShapeEdge()
		for i in range(4): edge.indices.append(i)
		edge.indices.append(0)
		for i,j in ((-1,-1),(1,-1),(1,1),(-1,1)):
			yl = i*self.sizey/2.0
			zl = j*self.sizez/2.0
			edge.vertices.append(Math.vertex(
				Math.vec3(
				self.dy.x*yl + self.dz.x*zl,
				self.dy.y*yl + self.dz.y*zl,
				self.dy.z*yl + self.dz.z*zl
				) + self.center, Math.vec3(0.0,0.0,0.0)))
		self.vrep.edges.append(edge)
		# make material for borders
		mat = FxMaterial()
		mat.lineWidth = 2.0
		mat.coloringMode = FxColoringMode.ColorizeUsingNormals
		mat.lighting = False
		mat.depthTest = False
		mat.transparency = 0.8
		mat.transparencyOnlyOnFaces = False
		self.vrep.material = mat
		self.vrep.commitChanges()
		# make face
		face = FxShapeFace()
		for v in edge.vertices: face.vertices.append(v)
		face.triangles.append(Math.triangle(0,1,2))
		face.triangles.append(Math.triangle(0,2,3))
		self.vrep_face.faces.append(face)
		mat = FxMaterial()
		mat.transparency = 0.4
		self.vrep_face.material = mat
		self.vrep_face.commitChanges()
		# add
		self.doc.addCustomDrawableEntity(self.vrep)
		self.doc.addCustomDrawableEntity(self.vrep_face)
		App.updateActiveView()

class _NonModalDialogBase(QDialog):
	def __init__(self,parent=None):
		super().__init__(parent=parent)
		self.setWindowFlags(self.windowFlags() | Qt.WindowStaysOnTopHint)
		l = QVBoxLayout()
		self.centralWidget = QWidget()
		bbox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		bbox.accepted.connect(self.accept)
		bbox.rejected.connect(self.reject)
		hline = QFrame()
		hline.setFrameShape(QFrame.HLine)
		hline.setFrameShadow(QFrame.Sunken)
		l.addWidget(self.centralWidget)
		l.addWidget(hline)
		l.addWidget(bbox)
		self.setLayout(l)
		self.loop = QEventLoop()
		self.outCode = 0
	def start(self):
		self.show()
		self.loop.exec_()
		return self.outCode
	def done(self, value):
		self.loop.quit()
		self.outCode = value
		return super().done(value)

class _CPlaneDialog(_NonModalDialogBase):
	def __init__(self, mesh, doc=None, parent=None):
		super().__init__(parent=parent)
		self.cplane = _CPlane(mesh, doc=doc)
		self.centralWidget.setLayout(QGridLayout())
		#self.centralWidget.layout().setColumnStretch(0,0)
		#self.centralWidget.layout().setColumnStretch(1,1)
		row = 0
		# infinite
		self.cbox_inf = QCheckBox('Infinite Plane')
		self.cbox_inf.setChecked(self.cplane.infinite)
		self.centralWidget.layout().addWidget(self.cbox_inf, row, 0, 1, 2)
		row += 1
		# labels
		self.label_center = QLabel()
		self.centralWidget.layout().addWidget(self.label_center, row, 0)
		self.label_dx = QLabel()
		self.centralWidget.layout().addWidget(self.label_dx, row, 1)
		row += 1
		self.label_size_y = QLabel()
		self.centralWidget.layout().addWidget(self.label_size_y, row, 0)
		self.label_dy = QLabel()
		self.centralWidget.layout().addWidget(self.label_dy, row, 1)
		row += 1
		self.label_size_z = QLabel()
		self.centralWidget.layout().addWidget(self.label_size_z, row, 0)
		self.label_dz = QLabel()
		self.centralWidget.layout().addWidget(self.label_dz, row, 1)
		row += 1
		# basic buttons
		self.btn_center = QPushButton('Set Center')
		self.centralWidget.layout().addWidget(self.btn_center,   row, 0)
		self.btn_dx = QPushButton('Set dX')
		self.centralWidget.layout().addWidget(self.btn_dx,   row, 1)
		row += 1
		self.btn_size_y = QPushButton('Set Y-Size')
		self.centralWidget.layout().addWidget(self.btn_size_y,   row, 0)
		self.btn_dy = QPushButton('Set dY')
		self.centralWidget.layout().addWidget(self.btn_dy,   row, 1)
		row += 1
		self.btn_size_z = QPushButton('Set Z-Size')
		self.centralWidget.layout().addWidget(self.btn_size_z,   row, 0)
		self.btn_dz = QPushButton('Set dZ')
		self.centralWidget.layout().addWidget(self.btn_dz,   row, 1)
		row += 1
		# spacer
		self.centralWidget.layout().addItem(QSpacerItem(1, 1, QSizePolicy.Minimum, QSizePolicy.Expanding), row, 0, 1, 2)
		row += 1
		# init labels
		self._updateCenter()
		self._updateSizeY()
		self._updateSizeZ()
		self._updateOrientation()
		# connections
		self.cbox_inf.toggled.connect(self.onInfiniteToggled)
		self.btn_center.clicked.connect(self.onCenterClicked)
		self.btn_size_y.clicked.connect(self.onSizeYClicked)
		self.btn_size_z.clicked.connect(self.onSizeZClicked)
		self.btn_dx.clicked.connect(self.onDxClicked)
		self.btn_dy.clicked.connect(self.onDyClicked)
		self.btn_dz.clicked.connect(self.onDzClicked)
		# done
		print('4. CUT-PLANE: Edit cut plane location, size and orientation.')
		print('   Press Ok to accept [or Cancel to abort]')
	def _lock(self):
		self.setEnabled(False)
		self.setWindowOpacity(0.5)
	def _unlock(self):
		self.setEnabled(True)
		self.setWindowOpacity(1.0)
	def _hitTest(self):
		doc = self.cplane.doc
		watcher = _MovingPointTracer(doc)
		watcher.begin()
		point, done, num_comp_inp = doc.scene.waitForHitPoint()
		watcher.end()
		doc.scene.unselectAll()
		return (point, done, num_comp_inp)
	def _updateCenter(self):
		self.label_center.setText('C = ({:.4g}, {:.4g}, {:.4g})'.format(
			self.cplane.center.x, self.cplane.center.y, self.cplane.center.z))
	def _updateSizeY(self):
		self.label_size_y.setText('LY = {:.4g}'.format(self.cplane.sizey))
	def _updateSizeZ(self):
		self.label_size_z.setText('LZ = {:.4g}'.format(self.cplane.sizez))
	def _updateOrientation(self):
		self.label_dx.setText('dX = ({:.4g}, {:.4g}, {:.4g})'.format(
			self.cplane.dx.x, self.cplane.dx.y, self.cplane.dx.z))
		self.label_dy.setText('dY = ({:.4g}, {:.4g}, {:.4g})'.format(
			self.cplane.dy.x, self.cplane.dy.y, self.cplane.dy.z))
		self.label_dz.setText('dZ = ({:.4g}, {:.4g}, {:.4g})'.format(
			self.cplane.dz.x, self.cplane.dz.y, self.cplane.dz.z))
	@Slot(bool)
	def onInfiniteToggled(self, checked):
		try:
			self.cplane.infinite = checked
		except Exception as ex:
			print(ex)
	@Slot()
	def onCenterClicked(self):
		try:
			self._lock()
			print('   4.1 CENTER: click on a point to edit the center of the cut-plane\n'
				'      or type X Y <Z> [+ Enter] coordinates.')
			point, done, num_comp_inp = self._hitTest()
			if done == MpcSceneHitCode.AcceptInput and num_comp_inp < 2:
				print('      > Invalid input')
				return None
			if done >= 0:
				self.cplane.center.x = point.x
				self.cplane.center.y = point.y
				self.cplane.center.z = point.z
				self.cplane.draw()
				self._updateCenter()
		except Exception as ex:
			print(ex)
		finally:
			self._unlock()
	def _onSizeInternal(self, yoz):
		try:
			if yoz == 1:
				vstep = 2
				label = 'Y'
				_updateFun = self._updateSizeY
				def _sizefun(x):
					self.cplane.setSize(x, self.cplane.sizez)
			else:
				vstep = 3
				label = 'Z'
				_updateFun = self._updateSizeZ
				def _sizefun(x):
					self.cplane.setSize(self.cplane.sizey, x)
			self._lock()
			print('   4.{} SIZE {}: click on the first point or type X Y <Z> [+ Enter] coordinates.\n'
				'      or type the distance L'.format(vstep, label))
			point, done, num_comp_inp = self._hitTest()
			if done < 0:
				return None
			if done == MpcSceneHitCode.AcceptInput:
				if num_comp_inp < 1:
					print('      > Invalid input')
					return None
				if num_comp_inp == 1:
					_sizefun(point.x)
					self.cplane.draw()
					_updateFun()
					return None
			P1 = point
			print('   4.{} SIZE {}: click on the second point or type X Y <Z> [+ Enter] coordinates.'.format(
				vstep, label))
			point, done, num_comp_inp = self._hitTest()
			if done < 0:
				return None
			if done == MpcSceneHitCode.AcceptInput and num_comp_inp < 2:
				print('      > Invalid input')
				return None
			P2 = point
			L = (P2 - P1).norm()
			if L == 0.0:
				print('      > Cannot set size = 0')
				return None
			_sizefun(L)
			self.cplane.draw()
			_updateFun()
		except Exception as ex:
			print(ex)
		finally:
			self._unlock()
	@Slot()
	def onSizeYClicked(self):
		self._onSizeInternal(1)
	@Slot()
	def onSizeZClicked(self):
		self._onSizeInternal(2)
	def _onDirInternal(self, d):
		try:
			vstep = d + 3
			label = ('X','Y','Z')[d-1]
			self._lock()
			print('   5.{} DIRECTION {}: click on the first point or type X Y <Z> [+ Enter] components of the direction vector.'.format(vstep, label))
			point, done, num_comp_inp = self._hitTest()
			direction = None
			if done < 0:
				return None
			if done == MpcSceneHitCode.AcceptInput:
				if num_comp_inp < 1:
					print('      > Invalid input')
					return None
				direction = point.normalized()
			if not direction:
				P1 = point
				print('   5.{} DIRECTION {}: click on the second point or type X Y <Z> [+ Enter] coordinates.'.format(vstep, label))
				point, done, num_comp_inp = self._hitTest()
				if done < 0:
					return None
				if done == MpcSceneHitCode.AcceptInput and num_comp_inp < 2:
					print('      > Invalid input')
					return None
				P2 = point
				direction = (P2-P1).normalized()
			if direction.norm() == 0.0:
				print('      > Invalid input: direction has 0 norm')
				return None
			return direction
		except Exception as ex:
			print(ex)
		finally:
			self._unlock()
	@Slot()
	def onDxClicked(self):
		dx = self._onDirInternal(1)
		angle = _acos(self.cplane.dx.dot(dx))
		if angle == 0.0:
			print('      No transformation (0° angle)')
			return None
		axis = self.cplane.dx.cross(dx)
		q = Math.quaternion.fromAxisAngle(axis, angle)
		self.cplane.dx = q.rotate(self.cplane.dx).normalized()
		self.cplane.dy = q.rotate(self.cplane.dy).normalized()
		self.cplane.dz = q.rotate(self.cplane.dz).normalized()
		self.cplane.draw()
		self._updateOrientation()
	@Slot()
	def onDyClicked(self):
		dy = self._onDirInternal(2)
		angle = _acos(self.cplane.dy.dot(dy))
		if angle == 0.0:
			print('      No transformation (0° angle)')
			return None
		axis = self.cplane.dy.cross(dy)
		q = Math.quaternion.fromAxisAngle(axis, angle)
		self.cplane.dx = q.rotate(self.cplane.dx).normalized()
		self.cplane.dy = q.rotate(self.cplane.dy).normalized()
		self.cplane.dz = q.rotate(self.cplane.dz).normalized()
		self.cplane.draw()
		self._updateOrientation()
	@Slot()
	def onDzClicked(self):
		dz = self._onDirInternal(3)
		angle = _acos(self.cplane.dz.dot(dz))
		if angle == 0.0:
			print('      No transformation (0° angle)')
			return None
		axis = self.cplane.dz.cross(dz)
		q = Math.quaternion.fromAxisAngle(axis, angle)
		self.cplane.dx = q.rotate(self.cplane.dx).normalized()
		self.cplane.dy = q.rotate(self.cplane.dy).normalized()
		self.cplane.dz = q.rotate(self.cplane.dz).normalized()
		self.cplane.draw()
		self._updateOrientation()

class _SelectedMeshGraphics:
	def __init__(self):
		self.vrep_glow = None
		self.vrep = None
	def _base_color():
		mat = FxMaterial()
		mat.lighting = False
		mat.depthTest = False
		mat.transparency = 0.99
		mat.transparencyOnlyOnFaces = False
		mat.edgeColor = FxColor(0.1, 0.8, 1.0)
		return mat
	def clear(self, doc):
		if self.vrep:
			doc.removeCustomDrawableEntity(self.vrep)
		if self.vrep_glow:
			doc.removeCustomDrawableEntity(self.vrep_glow)
		App.updateActiveView()
	def draw(self, doc, mesh):
		self.clear(doc)
		mesh.buildBoundaries()
		vrep = mesh.makeVisualRepresentation()
		self.vrep_glow = FxShape()
		self.vrep = FxShape()
		for edge in vrep.edges:
			self.vrep_glow.edges.append(edge)
			self.vrep.edges.append(edge)
		for edge in vrep.meshEdges:
			self.vrep_glow.meshEdges.append(edge)
			self.vrep.meshEdges.append(edge)
		# glow material
		self.vrep.material = _SelectedMeshGraphics._base_color()
		self.vrep_glow.material = _SelectedMeshGraphics._base_color()
		self.vrep_glow.material.lineWidth = 4.0
		self.vrep_glow.material.transparency = 0.4
		doc.addCustomDrawableEntity(self.vrep_glow)
		doc.addCustomDrawableEntity(self.vrep)
		App.updateActiveView()

class _CutElementsGraphics:
	def __init__(self):
		self.vrep_glow = None
		self.vrep = None
	def _base_color():
		mat = FxMaterial()
		mat.lighting = False
		mat.depthTest = False
		mat.transparency = 0.99
		mat.transparencyOnlyOnFaces = False
		mat.edgeColor = FxColor(1.0, 0.0, 0.0)
		mat.lineWidth = 2.0
		return mat
	def _glow_color():
		mat = _CutElementsGraphics._base_color()
		mat.lineWidth = 6.0
		mat.transparency = 0.4
		return mat
	def clear(self, doc):
		if self.vrep:
			doc.removeCustomDrawableEntity(self.vrep)
		if self.vrep_glow:
			doc.removeCustomDrawableEntity(self.vrep_glow)
		App.updateActiveView()
	def draw(self, doc, cutelements):
		self.clear(doc)
		edge_glow = FxShapeEdge()
		edge_glow.edgeType = FxShapeEdgeType.Lines
		edge = FxShapeEdge()
		edge.edgeType = FxShapeEdgeType.Lines
		self.vrep_glow = FxShape()
		self.vrep_glow.edges.append(edge_glow)
		self.vrep = FxShape()
		self.vrep.edges.append(edge)
		# upoints: used to merge points from all separated segments
		upoints_counter = 1
		upoints = defaultdict(_TPID)
		# uedges: used to track the counter of edges (> 1 = internal)
		uedges = defaultdict(_TPID)
		for _, cele in cutelements.items():
			ncp = len(cele.cpoints)
			uids = [] # remap cpoints to unique ids (this vector should have len == len(cpoints)
			for ip in cele.cpoints:
				index = upoints[_TP(ip)]
				if index.i == 0:
					index.i = upoints_counter
					upoints_counter += 1
				uids.append(index.i)
			uids_ucounter = len(set(uids)) # count the number of unique ids (due to merging)
			if ncp == 1:
				# directly append vertices
				pass # todo: line -> vertex
			elif ncp > 1:
				if ncp == 2: # cut on surfaces
					min_ucounter = 1
					edge_generator = zip(range(ncp-1), range(1,ncp))
				else:
					min_ucounter = 2
					edge_generator = zip(range(ncp), itertools.chain(range(1,ncp), [0]))
				if uids_ucounter > min_ucounter: # this may happen due to rounding
					for ii,jj in edge_generator:
						i = uids[ii]
						j = uids[jj]
						if i != j: # this may happen due to rounding
							uedges[(j,i) if j < i else (i,j)].i += 1
		# keep unique edges
		uedges = [i for i,j in uedges.items() if j.i == 1]
		# inverse map, index in value, 1-based, due to dict order-by-insertion they are continuous
		upoints_inv = {b.i:a.p for a,b in upoints.items()}
		# save current offset in vrep edges
		edge_offset = len(edge.vertices) # same for glow edge
		# add unique points
		for _, ip in upoints_inv.items():
			for ie in (edge, edge_glow):
				ie.vertices.append(Math.vertex(ip))
		# add unique edges using offset - 1
		for ue in uedges:
			for i in ue:
				for ie in (edge, edge_glow):
					ie.indices.append(i-1+edge_offset)
		# glow material
		self.vrep.material = _CutElementsGraphics._base_color()
		self.vrep_glow.material = _CutElementsGraphics._glow_color()
		doc.addCustomDrawableEntity(self.vrep_glow)
		doc.addCustomDrawableEntity(self.vrep)
		App.updateActiveView()

class _PlotScaleBackup:
	def __init__(self):
		self.scales = {}
		self.begin()
	def begin(self):
		doc = App.postDocument()
		if doc is None: return
		if doc.activePlotGroup is None: return
		for pid, plot in doc.activePlotGroup.plots.items():
			self.scales[pid] = plot.plotData.scale
			plot.plotData.scale = 0.0
		doc.activePlotGroup.updateAndWait(MpcOdpRegeneratorUpdateInfo())
		App.setBusy(True)
	def end(self):
		doc = App.postDocument()
		if doc is None: return
		if doc.activePlotGroup is None: return
		if len(self.scales) == 0: return
		for pid, plot in doc.activePlotGroup.plots.items():
			plot.plotData.scale = self.scales.get(pid, plot.plotData.scale)
		doc.activePlotGroup.updateAndWait(MpcOdpRegeneratorUpdateInfo())
		App.setBusy(True)

class _MovingPointTracer(QObject):
	def __init__(self, doc, parent=None):
		super().__init__(parent=parent)
		self.doc = doc
		self.timer = QTimer()
		self.timer.setInterval(10)
		self.timer.timeout.connect(self.onTimeout)
		self.vrep = None
	def begin(self):
		self.timer.start()
	def end(self):
		self.timer.stop()
		if self.vrep:
			self.doc.removeCustomDrawableEntity(self.vrep)
	@Slot()
	def onTimeout(self):
		# todo: improve graphics
		# todo: edit point in space constraints?
		try:
			if self.vrep:
				self.doc.removeCustomDrawableEntity(self.vrep)
			self.vrep = FxShape()
			p = self.doc.scene.pointInSpace()
			self.vrep.vertices.vertices.append(Math.vertex(p))
			self.vrep.vertices.indices.append(0)
			mat = FxMaterial()
			mat.lighting = False
			mat.depthTest = False
			mat.transparency = 0.5
			mat.transparencyOnlyOnFaces = False
			mat.pointColor = FxColor(0.1, 0.8, 1.0)
			mat.pointSize = 10
			mat.visibilityOptionsOverride = FxMaterialVisibilityOptions(True, True, True, True)
			self.vrep.material = mat
			self.vrep.commitChanges()
			self.doc.addCustomDrawableEntity(self.vrep)
			App.updateActiveView()
		except Exception as ex:
			print(ex)

class _CutElement:
	# edges for different mesh types
	# assumptions: quadratic elements treated as linear: use only corner nodes
	_L2 = [(0,1)] # line
	_T3 = [(0,1),(1,2),(2,0)] # triangle
	_Q4 = [(0,1),(1,2),(2,3),(3,0)] # quadrilateral
	_T4 = [
		(0,1),(1,2),(2,0),
		(0,3),(1,3),(2,3)] # tetrahedron
	_H8 = [
		(0,1),(1,2),(2,3),(3,0),
		(4,5),(5,6),(6,7),(7,4),
		(0,4),(1,5),(2,6),(3,7)] # hexaedron
	def __init__(self, ele, cplane, db):
		self.ele = ele # the souce element
		self.csize = 0.0 # the cut size (1 for L2, length for T3 Q4, area for T4 H8)
		self.cpoints = [] # polygon points
		self.centroid = None # cut centroid
		self.cplane = cplane # ref to cut plane
		self.db = db # ref to database
		self.valid = False
		self._determineCut()
	def _processPolygon(self):
		# order cpoints to form a ccw polygon on the cut-plane
		# calculate cpoints in local coordinates (cut-normal = X)
		dx = self.cplane.dy
		dy = self.cplane.dz
		dz = self.cplane.dx
		p = [[ip.dot(dx), ip.dot(dy)] for ip in [jp - self.cplane.center for jp in self.cpoints]]
		# sort them to be ccw (use a temporary centroid) (, ip.dot(dz) is useless)
		n = len(p)
		cx = 0.0
		cy = 0.0
		for ip in p:
			cx += ip[0]
			cy += ip[1]
		cx /= n
		cy /= n
		p.sort(key = lambda ip : math.atan2(ip[1] - cy, ip[0] - cx))
		# compute area triangulation and centroid
		x0,y0 = p[0]
		A = 0.0
		cx = 0.0
		cy = 0.0
		for i in range(1, n-1):
			j = i+1
			pi = p[i]
			pj = p[j]
			iA = (-x0 + pi[0])*(-y0 + pj[1]) - (-x0 + pj[0])*(-y0 + pi[1])
			cx += (x0+pi[0]+pj[0])/3.0*iA
			cy += (y0+pi[1]+pj[1])/3.0*iA
			A += iA
		cx /= A
		cy /= A
		A /= 2.0
		# save data and come back to global coordinates
		R = Math.mat3(dx, dy, dz)
		self.csize = A
		self.centroid = R*Math.vec3(cx, cy, 0.0) + self.cplane.center
		self.cpoints = [R*Math.vec3(ip[0], ip[1], 0.0) + self.cplane.center for ip in p]
	def _determineCut(self):
		# references to cut-plane
		X = self.cplane.dx
		C = self.cplane.center
		T = self.cplane.tolerance
		if not self.cplane.infinite:
			Y = self.cplane.dy
			Z = self.cplane.dz
			HY = self.cplane.sizey/2.0
			HZ = self.cplane.sizez/2.0
		# compute signed distances for each point
		distances = [(node.position - C).dot(X)+T for node in self.ele.nodes]
		# get edges based on element type
		fam = self.ele.geometryFamilyType()
		geom_dim = 0
		if fam == MpcElementGeometryFamilyType.Line:
			edges = _CutElement._L2
			geom_dim = 1
		elif fam == MpcElementGeometryFamilyType.Triangle:
			edges = _CutElement._T3
			geom_dim = 2
		elif fam == MpcElementGeometryFamilyType.Quadrilateral:
			edges = _CutElement._Q4
			geom_dim = 2
		elif fam == MpcElementGeometryFamilyType.Tetrahedron:
			edges = _CutElement._T4
			geom_dim = 3
		elif fam == MpcElementGeometryFamilyType.Hexahedron:
			edges = _CutElement._H8
			geom_dim = 3
		else:
			# invalid
			return None
		# process each edge, if edge node have distances
		# with different signes there is an intersection
		for i,j in edges:
			di = distances[i]
			dj = distances[j]
			if di*dj < 0.0:
				alpha = di/(di-dj)
				pi = self.ele.nodes[i].position
				pj = self.ele.nodes[j].position
				palpha = pj*alpha + pi*(1.0-alpha)
				if not self.cplane.infinite:
					d_palpha = palpha - C
					dy, dysign = _asign(d_palpha.dot(Y))
					dz, dzsign = _asign(d_palpha.dot(Z))
					if dy > HY: palpha -= dysign*(dy-HY)*Y
					if dz > HZ: palpha -= dzsign*(dz-HZ)*Z
				self.cpoints.append(palpha)
		# verify cut
		if geom_dim == 1 and len(self.cpoints) == 1:
			self.centroid = self.cpoints[0]
			self.csize = 1.0
			self.valid = True
		elif geom_dim == 2 and len(self.cpoints) == 2:
			self.centroid = (self.cpoints[1] + self.cpoints[0]) / 2.0
			self.csize = (self.cpoints[1] - self.cpoints[0]).norm()
			if self.csize > 10.0*T: self.valid = True # skip small cuts
		elif geom_dim == 3 and len(self.cpoints) > 2:
			self._processPolygon()
			if math.sqrt(self.csize) > 10.0*T: self.valid = True # skip small cuts
	def _getTractionSolid(self, results):
		'''
		if extrapolated to nodes: interpolate at point cuts, 
		sum F at element cut centroid
		sum M local about element cut centroid
		'''
		if results.tensor is None:
			return (None, None)
		s = None
		wsum = 0.0
		for dpid, gp in enumerate(self.ele.integrationRule.integrationPoints):
			row = MpcOdbResultField.gauss(self.ele.id, dpid)
			w = gp.w
			try:
				if s is None:
					s = results.tensor_field[row] * w
				else:
					s += results.tensor_field[row] * w
				wsum += w
			except:
				return (None, None) # no result available
		if wsum == 0.0:
			return (None, None)
		s /= wsum
		cplane = self.cplane
		n = cplane.dx
		# traction vector in global coordinates
		T = Math.vec3(
			n.x*s[0] + n.y*s[3] + n.z*s[5],
			n.x*s[3] + n.y*s[1] + n.z*s[4],
			n.x*s[5] + n.y*s[4] + n.z*s[2])
		# lumped traction vector in cut-plane local coordinates
		F = Math.vec3(T.dot(cplane.dx)*self.csize, T.dot(cplane.dy)*self.csize, T.dot(cplane.dz)*self.csize)
		# distance from cut-center in local coordinates
		DP = self.centroid - cplane.center
		LP = Math.vec3(DP.dot(cplane.dx), DP.dot(cplane.dy), DP.dot(cplane.dz))
		# compute lumped moments due to lumped force eccentricity in cut-plane local coordinates (LP.X is zero in local coords!)
		M = Math.vec3(LP.y*F.z - LP.z*F.y, LP.z*F.x, -LP.y*F.x)
		# done
		return (F, M)
	def _getTractionShell(self, results):
		# valid for: 3D shells: 8 components:
		# Fxx Fyy Fxy Mxx Myy Mxy Vxz Vyz
		if results.shell is None:
			return (None, None)
		sf = None
		wsum = 0.0
		for dpid, gp in enumerate(self.ele.integrationRule.integrationPoints):
			row = MpcOdbResultField.gauss(self.ele.id, dpid)
			w = gp.w
			try:
				if sf is None:
					sf = results.shell_field[row] * w
				else:
					sf += results.shell_field[row] * w
				wsum += w
			except:
				return (None, None) # no result available
		if wsum == 0.0:
			return (None, None)
		sf /= wsum
		# cplane
		cplane = self.cplane
		# element orientation and its inverse
		Re = self.ele.orientation.computeOrientation()
		ReT = Re.transpose()
		# closest local axes aligned with plane: project cutplane on shell plane
		# explanation: if we use the real normal we need to compute the scale factor on the cut face
		# which will be (in general) different from the thicness (normal to the shell plane)
		ele_normal = Re.col(2)
		n = (cplane.dx - ele_normal*cplane.dx.dot(ele_normal)).normalized()
		# express it in element local coordinates
		n = ReT*n
		# extract membrane part + oop shear (constant part) in element local coordinates
		s = Math.vec(6)
		s[0] = sf[0] # Fxx
		s[1] = sf[1] # Fyy
		s[2] = 0.0   # Fzz = 0
		s[3] = sf[2] # Fxy
		s[4] = sf[7] # Fyz = |s23
		s[5] = sf[6] # Fxz = |s13
		# traction vector in element local coordinates
		T = Math.vec3(
			n.x*s[0] + n.y*s[3] + n.z*s[5],
			n.x*s[3] + n.y*s[1] + n.z*s[4],
			n.x*s[5] + n.y*s[4] + n.z*s[2])
		# traction vector in global coordinates
		T = Re*T
		# lumped traction vector in cut-plane local coordinates
		F = Math.vec3(T.dot(cplane.dx)*self.csize, T.dot(cplane.dy)*self.csize, T.dot(cplane.dz)*self.csize)
		# extract bending part in element local coordinates
		s[0] = sf[3] # Mxx
		s[1] = sf[4] # Myy
		s[2] = 0.0   # Mzz = 0
		s[3] = sf[5] # Mxy
		s[4] = 0.0   # Myz = 0
		s[5] = 0.0   # Mxz = 0
		# traction couple vector in element local coordinates
		T = Math.vec3(
			n.x*s[0] + n.y*s[3] + n.z*s[5],
			n.x*s[3] + n.y*s[1] + n.z*s[4],
			n.x*s[5] + n.y*s[4] + n.z*s[2])
		# swap first and second column:
		# assume Mxx = Sxx*H, Myy = Syy*H, Mxy = Sxy*H (2 equal and opposite stress tensor at dist = H)
		# so that S1 = [Sxx,Syy,Sxy] and S2 = -S1
		# T1 = S1*n and T2=S2*n = -S1*n
		# from traction vector to local couple vector about the shell mid-plane
		# d = [0 -H/2 0; H/2 0 0; 0 0 0]
		# M = d*S1 - d*S2 = d*S1 + d*S1 = S1*H swapping first and second column, changing size of second col.
		# but S1*H is M... so we can perform the same operations we did on the memebrane forces, just doing the swap
		T[0],T[1] = -T[1],T[0]
		# traction couple vector in global coordinates
		T = Re*T
		# lumped traction couple vector in cut-plane local coordinates
		M = Math.vec3(T.dot(cplane.dx)*self.csize, T.dot(cplane.dy)*self.csize, T.dot(cplane.dz)*self.csize)
		# distance from cut-center in local coordinates
		DP = self.centroid - cplane.center
		LP = Math.vec3(DP.dot(cplane.dx), DP.dot(cplane.dy), DP.dot(cplane.dz))
		# add lumped moments due to lumped force eccentricity in cut-plane local coordinates (LP.X is zero in local coords!)
		M += Math.vec3(LP.y*F.z - LP.z*F.y, LP.z*F.x, -LP.y*F.x)
		# done
		return (F, M)
	def _getTractionBeam(self, results):
		...
	def getTraction(self, results):
		family = self.ele.geometryFamilyType()
		dimension = self.db.info.spatialDimension
		if dimension == MpcOdbSpatialDimension.D1:
			return (None, None)
		if family == MpcElementGeometryFamilyType.Line:
			return self._getTractionBeam(results)
		elif family == MpcElementGeometryFamilyType.Triangle or family == MpcElementGeometryFamilyType.Quadrilateral:
			if dimension == MpcOdbSpatialDimension.D2:
				return self._getTractionSolid(results)
			else:
				return self._getTractionShell(results)
		elif family == MpcElementGeometryFamilyType.Tetrahedron or family == MpcElementGeometryFamilyType.Hexahedron:
			return self._getTractionSolid(results)
		else:
			return (None, None)

class _ResultCollection:
	def __init__(self):
		# tensorial result for both 2D and 3D models
		# like material.strain or material.stress
		self.tensor = None
		self.tensor_field = None
		# shell result for 3D models
		# like section.force
		self.shell = None
		self.shell_field = None
		# beam result, end-node forces in global coordinates
		self.beam = None
		self.beam_field = None
		# todo: trusses?
	def evaluate(self, opt):
		try:
			orient = opt.orientation
			if self.tensor and self.tensor.size() > 0:
				opt.orientation = MpcOdbResultOrientationType.Global
				self.tensor_field = self.tensor.evaluate(opt)
				if self.tensor_field is None:
					raise Exception('Cannot evaluate tensor field')
			if self.shell and self.shell.size() > 0:
				opt.orientation = MpcOdbResultOrientationType.Local
				self.shell_field = self.shell.evaluate(opt)
				if self.shell_field is None:
					raise Exception('Cannot evaluate shell field')
			if self.beam and self.beam.size() > 0:
				opt.orientation = MpcOdbResultOrientationType.Global
				self.beam_field = self.beam.evaluate(opt)
				if self.beam_field is None:
					raise Exception('Cannot evaluate beam field')
		finally:
			opt.orientation = orient

def _getDatabase(doc):
	print('1. DATABASE: Select the source database')
	data = {}
	db = None
	for id, db in doc.databases.items():
		key = '[{}] {}'.format(id, os.path.basename(db.fileName))
		data[key] = id
	if len(data) == 1:
		db_id = list(data.values())[0]
		db = doc.getDatabase(db_id)
	elif len(data) > 1:
		items = list(data.values())
		result, ok = QInputDialog.getItem(
			QApplication.activeWindow(), QApplication.applicationName(), 
			'Select DB', items, editable=False)
		if ok:
			db_id = data[result]
			db = doc.getDatabase(db_id)
	if db:
		print('   Selected database = "{}"\n   in "{}"'.format(
			os.path.basename(db.fileName),
			os.path.dirname(db.fileName)
			))
	else:
		print('   No database selected')
	return db

def _getMesh(db):
	print('2. MESH: Obtaining base mesh')
	U = db.getNodalResult('Displacement')
	if U is None:
		print('   Error: cannot find "Displacement" nodal result')
		return None
	all_stages = db.getStageIDs()
	if len(all_stages) == 0:
		print('   Error: "Displacement" nodal result has not stage')
		return None
	last_stage = all_stages[-1]
	all_steps = db.getStepIDs(last_stage)
	if len(all_steps) == 0:
		print('   Error: "Displacement" nodal result has not step')
		return None
	last_step = all_steps[-1]
	opt = MpcOdbVirtualResultEvaluationOptions()
	opt.stage = last_stage
	opt.step = last_step
	field = U.evaluate(opt)
	if field is None:
		print('   Error: Field from "Displacement" cannot be obtained')
		return None
	print('   Found a base mesh with {} nodes and {} elements'.format(
		len(field.mesh.nodes), len(field.mesh.elements)))
	return field.mesh

def _reduceMeshFromSelection(doc, mesh):
	print('3. MESH: Reducing mesh from selection')
	print('   Please perform a selection to consider only a portion of the whole mesh,\n'
		'   Press ENTER [or Right-Click] to accept the selection,\n'
		'   or press ESC to abort and use the whole mesh')
	code, value = doc.scene.waitForInput(
		mb = [int(Qt.RightButton)], 
		mk = [int(Qt.Key_Escape), int(Qt.Key_Enter), int(Qt.Key_Return)],
		timeout = _TIMEOUT,
		)
	if (code == MpcSceneKeyCode.Timeout) or (code == MpcSceneKeyCode.Key and value == int(Qt.Key_Escape)):
		print('   Result: using the whole mesh')
		return mesh
	t0 = time.time()
	source_nodes = dict((key,value) for key,value in mesh.nodes.items())
	source_elements = dict((key,value) for key,value in mesh.elements.items())
	mesh.nodes.clear()
	mesh.elements.clear()
	for plot_id, data in doc.scene.plotSelection.items():
		for id in data.info.elements:
			element = source_elements[id]
			mesh.addElement(element)
			for node in element.nodes:
				mesh.addNode(node)
	print('   Result: using reduced mesh with {} nodes and {} elements'.format(
		len(mesh.nodes), len(mesh.elements)))
	t1 = time.time()
	print('      Elapsed: {:8.4f} s. (Mesh reduce - selection)'.format(t1-t0))
	return mesh

def _reduceMeshFromIntersection(doc, db, mesh, cplane):
	t0 = time.time()
	source_elements = dict((key,value) for key,value in mesh.elements.items())
	mesh.nodes.clear()
	mesh.elements.clear()
	for _, ele in source_elements.items():
		dx_min = float('inf')
		dx_max = -dx_min
		dy_min = dx_min
		dz_min = dx_min
		for node in ele.nodes:
			# Test distance in the normal (X) direction
			# use a tolerance to keep elements on the -X side when cut is on edge!
			dpos = node.position - cplane.center
			dist = dpos.dot(cplane.dx)+cplane.tolerance
			dx_max = max(dx_max, dist)
			dx_min = min(dx_min, dist)
			size_check = True
			if not cplane.infinite:
				dy = abs(dpos.dot(cplane.dy))
				dz = abs(dpos.dot(cplane.dz))
				dy_min = min(dy_min, dy)
				dz_min = min(dz_min, dz)
				size_check = dy_min <= cplane.sizey/2.0 and dz_min <= cplane.sizez/2.0
			if dx_min < 0.0 and dx_max > 0.0 and size_check:
				# the element crosses the plane
				mesh.addElement(ele)
				for inode in ele.nodes:
					mesh.addNode(inode)
				break
	t1 = time.time()
	print('      Elapsed: {:8.4f} s. (Mesh reduce - intersection)'.format(t1-t0))
	return mesh

def _intersect(doc, db, mesh, cplane):
	t0 = time.time()
	source_elements = dict((key,value) for key,value in mesh.elements.items())
	mesh.nodes.clear()
	mesh.elements.clear()
	cutelements = {}
	for _, ele in source_elements.items():
		cutele = _CutElement(ele, cplane, db)
		if cutele.valid:
			cutelements[ele.id] = cutele
			mesh.addElement(ele)
			for inode in ele.nodes:
				mesh.addNode(inode)
	t1 = time.time()
	print('      Elapsed: {:8.4f} s. (Mesh cut)'.format(t1-t0))
	return (cutelements, mesh)

def _makeSectionCut():
	# start
	doc = App.postDocument()
	doc.clearCustomDrawableEntities()
	App.updateActiveView()
	# get database
	db = _getDatabase(doc)
	if db is None:
		QMessageBox.critical(_WIN, _TITLE, 'Abort:\nNo Database selected')
		return False
	# get the whole mesh
	mesh = _getMesh(db)
	if mesh is None:
		QMessageBox.critical(_WIN, _TITLE, 'Abort:\nCannot obtain the mesh')
		return False
	# reduce mesh using only selected elements
	mesh = _reduceMeshFromSelection(doc, mesh)
	doc.scene.unselectAll()
	selected_mesh_graphics = _SelectedMeshGraphics()
	selected_mesh_graphics.draw(doc, mesh)
	# get cut plane
	dia = _CPlaneDialog(mesh, doc=doc)
	result = dia.start()
	if result != QDialog.Accepted:
		QMessageBox.critical(_WIN, _TITLE, 'Abort:\nUser aborted cut-plane')
		return False
	cplane = dia.cplane
	# further reduction before intersection
	mesh = _reduceMeshFromIntersection(doc, db, mesh, cplane)
	selected_mesh_graphics.draw(doc, mesh)
	# perform intersections
	cutelements, mesh = _intersect(doc, db, mesh, cplane)
	selected_mesh_graphics.clear(doc)
	cut_graphics = _CutElementsGraphics()
	cut_graphics.draw(doc, cutelements)
	# get results
	# todo user GUI <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	results = _ResultCollection()
	#results.tensor = db.getElementalResult('material.stress', match=MpcOdb.StartsWith)
	results.tensor = db.getElementalResult('stress', match=MpcOdb.StartsWith)
	results.shell = db.getElementalResult('section.force (Surfaces; 8 Components', match=MpcOdb.StartsWith)
	# evaluate all results
	# todo ON ALL STEPS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	opt = MpcOdbVirtualResultEvaluationOptions()
	opt.mesh = mesh
	opt.stage = db.getStageIDs()[-1]
	opt.step = db.getStepIDs(opt.stage)[-1]
	results.evaluate(opt)
	# integrate
	F = Math.vec3(0.0, 0.0, 0.0)
	M = Math.vec3(0.0, 0.0, 0.0)
	S = []
	for _, ele in cutelements.items():
		iF, iM = ele.getTraction(results)
		if iF: F += iF
		if iM: M += iM
		# vrep
		if iF:
			S.append((ele.centroid, iF/ele.csize))
	# VREP <<<<<<<<<<<<<<
	T = Math.mat4()
	T[3,3] = 1.0
	dx_arrow = FxShapeFactory.arrow(Math.vec3(0.0,0.0,0.0), cplane.dx, 1.0, 0.02, 6, False)
	Fmax = 0.0
	for ip, iF in S:
		Fmax = max(Fmax, abs(iF.x))
	if Fmax > 0:
		scale = ((cplane.sizey + cplane.sizez)/2.0 * 0.1) / Fmax
		face = FxShapeFace()
		for ip, iF in S:
			iscale = scale * iF.x
			for i in range(3):
				T[i,i] = iscale
				T[i,3] = ip[i]
			icopy = FxShapeFace(dx_arrow)
			icopy.transform(T)
			face.merge(icopy)
		vrep = FxShape()
		vrep.faces.append(face)
		vrep.material = FxMaterial()
		vrep.material.diffuse = FxColor(1,0,0)
		vrep.material.depthTest = False
		vrep.material.transparency = 0.6
		vrep.material.transparencyOnlyOnFaces = False
		doc.addCustomDrawableEntity(vrep)
		App.updateActiveView()
	# VREP <<<<<<<<<<<<<<
	# done
	print(''.join(['{:12.3g}']*6).format(*F, *M))

# main call
try:
	App.clearTerminal()
	pscale_backup = _PlotScaleBackup()
	App.setBusy(True)
	_makeSectionCut()
finally:
	pscale_backup.end()
	App.setBusy(False)
	