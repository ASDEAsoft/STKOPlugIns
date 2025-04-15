from tclimport.mathutils import *

# some global settings
class settings:
	# supported objects
	BEAMS = ['dispBeamColumn', 'forceBeamColumn',       'elasticBeamColumn', 'ElasticTimoshenkoBeam']
	BEAMS_WITH_SEC = ['dispBeamColumn', 'forceBeamColumn']
	BEAMS_ELA = ['elasticBeamColumn', 'ElasticTimoshenkoBeam']
	SHELLS = ['ShellDKGQ', 'ShellNLDKGQ', 'ShellMITC4', 'ASDShellQ4']
	WALLS = ['MVLEM_3D']
	LINKS = ['zeroLength', 'zeroLengthSection', 'twoNodeLink', 'TripleFrictionPendulum']
	LINKS_WITH_SEC = ['zeroLengthSection']

# node
class node_t:
	def __init__(self, id, pos, mass):
		self.id = id
		self.pos = pos
		self.mass = mass
		self.cae_source = None # contains the (geom, subshape_id, and type)
	def __str__(self):
		return '[{}] {} -mass {}; CAE: {}'.format(self.id, self.pos, self.mass, self.cae_source)

# limitCurve
class limit_curve_t:
	def __init__(self, id, name, params):
		self.id = id
		self.name = name
		self.params = params
		self.original_id = id
	def __str__(self):
		return '[{}] - {} {} ({})'.format(self.id, self.name, ' '.join(str(i) for i in self.params), self.original_id)

# friction models
class friction_model_t:
	def __init__(self, id, name, params):
		self.id = id
		self.name = name
		self.params = params
		self.original_id = id
	def __str__(self):
		return '[{}] - {} {} ({})'.format(self.id, self.name, ' '.join(str(i) for i in self.params), self.original_id)

# region
class region_t:
	def __init__(self, id, the_nodes = None, the_elements = None):
		self.id = id
		self.nodes = the_nodes
		self.elements = the_elements
		self.cae_source = None
	def __str__(self):
		if self.nodes is not None:
			return '[{}] -N {} - CAE: {}'.format(self.id, self.nodes, self.cae_source)
		else:
			return '[{}] -E {} - CAE: {}'.format(self.id, self.elements, self.cae_source)

# recorder
class recorder_t:
	def __init__(self, pre_data, region, post_data):
		self.pre_data = pre_data
		self.region = region
		self.post_data = post_data
	def __str__(self):
		return 'recorder {} -region {} {}'.format(
			' '.join(self.pre_data),
			self.region,
			' '.join(self.post_data))

# trans
class geom_trans_t:
	def __init__(self, id, type, vz, offi, offj):
		self.id = id
		self.type = type
		self.vz = vz
		self.vx = None
		self.vy = None
		self.offi = offi
		self.offj = offj
		self.cae_source = None # contains the MpcLocalAxes if not equal to the STKO default
	def setxy(self, vx, vy):
		# for links we use vx and vy instead of vz
		vz = None
		self.vx = vx
		self.vy = vy
	def __str__(self):
		if self.vz:
			return '[{}] type: {}; vz: {}; offi: {}; offj: {}'.format(self.id, self.type, self.vz, self.offi, self.offj)
		else:
			return '[{}] type: {}; vx: {}; vy: {}'.format(self.id, self.type, self.vx, self.vy)

# element
class element_t:
	UNKNOWN = 0
	EDGE = 1
	FACE = 2
	LINK = 3
	def __init__(self, id, nodes, type, name, params):
		self.id = id
		self.nodes = nodes
		self.type = type
		self.name = name
		self.params = params
		self.cae_source = None # contains the (geom, subshape_id, and type) if type != LINK, otherwise the interaction
	def __str__(self):
		return '[{}] nodes: {}; type: {}; name: {}; params: {}; CAE: {}'.format(self.id, [n.id for n in self.nodes], self.type, self.name, self.params, self.cae_source)

# geometry
class geometry_t:
	def __init__(self, id):
		self.id = id
		self.elements = []
		self.type = element_t.UNKNOWN
	def computeType(self):
		self.type = element_t.UNKNOWN
		itype = None
		for e in self.elements:
			if itype is None:
				itype = e.type
			else:
				if itype != e.type:
					return
		if itype is not None:
			self.type = itype
	def __str__(self):
		return 'Geometry {} (Type: {}) {{\n{}\n}}\n'.format(
			self.id, 
			self.type, 
			'\n'.join(['\t{}'.format(i) for i in self.elements]))

# a fiber path
class fiber_patch_t:
	def __init__(self, mat, nx, ny, p1, p2, p3, p4):
		self.mat = mat
		self.nx = nx
		self.ny = ny
		self.p = (p1, p2, p3, p4)
	def copy(self):
		return fiber_patch_t(self.mat, self.nx, self.ny, self.p[0], self.p[1], self.p[2], self.p[3])
	def __str__(self):
		return '{} - {} {} {} {} - {} x {}'.format(self.mat, *self.p, self.nx, self.ny)
	def __repr__(self):
		return self.__str__()

# a fiber layer
class fiber_layer_t:
	def __init__(self, mat, num, area, p1, p2):
		self.mat = mat
		self.num = num
		self.area = area
		self.p = (p1, p2)
	def copy(self):
		return fiber_layer_t(self.mat, self.num, self.area, self.p[0], self.p[1])
	def __str__(self):
		return '{} - {} @ {}'.format(self.mat, self.num, self.area)
	def __repr__(self):
		return self.__str__()

# fiber point
class fiber_point_t:
	def __init__(self, mat, y, z, area):
		self.mat = mat
		self.y = y
		self.z = z
		self.a = area
	def copy(self):
		return fiber_point_t(self.mat, self.y, self.z, self.a)
	def __str__(self):
		return '{} - ({}, {}) A = {}'.format(self.mat, self.y, self.z, self.a)
	def __repr__(self):
		return self.__str__()

# fiber sections
class fibers_t:
	def __init__(self, id):
		self.id = id
		self.gj = 0.0
		self.patches = []
		self.layers = []
		self.points = []
		self.offy = 0.0
		self.offz = 0.0
	def copy(self, id):
		c = fibers_t(id)
		c.gj = self.gj
		c.offy = self.offy
		c.offz = self.offz
		c.patches = [i.copy() for i in self.patches]
		c.layers = [i.copy() for i in self.layers]
		c.points = [i.copy() for i in self.points]
		return c
	def __str__(self):
		return '[{}], GJ = {}, patches = {}, layers = {}, points = {}, offset = {}'.format(self.id, self.gj, self.patches, self.layers, self.points, (self.offy, self.offz))

# aggregator
class aggregator_t:
	def __init__(self, id):
		self.id = id
		self.materials = []
		self.section = None
	def __str__(self):
		return '[{}], materials = {}, section = {}'.format(self.id, self.materials, self.section)

# beam property
class beam_prop_t:
	def __init__(self, id):
		self.id = id
		self.section = None
		self.nump = 5
	def __str__(self):
		return '[{}], section = {}, num points = {}'.format(self.id, self.section, self.nump)

# plate section
class plate_section_t:
	def __init__(self, id, E, nu, h, rho, ep_mod):
		self.id = id
		self.E = E
		self.nu = nu
		self.h = h
		self.rho = rho
		self.ep_mod = ep_mod
	def __str__(self):
		return '[{}], E = {}, nu = {}, h = {}, rho = {}, Ep_mod = {}'.format(self.id, self.E, self.nu, self.h, self.rho, self.ep_mod)

# elastic beam section
class section_t:
	def __init__(self, id, params):
		self.id = id
		self.E = params[0]
		self.G = params[1]
		self.A = params[2]
		self.Jx = params[3]
		self.Jy = params[4]
		self.Jz = params[5]
		self.Avy = self.A*5.0/6.0
		self.Avz = self.A*5.0/6.0
		self.offy = 0.0
		self.offz = 0.0
		if len(params) == 8:
			self.Avy = params[6]
			self.Avz = params[7]
	def copy(self, id):
		c = section_t(id, [self.E,self.G,self.A,self.Jx,self.Jy,self.Jz])
		c.Avy = self.Avy
		c.Avz = self.Avz
		c.offy = self.offy
		c.offz = self.offz
		return c
	def __eq__(self, b):
		def cmp(x, y):
			if x is None and y is None: return True
			return abs(x-y)<1.0e-12
		return cmp(self.E, b.E) and cmp(self.G, b.G) and cmp(self.A, b.A) and cmp(self.Jx, b.Jx) and  cmp(self.Jy, b.Jy) and cmp(self.Jz, b.Jz) and cmp(self.Avy, b.Avy) and cmp(self.Avz, b.Avz)
	def __hash__(self):
		return hash((self.E, self.G, self.A, self.Jx, self.Jy, self.Jz, self.Avy, self.Avz))
	def __str__(self):
		return '[{}], {}, offset = {}'.format(self.id, [self.E, self.G, self.A, self.Jx, self.Jy, self.Jz, self.Avy, self.Avz], (self.offy, self.offz))

# uniaxial materials
class material_1d_t:
	def __init__(self, id, name, params):
		self.id = id
		self.name = name
		self.params = params
	def __str__(self):
		return '[{}], name = {}, params = {}'.format(self.id, self.name, self.params)

# link material
class link_material_t:
	ZLEN = 1
	LINK = 2
	TFP  = 3
	def __init__(self, id, mats, dirs):
		self.id = id
		self.mats = mats
		self.dirs = dirs
		self.type = link_material_t.ZLEN
	def __str__(self):
		return '[{}], mats = {}, dirs = {} ({})'.format(self.id, self.mats, self.dirs, self.type)

# fix
class fix_t:
	all_labels = ('Ux', 'Uy', 'Uz', 'Rx', 'Ry', 'Rz/3D')
	def __init__(self, dofs):
		self.nodes = []
		self.dofs = dofs # a int bitflag
	def dofstr(self):
		return ''.join([str(int(bool(self.dofs & (1 << j)))) for j in range(6)])
	def doflabels(self):
		return [self.all_labels[j] for j in range(6) if bool(self.dofs & (1 << j))]
	def __str__(self):
		return 'dofs: {}; nodes: {}'.format(
			self.dofstr(),
			self.nodes)

# diaphragm
class diaphragm_t:
	def __init__(self, dir, master):
		self.dir = dir
		self.master = master
		self.slaves = []
		self.cae_source = None
	def __str__(self):
		return 'dir: {}; master: {}; nodes: {}; CAE: {}'.format(self.dir, self.master, self.slaves, self.cae_source)

#edof_t
class edof_t:
	all_labels = ('Ux', 'Uy', 'Uz', 'Rx', 'Ry', 'Rz/3D')
	def __init__(self, master, dofs):
		self.master = master
		self.dofs = dofs # a int bitflag
		self.slaves = []
		self.cae_source = None
	def dofstr(self):
		return ''.join([str(int(bool(self.dofs & (1 << j)))) for j in range(6)])
	def doflabels(self):
		return [self.all_labels[j] for j in range(6) if bool(self.dofs & (1 << j))]
	def __str__(self):
		return 'dofs: {}; master: {}; slaves: {}; CAE: {}'.format(
			self.dofstr(),
			self.master, self.slaves, self.cae_source)

# load pattern
class pattern_t:
	def __init__(self, id):
		self.id = id
		self.forces = {}
		self.couples = {}
		self.ele_loads = {}
		self.__aux = {} #map node to load
	def add_load(self, load, node):
		if node in self.__aux:
			value = self.__aux[node]
		else:
			value = [0.0]*6
			self.__aux[node] = value
		for i in range(6):
			value[i] += load[i]
	def add_ele_load(self, load, ele):
		eles = self.ele_loads.get(load, None)
		if eles is None:
			eles = []
			self.ele_loads[load] = eles
		eles.append(ele)
	def finalize(self):
		for node, load in self.__aux.items():
			F = (load[0], load[1], load[2])
			M = (load[3], load[4], load[5])
			if v3norm(F) > 1.0e-12:
				if F in self.forces:
					nodes = self.forces[F]
				else:
					nodes = []
					self.forces[F] = nodes
				nodes.append(node)
			if v3norm(M) > 1.0e-12:
				if M in self.couples:
					nodes = self.couples[M]
				else:
					nodes = []
					self.couples[M] = nodes
				nodes.append(node)
	def __str__(self):
		return 'pattern: {}; forces: {}; couples: {}; ele_loads: {}'.format(self.id, self.forces, self.couples, self.ele_loads)

# stage
class stage_t:
	def __init__(self, id):
		self.id = id
		self.fix = {}
		self.diaphragms = {}
		self.edofs = {}
		self.patterns = {}
	# return a string representation of self
	def __str__(self):
		return (
			'[{}] STAGE {{\n'
			'\t\t\tFIX {{\n{}\n\t\t\t}}\n'
			'\t\t\tDIAPHRAGMS {{\n{}\n\t\t\t}}\n'
			'\t\t\tEQUAL_DOFS {{\n{}\n\t\t\t}}\n'
			'\t\t\tPATTERNS {{\n{}\n\t\t\t}}\n'
			'\t\t}}'
			).format(
				self.id,
				'\n'.join(['\t\t\t\t{}'.format(i) for i in self.fix.values()]),
				'\n'.join(['\t\t\t\t{}'.format(i) for i in self.diaphragms.values()]),
				'\n'.join(['\t\t\t\t{}'.format(i) for i in self.edofs.values()]),
				'\n'.join(['\t\t\t\t{}'.format(i) for i in self.patterns.values()]),
				)