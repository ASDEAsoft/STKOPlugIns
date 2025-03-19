import importlib
import tclimport.tclobjects
import tclimport.graph
importlib.reload(tclimport.tclobjects)
importlib.reload(tclimport.graph)
from tclimport.tclobjects import *
from tclimport.graph import *
from PyMpc import *
import math

#############################################################################
# 
#############################################################################

class _support:
	def __init__(self):
		# the current step id
		self.step_id = 1
		# the current condition id
		self.condition_id = 1
		# the current stage to build
		self.stage = None

# makes a sorted tuple
def _sorted_pair(a, b):
	if a < b:
		return (a, b)
	else:
		return (b, a)

# as taken from the elasticBeamColumn XObject
class _utils_release:
	release_code_map = {1 : 'I-End only', 2 : 'J-End only', 3 : 'I-End & J-End'}

# build: group nodes using graph
def _build_group_nodes(doc):
	the_graph = Graph(len(doc.nodes))
	for ele_id, ele in doc.elements.items():
		for i in range((len(ele.nodes))-1):
			n1 = ele.nodes[i].id
			n2 = ele.nodes[i+1].id
			the_graph.addSub(n1, n2)
	return the_graph.connectedComponents()

# build: makes geometries grouping
# elements by connectivity
def _build_make_geoms(doc):
	node_groups = _build_group_nodes(doc)
	geometries = []
	geom_id = 0
	for node_ids in node_groups:
		# find all elements sharing at least one node with node_ids
		igeom = geometry_t(0)
		for ele_id, ele in doc.elements.items():
			for node in ele.nodes:
				if node.id in node_ids:
					igeom.elements.append(ele)
					break # break node loop
		if len(igeom.elements) > 0:
			geom_id += 1
			igeom.id = geom_id
			igeom.computeType()
			geometries.append(igeom)
	return geometries

# build: makes a wire from a geometry
def _build_make_wire(doc, geom, get_edge):
	source = geom.elements[:]
	added_nodes = []
	edges = []
	iter = 0
	while(True):
		iter += 1
		if iter > len(geom.elements)+100:
			raise Exception('make-wire iteration failed') # this should never happen
		#print('iter {} - {} / {}\n'.format(iter, len(source), len(geom.elements)))
		# process all elements in source.
		# if an element cannot be added, put it in source_ex
		source_aux = []
		for ele in source:
			n1 = ele.nodes[0].id
			n2 = ele.nodes[1].id
			if (len(edges) == 0) or (n1 in added_nodes) or (n2 in added_nodes):
				edge = get_edge(n1, n2)
				ele.cae_source = edge # temporary store the fxocc shape
				edges.append(get_edge(n1, n2))
				added_nodes.append(n1)
				added_nodes.append(n2)
			else:
				source_aux.append(ele)
		if len(source_aux) == len(source):
			break # nothing added in this trial
		if len(source_aux) == 0:
			break # all edges added
		source = source_aux
	# make the wire
	return FxOccBuilder.makeWire(edges)

# build: makes a shell from a geometry
def _build_make_shell(doc, geom, get_edge):
	# build all faces
	faces = []
	for ele in geom.elements:
		e1 = get_edge(ele.nodes[0].id, ele.nodes[1].id)
		e2 = get_edge(ele.nodes[1].id, ele.nodes[2].id)
		e3 = get_edge(ele.nodes[2].id, ele.nodes[3].id)
		e4 = get_edge(ele.nodes[3].id, ele.nodes[0].id)
		wire = FxOccBuilder.makeWire([e1, e2, e3, e4])
		face = FxOccBuilder.makeFace(wire)
		ele.cae_source = face # temporary store the fxocc shape
		faces.append(face)
	return FxOccBuilder.makeShell(faces)

# build: makes a generic compound
def _build_make_comp(doc, geom, get_edge):
	sources = []
	for ele in geom.elements:
		if ele.type == element_t.EDGE:
			edge = get_edge(ele.nodes[0].id, ele.nodes[1].id)
			ele.cae_source = edge # temporary store the fxocc shape
			sources.append(edge)
		elif ele.type == element_t.FACE:
			e1 = get_edge(ele.nodes[0].id, ele.nodes[1].id)
			e2 = get_edge(ele.nodes[1].id, ele.nodes[2].id)
			e3 = get_edge(ele.nodes[2].id, ele.nodes[3].id)
			e4 = get_edge(ele.nodes[3].id, ele.nodes[0].id)
			wire = FxOccBuilder.makeWire([e1, e2, e3, e4])
			face = FxOccBuilder.makeFace(wire)
			ele.cae_source = face # temporary store the fxocc shape
			sources.append(face)
	return FxOccBuilder.makeCompound(sources)

# build cae: model
def _build_cae_model(doc, callbacks):
	
	# document
	cae_doc = App.caeDocument()
	
	# some id counters
	interaction_id = 0
	selection_set_id = cae_doc.selectionSets.getlastkey(0)
	
	# re-number nodes in 0-based sequential indexing
	# original indices are in k
	counter = 0
	for k, v in doc.nodes.items():
		v.id = counter
		counter += 1
	
	# creates geometries as compound of elements
	geometries = _build_make_geoms(doc)
	
	# now we can restore original node ids
	for k, v in doc.nodes.items():
		v.id = k
	
	# make vertex map (maps a node_id to a vertex)
	# MAP<int, vertex>
	vertex_map = {}
	for node_id, node in doc.nodes.items():
		p = node.pos
		vertex_map[node.id] = FxOccBuilder.makeVertex(p[0], p[1], p[2])
	
	# make edge map (maps 2 node_ids to an edge)
	# MAP<tuple(2-sorted), edge>
	edge_map = {}
	def get_edge(n1, n2):
		sip = _sorted_pair(n1, n2)
		try:
			edge = edge_map[sip]
			return edge
		except:
			v1 = vertex_map[n1]
			v2 = vertex_map[n2]
			edge = FxOccBuilder.makeEdge(v1, v2)
			edge_map[sip] = edge
			return edge
	
	# build geometries from elements and add them to the document
	for geom in geometries:
		if len(geom.elements) == 0:
			continue
		if geom.type == element_t.EDGE:
			wire = _build_make_wire(doc, geom, get_edge)
			callbacks.addGeometry(MpcGeometry(geom.id, "Wire_{}".format(geom.id), wire))
		elif geom.type == element_t.FACE:
			shell = _build_make_shell(doc, geom, get_edge)
			callbacks.addGeometry(MpcGeometry(geom.id, "Shell_{}".format(geom.id), shell))
		else:
			comp = _build_make_comp(doc, geom, get_edge)
			callbacks.addGeometry(MpcGeometry(geom.id, "Compound_{}".format(geom.id), comp))
		# override the element_t.cae_source here that the geom has been added
		for ele in geom.elements:
			if ele.cae_source is not None:
				ele.cae_source = callbacks.findGeometry(ele.cae_source)[0]
	
	# build geometries for hanging nodes
	node_usage = {node_id:False for node_id in doc.nodes} # key = node_id, value = True if used
	for geom in geometries:
		for ele in geom.elements:
			for node in ele.nodes:
				node_usage[node.id] = True
	hnodes = []
	for node_id, used in node_usage.items():
		if not used:
			hnodes.append(doc.nodes[node_id])
	# make node
	geom_id = len(geometries) + 1
	for ihnode in hnodes:
		hvertex = vertex_map[ihnode.id]
		callbacks.addGeometry(MpcGeometry(geom_id, "Vertex_{}".format(geom_id), hvertex))
		geom_id += 1
	
	# compute cae_source of nodes
	for _, node in doc.nodes.items():
		node.cae_source = callbacks.findGeometry(vertex_map[node.id])[0]
	
	# link 2-node elements that become interactions
	for _, ele in doc.links.items():
		inter = MpcInteraction(ele.id, "Link_{}".format(ele.id))
		inter.type = MpcInteractionType.NodeToNode
		geom, subshape_id, subshape_type = callbacks.findGeometry(vertex_map[ele.nodes[0].id])[0] 
		inter.items.masters.append(MpcInteractionItem(geom, subshape_type, subshape_id))
		geom, subshape_id, subshape_type = callbacks.findGeometry(vertex_map[ele.nodes[1].id])[0]
		inter.items.slaves.append(MpcInteractionItem(geom, subshape_type, subshape_id))
		callbacks.addInteraction(inter)
		ele.cae_source = inter # for interaction we just need the interaction
		if ele.id > interaction_id:
			interaction_id = ele.id
	interaction_id += 1
	
	# conditions that become interactions
	for stage in doc.stages:
		# equal dofs
		for _, edof in stage.edofs.items():
			interaction_id += 1
			inter = MpcInteraction(interaction_id, "EqualDOF_{}".format(edof.dofstr()))
			inter.type = MpcInteractionType.NodeToNode
			geom, subshape_id, subshape_type = callbacks.findGeometry(vertex_map[edof.master])[0]
			inter.items.masters.append(MpcInteractionItem(geom, subshape_type, subshape_id))
			for slave in edof.slaves:
				geom, subshape_id, subshape_type = callbacks.findGeometry(vertex_map[slave])[0]
				inter.items.slaves.append(MpcInteractionItem(geom, subshape_type, subshape_id))
			callbacks.addInteraction(inter)
			edof.cae_source = inter
		# rigid diaphragms
		for _, dia in stage.diaphragms.items():
			interaction_id += 1
			inter = MpcInteraction(interaction_id, "RigidDiaphragm_[{}]({})".format(dia.master, dia.dir))
			inter.type = MpcInteractionType.NodeToNode
			geom, subshape_id, subshape_type = callbacks.findGeometry(vertex_map[dia.master])[0]
			inter.items.masters.append(MpcInteractionItem(geom, subshape_type, subshape_id))
			for slave in dia.slaves:
				geom, subshape_id, subshape_type = callbacks.findGeometry(vertex_map[slave])[0]
				inter.items.slaves.append(MpcInteractionItem(geom, subshape_type, subshape_id))
			callbacks.addInteraction(inter)
			dia.cae_source = inter
	
	# regions: make selection sets.
	# for each region we make an a selection set and a new region analysis step
	for region_id, region in doc.regions.items():
		selection_set_id += 1
		sset = MpcSelectionSet(selection_set_id, "Selection_Set_REGION_{}".format(region_id))
		if region.nodes is not None:
			for node_id in region.nodes:
				vertex = vertex_map[node_id]
				geom, subshape_id, subshape_type = callbacks.findGeometry(vertex)[0]
				if geom.id in sset.geometries:
					item = sset.geometries[geom.id]
				else:
					item = MpcSelectionSetItem()
					item.wholeGeometry = False
					sset.geometries[geom.id] = item
				item.vertices.append(subshape_id)
		else:
			for ele_id in region.elements:
				if ele_id in doc.elements:
					# geometry
					ele = doc.elements[ele_id]
					geom, subshape_id, subshape_type = ele.cae_source
					if geom.id in sset.geometries:
						item = sset.geometries[geom.id]
					else:
						item = MpcSelectionSetItem()
						item.wholeGeometry = False
						sset.geometries[geom.id] = item
					if ele.type == element_t.EDGE:
						item.edges.append(subshape_id)
					elif ele.type == element_t.FACE:
						item.faces.append(subshape_id)
				elif ele_id in doc.links:
					# interaction
					ele = doc.links[ele_id]
					sset.interactions.append(ele.cae_source.id)
		callbacks.addSelectionSet(sset)
		region.cae_source = sset
	
	# done
	return (vertex_map, edge_map)


# build cae: definitions
def _build_cae_definitions(doc, callbacks, vertex_map):
	
	cae_doc = App.caeDocument()
	
	def make_sset_node(node_id, curve_id, name):
		sset_id = cae_doc.selectionSets.getlastkey(0) + 1
		sset = MpcSelectionSet(sset_id, "Selection_Set_LimitCurve_{}_{}".format(curve_id, name))
		vertex = vertex_map[node_id]
		geom, subshape_id, subshape_type = callbacks.findGeometry(vertex)[0]
		if geom.id in sset.geometries:
			item = sset.geometries[geom.id]
		else:
			item = MpcSelectionSetItem()
			item.wholeGeometry = False
			sset.geometries[geom.id] = item
		item.vertices.append(subshape_id)
		callbacks.addSelectionSet(sset)
		return sset_id
	
	def make_sset_elem(ele_id, curve_id, name):
		sset_id = cae_doc.selectionSets.getlastkey(0) + 1
		sset = MpcSelectionSet(sset_id, "Selection_Set_LimitCurve_{}_{}".format(curve_id, name))
		if ele_id in doc.elements:
			# geometry
			ele = doc.elements[ele_id]
			geom, subshape_id, subshape_type = ele.cae_source
			if geom.id in sset.geometries:
				item = sset.geometries[geom.id]
			else:
				item = MpcSelectionSetItem()
				item.wholeGeometry = False
				sset.geometries[geom.id] = item
			if ele.type == element_t.EDGE:
				item.edges.append(subshape_id)
			elif ele.type == element_t.FACE:
				item.faces.append(subshape_id)
		elif ele_id in doc.links:
			# interaction
			ele = doc.links[ele_id]
			sset.interactions.append(ele.cae_source.id)
		callbacks.addSelectionSet(sset)
		return sset_id
	
	meta = callbacks.getMetaDefinition('limitCurves.Rotation')
	for _, icurve in doc.limit_curves.items():
		defi = MpcDefinition()
		defi.name = 'LimitCurve Rotation ({})'.format(icurve.original_id)
		defi.id = icurve.id
		xobj = MpcXObject.createInstanceOf(meta)
		xobj.getAttribute('eleTag').index = make_sset_elem(icurve.params[0], icurve.original_id, 'eleTag')
		xobj.getAttribute('dofl').integer = icurve.params[1]
		xobj.getAttribute('dofv').integer = icurve.params[2]
		xobj.getAttribute('iNodeTag').index = make_sset_node(icurve.params[3], icurve.original_id, 'iNodeTag')
		xobj.getAttribute('jNodeTag').index = make_sset_node(icurve.params[4], icurve.original_id, 'jNodeTag')
		xobj.getAttribute('fpc').quantityScalar.value = icurve.params[5]
		xobj.getAttribute('fyt').quantityScalar.value = icurve.params[6]
		xobj.getAttribute('Ag').quantityScalar.value = icurve.params[7]
		xobj.getAttribute('rhot').real = icurve.params[8]
		xobj.getAttribute('thetay').real = icurve.params[9]
		xobj.getAttribute('VColOE').quantityScalar.value = icurve.params[10]
		xobj.getAttribute('Kunload').quantityScalar.value = icurve.params[11]
		if icurve.params[12]:
			xobj.getAttribute('-VyE').boolean = True
			xobj.getAttribute('VyE').quantityScalar.value = icurve.params[12]
		defi.XObject = xobj
		callbacks.addDefinition(defi)
		defi.commitXObjectChanges()


# build cae: local axes
def _build_cae_locax(doc, callbacks):
	# add local axes
	for _, t in doc.geom_trans.items():
		callbacks.addLocalAxes(t)
	# assign local axes
	for _, ele in doc.elements.items():
		if ele.name in tclimport.tclobjects.settings.BEAMS:
			id = ele.params[-2]# the last is the transf type from end of pars (sec, releasecodes, tid, ttype)
			if id is not None:
				callbacks.select(ele.cae_source)
				callbacks.runCommand('AssignLocalAxes', str(id))
				callbacks.unselect()
	for _, ele in doc.links.items():
		id = ele.params[1]
		if id is not None:
			callbacks.select(ele.cae_source)
			callbacks.runCommand('AssignLocalAxes', str(id))
			callbacks.unselect()

# build cae: uniaxial
def _build_cae_pprop_uniaxial(doc, callbacks):
	for id, mat in doc.materials_1d.items():
		meta = callbacks.getMetaPhysicalProperty('materials.uniaxial.{}'.format(mat.name))
		xobj = MpcXObject.createInstanceOf(meta)
		if mat.name == 'Elastic' or mat.name == 'ENT':
			xobj.getAttribute('E').quantityScalar.value = mat.params[0]
		elif mat.name == 'Concrete01':
			xobj.getAttribute('fpc').quantityScalar.value = mat.params[0]
			xobj.getAttribute('epsc0').real = mat.params[1]
			xobj.getAttribute('fpcu').quantityScalar.value = mat.params[2]
			xobj.getAttribute('epscu').real = mat.params[3]
		elif mat.name == 'Concrete02':
			xobj.getAttribute('fpc').quantityScalar.value = mat.params[0]
			xobj.getAttribute('epsc0').real = mat.params[1]
			xobj.getAttribute('fpcu').quantityScalar.value = mat.params[2]
			xobj.getAttribute('epsU').real = mat.params[3]
			xobj.getAttribute('lambda').real = mat.params[4]
			xobj.getAttribute('ft').quantityScalar.value = mat.params[5]
			xobj.getAttribute('Ets').quantityScalar.value = mat.params[6]
		elif mat.name == 'Steel02':
			xobj.getAttribute('Fy').quantityScalar.value = mat.params[0]
			xobj.getAttribute('E0').quantityScalar.value = mat.params[1]
			xobj.getAttribute('b').real = mat.params[2]
			xobj.getAttribute('R0').real = mat.params[3]
			xobj.getAttribute('CR1').real = mat.params[4]
			xobj.getAttribute('CR2').real = mat.params[5]
		elif mat.name == 'Hysteretic':
			# (18) $matTag $s1p $e1p $s2p $e2p <$s3p $e3p> $s1n $e1n $s2n $e2n <$s3n $e3n> $pinchX $pinchY $damage1 $damage2 <$beta>
			# (17) $matTag $s1p $e1p $s2p $e2p <$s3p $e3p> $s1n $e1n $s2n $e2n <$s3n $e3n> $pinchX $pinchY $damage1 $damage2
			# (13) $matTag $s1p $e1p $s2p $e2p $s1n $e1n $s2n $e2n $pinchX $pinchY $damage1 $damage2
			n = len(mat.params)
			if n == 12:
				par_names = ('s1p', 'e1p', 's2p', 'e2p', 's1n', 'e1n', 's2n', 'e2n', 'pinchx', 'pinchy', 'damage1', 'damage2')
			elif n == 16:
				par_names = ('s1p', 'e1p', 's2p', 'e2p', 's3p', 'e3p', 's1n', 'e1n', 's2n', 'e2n', 's3n', 'e3n', 'pinchx', 'pinchy', 'damage1', 'damage2')
				xobj.getAttribute('Optional').boolean = True
			elif n == 17:
				par_names = ('s1p', 'e1p', 's2p', 'e2p', 's3p', 'e3p', 's1n', 'e1n', 's2n', 'e2n', 's3n', 'e3n', 'pinchx', 'pinchy', 'damage1', 'damage2', 'beta')
				xobj.getAttribute('Optional').boolean = True
				xobj.getAttribute('use_beta').boolean = True
			else:
				raise Exception('Hysteretic format not supported')
			for par, val in zip(par_names, mat.params):
				atr = xobj.getAttribute(par)
				if atr.type == MpcAttributeType.Real:
					atr.real = val
				else:
					atr.quantityScalar.value = val
		elif mat.name == 'MinMax':
			xobj.getAttribute('otherTag').index = mat.params[0]
			mmin = mat.params[1]
			mmax = mat.params[2]
			if mmin is not None:
				xobj.getAttribute('-min').boolean = True
				xobj.getAttribute('minStrain').real = mmin
			if mmax is not None:
				xobj.getAttribute('-max').boolean = True
				xobj.getAttribute('maxStrain').real = mmax
		elif mat.name == 'Pinching4':
			eP = mat.params[0]
			eN = mat.params[1]
			rP = mat.params[2]
			rN = mat.params[3]
			gK = mat.params[4]
			gD = mat.params[5]
			gF = mat.params[6]
			gE = mat.params[7]
			opt = eN[0] is not None
			dType = mat.params[8]
			ePlabels = ('ePf1', 'ePd1', 'ePf2', 'ePd2', 'ePf3', 'ePd3', 'ePf4', 'ePd4')
			eNlabels = ('eNf1', 'eNd1', 'eNf2', 'eNd2', 'eNf3', 'eNd3', 'eNf4', 'eNd4')
			rPlabels = ('rDispP', 'rForceP', 'uForceP')
			rNlabels = ('rDispN', 'rForceN', 'uForceN')
			gKlabels = ('gK1', 'gK2', 'gK3', 'gK4', 'gKLim')
			gDlabels = ('gD1', 'gD2', 'gD3', 'gD4', 'gDLim')
			gFlabels = ('gF1', 'gF2', 'gF3', 'gF4', 'gFLim')
			if opt:
				xobj.getAttribute('Optional').boolean = True
			for i in range(len(eP)):
				xobj.getAttribute(ePlabels[i]).real = eP[i]
			if opt:
				for i in range(len(eN)):
					xobj.getAttribute(eNlabels[i]).real = eN[i]
			for i in range(len(rP)):
				xobj.getAttribute(rPlabels[i]).real = rP[i]
			if opt:
				for i in range(len(rN)):
					xobj.getAttribute(rNlabels[i]).real = rN[i]
			for i in range(len(gK)):
				xobj.getAttribute(gKlabels[i]).real = gK[i]
			for i in range(len(gD)):
				xobj.getAttribute(gDlabels[i]).real = gD[i]
			for i in range(len(gF)):
				xobj.getAttribute(gFlabels[i]).real = gF[i]
			xobj.getAttribute('gE').real = gE
			xobj.getAttribute('dmgType').string = dType
		elif mat.name == 'ViscousDamper':
			p1 = mat.params[0]
			p2 = mat.params[1]
			p3 = mat.params[2]
			xobj.getAttribute('K').quantityScalar.value = p1[0]
			xobj.getAttribute('Cd').real = p1[1]
			xobj.getAttribute('alpha').real = p1[2]
			if p2:
				xobj.getAttribute('use_LGap').boolean = True
				xobj.getAttribute('LGap').real = p2
			if p3:
				xobj.getAttribute('Optional').boolean = True
				xobj.getAttribute('NM').integer = int(p3[0])
				xobj.getAttribute('RelTol').real = p3[1]
				xobj.getAttribute('AbsTol').real = p3[2]
				xobj.getAttribute('MaxHalf').integer = int(p3[3])
		elif mat.name == 'LimitState':
			xobj.getAttribute('s1p').quantityScalar.value = mat.params[0]
			xobj.getAttribute('e1p').real = mat.params[1]
			xobj.getAttribute('s2p').quantityScalar.value = mat.params[2]
			xobj.getAttribute('e2p').real = mat.params[3]
			xobj.getAttribute('s3p').quantityScalar.value = mat.params[4]
			xobj.getAttribute('e3p').real = mat.params[5]
			xobj.getAttribute('s1n').quantityScalar.value = mat.params[6]
			xobj.getAttribute('e1n').real = mat.params[7]
			xobj.getAttribute('s2n').quantityScalar.value = mat.params[8]
			xobj.getAttribute('e2n').real = mat.params[9]
			xobj.getAttribute('s3n').quantityScalar.value = mat.params[10]
			xobj.getAttribute('e3n').real = mat.params[11]
			xobj.getAttribute('pinchX').real = mat.params[12]
			xobj.getAttribute('pinchY').real = mat.params[13]
			xobj.getAttribute('damage1').real = mat.params[14]
			xobj.getAttribute('damage2').real = mat.params[15]
			xobj.getAttribute('beta').real = mat.params[16]
			xobj.getAttribute('curveTag').index = mat.params[17]
			xobj.getAttribute('curveType').integer = mat.params[18]
			xobj.getAttribute('Optional_1').boolean = True
			xobj.getAttribute('Optional_2').boolean = True
			xobj.getAttribute('use_beta').boolean = True
		else:
			continue
		prop = MpcProperty()
		prop.id = id
		prop.name = '{}_{}'.format(mat.name, id)
		prop.XObject = xobj
		callbacks.addPhysicalProperty(prop)
		prop.commitXObjectChanges()

# build cae: sections
def _build_cae_pprop_section(doc, callbacks):
	# fiber sections
	meta = callbacks.getMetaPhysicalProperty('sections.Fiber')
	for id, sec in doc.fibers.items():
		xobj = MpcXObject.createInstanceOf(meta)
		xobj.getAttribute('-GJ').boolean = True
		xobj.getAttribute('GJ').quantityScalar.value = sec.gj
		xobj.getAttribute('Y/section_offset').quantityScalar.value = sec.offy
		xobj.getAttribute('Z/section_offset').quantityScalar.value = sec.offz
		xobj.getAttribute('Dimension').string = '3D'
		fsec = MpcBeamFiberSection('Fiber section')
		for item in sec.patches:
			v1 = FxOccBuilder.makeVertex(item.p[0][0], item.p[0][1], 0.0)
			v2 = FxOccBuilder.makeVertex(item.p[1][0], item.p[1][1], 0.0)
			v3 = FxOccBuilder.makeVertex(item.p[2][0], item.p[2][1], 0.0)
			v4 = FxOccBuilder.makeVertex(item.p[3][0], item.p[3][1], 0.0)
			e1 = FxOccBuilder.makeEdge(v1, v2)
			e2 = FxOccBuilder.makeEdge(v2, v3)
			e3 = FxOccBuilder.makeEdge(v3, v4)
			e4 = FxOccBuilder.makeEdge(v4, v1)
			w1 = FxOccBuilder.makeWire([e1, e2, e3, e4])
			f1 = FxOccBuilder.makeFace(w1)
			surf = MpcBeamFiberSectionSurfaceFiberGroup('Patch', f1)
			surf.material = callbacks.getPhysicalProperty(item.mat)
			surf.meshControls.seed = MpcBeamFiberSurfaceMeshSeedType.ByNumber
			surf.meshControls.type = MpcBeamFiberSurfaceMeshType.Quadrilateral
			nx = item.nx
			ny = item.ny
			if nx != ny:
				p1 = Math.vec3(item.p[0][0], item.p[0][1], 0.0)
				p2 = Math.vec3(item.p[1][0], item.p[1][1], 0.0)
				dx = (p2-p1).normalized()
				if abs(dx.y) > abs(dx.x):
					nx, ny = ny, nx
			surf.meshControls.nx = nx
			surf.meshControls.ny = ny
			surf.makeMesh()
			fsec.addSurfaceFiber(surf)
		for item in sec.layers:
			v1 = FxOccBuilder.makeVertex(item.p[0][0], item.p[0][1], 0.0)
			v2 = FxOccBuilder.makeVertex(item.p[1][0], item.p[1][1], 0.0)
			e1 = FxOccBuilder.makeEdge(v1, v2)
			lin = MpcBeamFiberSectionPunctualFiberGroup('Layer', e1)
			lin.material = callbacks.getPhysicalProperty(item.mat)
			lin.edgeData.inputType = MpcBeamFiberPunctualEdgeDataInputType.ByNumber
			lin.edgeData.diameter = 2.0*math.sqrt(item.area/math.pi)
			lin.edgeData.number = item.num
			lin.generateRebarsLocations()
			lin.generateFibers()
			fsec.addPunctualFiber(lin)
		for item in sec.points:
			v1 = FxOccBuilder.makeVertex(item.y, item.z, 0.0)
			lin = MpcBeamFiberSectionPunctualFiberGroup('Point', v1)
			lin.material = callbacks.getPhysicalProperty(item.mat)
			lin.edgeData.inputType = MpcBeamFiberPunctualEdgeDataInputType.ByNumber
			lin.edgeData.diameter = 2.0*math.sqrt(item.a/math.pi)
			lin.edgeData.number = 1
			lin.generateRebarsLocations()
			lin.generateFibers()
			fsec.addPunctualFiber(lin)
		fsec.regenerateVisualRepresentation()
		fsec.commitChanges()
		xobj.getAttribute('Fiber section').customObject = fsec
		xobj.getAttribute('Y/section_offset').quantityScalar.value = sec.offy
		xobj.getAttribute('Z/section_offset').quantityScalar.value = sec.offz
		prop = MpcProperty()
		prop.id = id
		prop.name = 'FiberSection_{}'.format(id)
		prop.XObject = xobj
		callbacks.addPhysicalProperty(prop)
		prop.commitXObjectChanges()
	# elastic sections
	meta = callbacks.getMetaPhysicalProperty('sections.Elastic')
	for id, sec in doc.sections.items():
		xobj = MpcXObject.createInstanceOf(meta)
		xobj.getAttribute('E').quantityScalar.value = sec.E
		xobj.getAttribute('G/3D').quantityScalar.value = sec.G
		co = MpcBeamSection(
			MpcBeamSectionShapeType.Custom, 'section_Custom', 'user', 'm',
			[sec.A, sec.Jy, sec.Jz, sec.Jx, sec.Avy/sec.A, sec.Avz/sec.A])
		xobj.getAttribute('Section').customObject = co
		xobj.getAttribute('Y/section_offset').quantityScalar.value = sec.offy
		xobj.getAttribute('Z/section_offset').quantityScalar.value = sec.offz
		prop = MpcProperty()
		prop.id = id
		prop.name = 'ElasticSection_{}'.format(id)
		prop.XObject = xobj
		callbacks.addPhysicalProperty(prop)
		prop.commitXObjectChanges()
	# aggregators
	meta = callbacks.getMetaPhysicalProperty('sections.Aggregator')
	for id, sec in doc.aggregators.items():
		xobj = MpcXObject.createInstanceOf(meta)
		if sec.section is not None:
			xobj.getAttribute('-section').string = 'UseSectionTag'
			xobj.getAttribute('UseSectionTag').boolean = True
			xobj.getAttribute('sectionTag').index = sec.section
		for i in sec.materials:
			xobj.getAttribute(i[0]).boolean = True
			xobj.getAttribute('matTag{}'.format(i[0])).index = i[1]
		prop = MpcProperty()
		prop.id = id
		prop.name = 'Aggregator_{}'.format(id)
		prop.XObject = xobj
		callbacks.addPhysicalProperty(prop)
		prop.commitXObjectChanges()
	# plate sections
	meta = callbacks.getMetaPhysicalProperty('sections.ElasticMembranePlateSection')
	for id, sec in doc.plate_sections.items():
		xobj = MpcXObject.createInstanceOf(meta)
		xobj.getAttribute('E').quantityScalar.value = sec.E
		xobj.getAttribute('nu').real = sec.nu
		xobj.getAttribute('h').quantityScalar.value = sec.h
		xobj.getAttribute('rho').quantityScalar.value = sec.rho
		xobj.getAttribute('Ep_mod').real = sec.ep_mod
		prop = MpcProperty()
		prop.id = id
		prop.name = 'PlateSection_{}'.format(id)
		prop.XObject = xobj
		callbacks.addPhysicalProperty(prop)
		prop.commitXObjectChanges()
	# beam properties
	meta = callbacks.getMetaPhysicalProperty('special_purpose.BeamSectionProperty')
	for id, prop in doc.beam_props.items():
		xobj = MpcXObject.createInstanceOf(meta)
		xobj.getAttribute('secTag/1').index = prop.section
		xobj.getAttribute('numIntPts/1').integer = prop.nump
		prop = MpcProperty()
		prop.id = id
		prop.name = 'BeamProperty_{}'.format(id)
		prop.XObject = xobj
		callbacks.addPhysicalProperty(prop)
		prop.commitXObjectChanges()

# build cae: link properties
def _build_cae_pprop_link(doc, callbacks):
	meta_zl = callbacks.getMetaPhysicalProperty('special_purpose.zeroLengthMaterial')
	meta_2l = callbacks.getMetaPhysicalProperty('special_purpose.twoNodeLinkMaterial')
	for id, prop in doc.link_materials.items():
		if prop.type == tclimport.tclobjects.link_material_t.LINK:
			xobj = MpcXObject.createInstanceOf(meta_2l)
		else:
			xobj = MpcXObject.createInstanceOf(meta_zl)
		xobj.getAttribute('Dimension').string = '3D'
		xobj.getAttribute('2D').boolean = False
		xobj.getAttribute('3D').boolean = True
		xobj.getAttribute('ModelType').string = 'U-R (Displacement+Rotation)'
		xobj.getAttribute('U (Displacement)').boolean = False
		xobj.getAttribute('U-R (Displacement+Rotation)').boolean = True
		for i in range(len(prop.mats)):
			imat = prop.mats[i]
			idir = prop.dirs[i]
			if idir == 3:
				xobj.getAttribute('matTag{}/3D'.format(idir)).index = imat
				xobj.getAttribute('dir{}/3D'.format(idir)).boolean = True
			else:
				xobj.getAttribute('matTag{}'.format(idir)).index = imat
				xobj.getAttribute('dir{}'.format(idir)).boolean = True
		prop = MpcProperty()
		prop.id = id
		prop.name = 'Link_MaterialType_{}'.format(id)
		prop.XObject = xobj
		callbacks.addPhysicalProperty(prop)
		prop.commitXObjectChanges()

# build cae: physical properties
def _build_cae_pprop(doc, callbacks):
	_build_cae_pprop_uniaxial(doc, callbacks)
	_build_cae_pprop_section(doc, callbacks)
	_build_cae_pprop_link(doc, callbacks)

# build cae: assign properties
def _build_cae_assign_pprop(doc, callbacks):
	from time import time
	pmap = {}
	def pmap_get(id):
		items = pmap.get(id, None)
		if items is None:
			items = []
			pmap[id] = items
		return items
	for _, ele in doc.elements.items():
		id = ele.params[0]
		items = pmap_get(id)
		items.append(ele.cae_source)
	for _, ele in doc.links.items():
		id = ele.params[0]
		items = pmap_get(id)
		items.append(ele.cae_source)
	for id, items in pmap.items():
		for cae_source in items:
			callbacks.select(cae_source)
		callbacks.runCommand('AssignPhysicalProperty', str(id))
		callbacks.unselect()

# build cae: element properties
def _build_cae_eprop(doc, callbacks):
	new_id = 1
	# the property map
	pmap = {}
	def pmap_get(id):
		items = pmap.get(id, None)
		if items is None:
			items = []
			pmap[id] = items
		return items
	# beams
	aux = {}
	for ele_id, ele in doc.elements.items():
		if ele.name in tclimport.tclobjects.settings.BEAMS:
			name = ele.name
			ttype = ele.params[-1]
			rcy, rcz = ele.params[-3] # (sec, (rcy,rcz), locax, ttype)
			key = (name, ttype, rcy, rcz)
			if key in aux:
				prop = aux[key]
			else:
				meta = callbacks.getMetaElementProperty('beam_column_elements.{}'.format(name))
				xobj = MpcXObject.createInstanceOf(meta)
				if name in tclimport.tclobjects.settings.BEAMS_WITH_SEC:
					attribute = 'transType'
				else:
					attribute = 'transfType'
				xobj.getAttribute(attribute).string = ttype
				xobj.getAttribute('Dimension').string = '3D'
				xobj.getAttribute('2D').boolean = False
				xobj.getAttribute('3D').boolean = True
				# handle release codes... now only in elasticBeamColumn
				if name == 'elasticBeamColumn':
					if rcy > 0:
						xobj.getAttribute('-releasey').boolean = True
						xobj.getAttribute('releaseyCode').string = _utils_release.release_code_map[rcy]
					if rcz > 0:
						xobj.getAttribute('-releasez').boolean = True
						xobj.getAttribute('releasezCode').string = _utils_release.release_code_map[rcz]
				prop = MpcElementProperty()
				prop.id = new_id
				prop.name = 'BeamType_{}_{}_RC({},{})'.format(name, ttype, rcy, rcz)
				prop.XObject = xobj
				callbacks.addElementProperty(prop)
				prop.commitXObjectChanges()
				aux[key] = prop
				new_id += 1
			# assign
			items = pmap_get(prop.id)
			items.append(ele.cae_source)
	# walls
	aux = {}
	aux_suffix = 0
	for ele_id, ele in doc.elements.items():
		if ele.name in tclimport.tclobjects.settings.WALLS:
			key = str((ele.name, ele.params))
			if key in aux:
				prop, suffix = aux[key]
			else:
				meta = callbacks.getMetaElementProperty('beam_column_elements.{}'.format(ele.name))
				xobj = MpcXObject.createInstanceOf(meta)
				xobj.getAttribute('m').integer = ele.params[0]
				xobj.getAttribute('Thicknesses').quantityVector.value = ele.params[1]
				xobj.getAttribute('Widths').quantityVector.value = ele.params[2]
				xobj.getAttribute('Reinforcing_ratios').quantityVector.value = ele.params[3]
				xobj.getAttribute('Concrete_tags').indexVector = ele.params[4]
				xobj.getAttribute('Steel_tags').indexVector = ele.params[5]
				xobj.getAttribute('Shear_tag').index = ele.params[6]
				if ele.params[7]:
					xobj.getAttribute('-CoR').boolean = True
					xobj.getAttribute('c').quantityScalar.value = ele.params[7]
				prop = MpcElementProperty()
				prop.id = new_id
				aux_suffix += 1
				prop.name = 'WallType_{}({})'.format(ele.name, aux_suffix)
				prop.XObject = xobj
				callbacks.addElementProperty(prop)
				prop.commitXObjectChanges()
				aux[key] = (prop, aux_suffix)
				new_id += 1
			# assign
			items = pmap_get(prop.id)
			items.append(ele.cae_source)
	# shells
	aux = {}
	for ele_id, ele in doc.elements.items():
		if ele.name in tclimport.tclobjects.settings.SHELLS:
			key = ele.name
			if key in aux:
				prop = aux[key]
			else:
				meta = callbacks.getMetaElementProperty('shell.{}'.format(key))
				xobj = MpcXObject.createInstanceOf(meta)
				prop = MpcElementProperty()
				prop.id = new_id
				prop.name = 'ShellType_{}'.format(key)
				prop.XObject = xobj
				callbacks.addElementProperty(prop)
				prop.commitXObjectChanges()
				aux[key] = prop
				new_id += 1
			# assign
			items = pmap_get(prop.id)
			items.append(ele.cae_source)
	# links
	if len(doc.links) > 0:
		prop_zl = None
		prop_zs = None
		prop_2l = None
		for ele_id, ele in doc.links.items():
			# check
			if ele.name == 'zeroLength':
				if prop_zl is None:
					meta = callbacks.getMetaElementProperty('zero_length_elements.zeroLength')
					xobj = MpcXObject.createInstanceOf(meta)
					xobj.getAttribute('Dimension').string = '3D'
					xobj.getAttribute('2D').boolean = False
					xobj.getAttribute('3D').boolean = True
					prop = MpcElementProperty()
					prop.id = new_id
					prop.name = 'LinkType_zeroLength'
					prop.XObject = xobj
					callbacks.addElementProperty(prop)
					prop.commitXObjectChanges()
					new_id += 1
					prop_zl = prop
				else:
					prop = prop_zl
			elif ele.name == 'zeroLengthSection':
				if prop_zs is None:
					meta = callbacks.getMetaElementProperty('zero_length_elements.zeroLengthSection')
					xobj = MpcXObject.createInstanceOf(meta)
					xobj.getAttribute('Dimension').string = '3D'
					xobj.getAttribute('2D').boolean = False
					xobj.getAttribute('3D').boolean = True
					prop = MpcElementProperty()
					prop.id = new_id
					prop.name = 'LinkType_zeroLengthSection'
					prop.XObject = xobj
					callbacks.addElementProperty(prop)
					prop.commitXObjectChanges()
					new_id += 1
					prop_zs = prop
				else:
					prop = prop_zs
			elif ele.name == 'twoNodeLink':
				if prop_2l is None:
					meta = callbacks.getMetaElementProperty('link_elements.twoNodeLink')
					xobj = MpcXObject.createInstanceOf(meta)
					xobj.getAttribute('Dimension').string = '3D'
					xobj.getAttribute('2D').boolean = False
					xobj.getAttribute('3D').boolean = True
					xobj.getAttribute('ModelType').string = 'U-R (Displacement+Rotation)'
					xobj.getAttribute('U (Displacement)').boolean = False
					xobj.getAttribute('U-R (Displacement+Rotation)').boolean = True
					lparams = ele.params[2]
					if lparams[4]:
						xobj.getAttribute('-orient').boolean = True
					pDelta = lparams[0]
					if pDelta:
						xobj.getAttribute('-pDelta').boolean = True
						xobj.getAttribute('Mratios').quantityVector.value = pDelta
					shearDist = lparams[1]
					if shearDist:
						xobj.getAttribute('-shearDist').boolean = True
						xobj.getAttribute('sDratios').quantityVector.value = shearDist
					if lparams[2]:
						xobj.getAttribute('-doRayleigh').boolean = True
					if lparams[3]:
						xobj.getAttribute('-mass').boolean = True
						xobj.getAttribute('m').quantityScalar.value = lparams[3]
					prop = MpcElementProperty()
					prop.id = new_id
					prop.name = 'LinkType_twoNodeLink'
					prop.XObject = xobj
					callbacks.addElementProperty(prop)
					prop.commitXObjectChanges()
					new_id += 1
					prop_2l = prop
				else:
					prop = prop_2l
			else:
				continue
			# assign
			items = pmap_get(prop.id)
			items.append(ele.cae_source)
	for id, items in pmap.items():
		for cae_source in items:
			callbacks.select(cae_source)
		callbacks.runCommand('AssignElementProperty', str(id))
		callbacks.unselect()

# build cae: mass
def _build_cae_mass(doc, callbacks, support):
	# map:
	# key = tuple(mass values)
	# value = tuple(MpcCondition, nodes)
	M = {}
	R = {}
	for _, node in doc.nodes.items():
		if node.mass is not None:
			# split M and R
			iM = (node.mass[0], node.mass[1], node.mass[2])
			iR = (node.mass[3], node.mass[4], node.mass[5])
			# process M
			if math.sqrt(iM[0]**2 + iM[1]**2 + iM[2]**2) > 1.0e-14:
				if iM in M:
					value = M[iM]
				else:
					con = MpcCondition()
					con.id = support.condition_id
					support.condition_id += 1
					con.name = 'Translational Mass {}'.format(iM)
					value = (con, [])
					meta = callbacks.getMetaCondition('Mass.NodeMass')
					xobj = MpcXObject.createInstanceOf(meta)
					xobj.getAttribute('mass').quantityVector3.value = Math.vec3(iM[0], iM[1], iM[2])
					con.XObject = xobj
					M[iM] = value
				value[1].append(node)
			# process R
			if math.sqrt(iR[0]**2 + iR[1]**2 + iR[2]**2) > 1.0e-14:
				if iR in R:
					value = R[iR]
				else:
					con = MpcCondition()
					con.id = support.condition_id
					support.condition_id += 1
					con.name = 'Rotational Mass {}'.format(iR)
					value = (con, [])
					meta = callbacks.getMetaCondition('Mass.NodeRotationalMass')
					xobj = MpcXObject.createInstanceOf(meta)
					xobj.getAttribute('mass').quantityVector3.value = Math.vec3(iR[0], iR[1], iR[2])
					con.XObject = xobj
					R[iR] = value
				value[1].append(node)
	# assign
	for MM in (M, R):
		for _, value in MM.items():
			con = value[0]
			nodes = value[1]
			# obtain MpcConditionIndexedSubSet (map vertices to geometries)
			gmap = {}
			for node in nodes:
				source = node.cae_source
				if len(source) != 3:
					raise Exception('Mass supported only on geometries')
				geom = source[0]
				vertex = source[1]
				if geom in gmap:
					vertices = gmap[geom]
				else:
					vertices = []
					gmap[geom] = vertices
				vertices.append(vertex)
			for geom, vertices in gmap.items():
				sset = MpcConditionIndexedSubSet()
				for i in vertices:
					sset.vertices.append(i)
				con.assignTo(geom, sset)
			# done
			callbacks.addCondition(con)
			con.commitXObjectChanges()

# build cae: fix
def _build_cae_fix(doc, callbacks, support):
	# map:
	# key = tuple(fix values)
	# value = tuple(MpcCondition, nodes)
	F = {}
	for _, fix in support.stage.fix.items():
		if fix.dofs in F:
			value = F[fix.dofs]
		else:
			con = MpcCondition()
			con.id = support.condition_id
			support.condition_id += 1
			con.name = 'Fix ({})'.format(fix.dofstr())
			value = (con, [])
			meta = callbacks.getMetaCondition('Constraints.sp.fix')
			xobj = MpcXObject.createInstanceOf(meta)
			xobj.getAttribute('Dimension').string = '3D'
			xobj.getAttribute('ModelType').string = 'U-R (Displacement+Rotation)'
			xobj.getAttribute('2D').boolean = False
			xobj.getAttribute('3D').boolean = True
			xobj.getAttribute('U (Displacement)').boolean = False
			xobj.getAttribute('U-P (Displacement+Pressure)').boolean = False
			xobj.getAttribute('U-R (Displacement+Rotation)').boolean = True
			for label in fix.doflabels():
				xobj.getAttribute(label).boolean = True
			con.XObject = xobj
			F[fix.dofs] = value
		for node in fix.nodes:
			value[1].append(doc.nodes[node])
	# assign
	for _, value in F.items():
		con = value[0]
		nodes = value[1]
		# obtain MpcConditionIndexedSubSet (map vertices to geometries)
		gmap = {}
		for node in nodes:
			source = node.cae_source
			if len(source) != 3:
				raise Exception('Fix supported only on geometries')
			geom = source[0]
			vertex = source[1]
			if geom in gmap:
				vertices = gmap[geom]
			else:
				vertices = []
				gmap[geom] = vertices
			vertices.append(vertex)
		for geom, vertices in gmap.items():
			sset = MpcConditionIndexedSubSet()
			for i in vertices:
				sset.vertices.append(i)
			con.assignTo(geom, sset)
		# done
		callbacks.addCondition(con)
		con.commitXObjectChanges()
	# return conditions
	return [value[0] for value in F.values()]

# build cae: edof
def _build_cae_edof(doc, callbacks, support):
	res = []
	for _, edof in support.stage.edofs.items():
		con = MpcCondition()
		con.id = support.condition_id
		support.condition_id += 1
		con.name = 'EDOF ({})'.format(edof.dofstr())
		meta = callbacks.getMetaCondition('Constraints.mp.equalDOF')
		xobj = MpcXObject.createInstanceOf(meta)
		xobj.getAttribute('Dimension').string = '3D'
		xobj.getAttribute('ModelType').string = 'U-R (Displacement+Rotation)'
		xobj.getAttribute('2D').boolean = False
		xobj.getAttribute('3D').boolean = True
		xobj.getAttribute('U (Displacement)').boolean = False
		xobj.getAttribute('U-P (Displacement+Pressure)').boolean = False
		xobj.getAttribute('U-R (Displacement+Rotation)').boolean = True
		for label in edof.doflabels():
				xobj.getAttribute(label).boolean = True
		con.XObject = xobj
		con.assignTo(edof.cae_source)
		callbacks.addCondition(con)
		con.commitXObjectChanges()
		res.append(con)
	return res

# build cae: rigid diaphragm
def _build_cae_rdia(doc, callbacks, support):
	res = []
	for _, dia in support.stage.diaphragms.items():
		con = MpcCondition()
		con.id = support.condition_id
		support.condition_id += 1
		con.name = 'RigidDiaphragm ({})'.format(dia.master)
		meta = callbacks.getMetaCondition('Constraints.mp.rigidDiaphragm')
		xobj = MpcXObject.createInstanceOf(meta)
		xobj.getAttribute('perpDirn').integer = 3
		con.XObject = xobj
		con.assignTo(dia.cae_source)
		callbacks.addCondition(con)
		con.commitXObjectChanges()
		res.append(con)
	return res

# build cae: loads
def _build_cae_load(doc, callbacks, support):
	# process all patterns
	for _, pat in support.stage.patterns.items():
		#
		# create pattern analysis step
		step = MpcAnalysisStep()
		step.id = support.step_id
		support.step_id += 1
		step.name = 'Load Pattern ({})'.format(pat.id)
		meta = callbacks.getMetaAnalysisStep('Patterns.addPattern.loadPattern')
		xobj = MpcXObject.createInstanceOf(meta)
		xobj.getAttribute('tsTag').index = 1
		indices = xobj.getAttribute('load').indexVector
		#
		# build loads/couples in this pattern
		aux_data = ( 
			(pat.forces, 'Node Force {}', 'Loads.Force.NodeForce', 'F'), 
			(pat.couples, 'Node Couple {}', 'Loads.Force.NodeCouple', 'M') )
		for fdata, name, metaname, atname in aux_data:
			for f, nodes in fdata.items():
				con = MpcCondition()
				con.id = support.condition_id
				support.condition_id += 1
				con.name = name.format(f)
				meta = callbacks.getMetaCondition(metaname)
				xobj_load = MpcXObject.createInstanceOf(meta)
				xobj_load.getAttribute(atname).quantityVector3.value = Math.vec3(f[0], f[1], f[2])
				con.XObject = xobj_load
				# obtain MpcConditionIndexedSubSet (map vertices to geometries)
				gmap = {}
				for node in nodes:
					source = doc.nodes[node].cae_source
					if len(source) != 3:
						raise Exception('Load supported only on geometries')
					geom = source[0]
					vertex = source[1]
					if geom in gmap:
						vertices = gmap[geom]
					else:
						vertices = []
						gmap[geom] = vertices
					vertices.append(vertex)
				for geom, vertices in gmap.items():
					sset = MpcConditionIndexedSubSet()
					for i in vertices:
						sset.vertices.append(i)
					con.assignTo(geom, sset)
				# done
				callbacks.addCondition(con)
				con.commitXObjectChanges()
				# add to pattern
				indices.append(con.id)
		#
		# add pattern
		step.XObject = xobj
		# done
		callbacks.addAnalysisStep(step)
		step.commitXObjectChanges()

# build cae: stages
def _build_cae_stages(doc, callbacks, support):
	# we can pre-build a linear time series
	meta = callbacks.getMetaDefinition('timeSeries.Linear')
	xobj = MpcXObject.createInstanceOf(meta)
	defi = MpcDefinition()
	defi.name = 'Linear Series'
	defi.id = 1
	defi.XObject = xobj
	callbacks.addDefinition(defi)
	defi.commitXObjectChanges()
	#
	# progressive step id
	support.step_id = 1
	#
	# process regions
	for _, region in doc.regions.items():
		step = MpcAnalysisStep()
		step.id = support.step_id
		support.step_id += 1
		step.name = 'Region ({})'.format(region.id)
		meta = callbacks.getMetaAnalysisStep('Misc_commands.region')
		xobj = MpcXObject.createInstanceOf(meta)
		is_nodal = (region.nodes is not None)
		if is_nodal:
			xobj.getAttribute('Region Type').string = 'Nodes'
		else:
			xobj.getAttribute('Region Type').string = 'Elements'
		xobj.getAttribute('Elements').boolean = not is_nodal
		xobj.getAttribute('Nodes').boolean = is_nodal
		xobj.getAttribute('SelectionSets').indexVector.append(region.cae_source.id)
		step.XObject = xobj
		# done
		callbacks.addAnalysisStep(step)
		step.commitXObjectChanges()
	#
	# process recorders
	# note: the original region id now becomes the ID of the region step 
	#(sequential and 1 based because regions are the first ones)
	for recorder in doc.recorders:
		step = MpcAnalysisStep()
		step.id = support.step_id
		support.step_id += 1
		step.name = 'Recorder ({})'.format(region.id)
		meta = callbacks.getMetaAnalysisStep('Misc_commands.customCommand')
		xobj = MpcXObject.createInstanceOf(meta)
		region_id = list(doc.regions.keys()).index(recorder.region) + 1
		xobj.getAttribute('TCLscript').string = 'recorder {} -region {} {}'.format(
			' '.join(recorder.pre_data),
			region_id,
			' '.join(recorder.post_data))
		step.XObject = xobj
		# done
		callbacks.addAnalysisStep(step)
		step.commitXObjectChanges()
	#
	# process stages
	for stage in doc.stages:
		support.stage = stage
		# add recorders
		# .. todo
		#
		# add constraints
		sp = _build_cae_fix(doc, callbacks, support)
		mp = _build_cae_edof(doc, callbacks, support) + _build_cae_rdia(doc, callbacks, support)
		#
		# make condition pattern
		step = MpcAnalysisStep()
		step.id = support.step_id
		support.step_id += 1
		step.name = 'Add Constraints ({})'.format(stage.id)
		meta = callbacks.getMetaAnalysisStep('Patterns.addPattern.constraintPattern')
		xobj = MpcXObject.createInstanceOf(meta)
		sp_list = xobj.getAttribute('sp').indexVector
		mp_list = xobj.getAttribute('mp').indexVector
		for item in sp:
			sp_list.append(item.id)
		for item in mp:
			mp_list.append(item.id)
		step.XObject = xobj
		#
		# add load patterns
		_build_cae_load(doc, callbacks, support)
		#
		# done
		callbacks.addAnalysisStep(step)
		step.commitXObjectChanges()

# build CAE model
def build_cae(doc, callbacks):
	support = _support()
	callbacks.sendLog(0.0, 'Building CAE model...')
	vmap, emap = _build_cae_model(doc, callbacks)
	callbacks.runCommand("Regenerate", "2")
	callbacks.sendLog(30.0, 'Building Definitions...')
	_build_cae_definitions(doc, callbacks, vmap)
	callbacks.sendLog(40.0, 'Building Local Axes...')
	_build_cae_locax(doc, callbacks)
	callbacks.sendLog(60.0, 'Building Physical Properties...')
	_build_cae_pprop(doc, callbacks)
	_build_cae_assign_pprop(doc, callbacks)
	callbacks.sendLog(80.0, 'Building Element Properties...')
	_build_cae_eprop(doc, callbacks)
	callbacks.sendLog(90.0, 'Building Analysis Steps...')
	_build_cae_mass(doc, callbacks, support)
	_build_cae_stages(doc, callbacks, support)