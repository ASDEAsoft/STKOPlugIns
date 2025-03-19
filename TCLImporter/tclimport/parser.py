import importlib
import tclimport.tclobjects
importlib.reload(tclimport.tclobjects)
from tclimport.tclobjects import *
import tclimport.mathutils
importlib.reload(tclimport.mathutils)
from tclimport.mathutils import *


class _parser_utils:
	prop_offset = 1000000

# get lines from file contents
def _parse_get_lines(doc, filecontents):
	# read all and convert multi-line commands into single-line ones
	# replace '{' with ' { ' and '}' with ' } ' just to make sure the split works
	lines = filecontents.replace(';','\n').replace('\\\n', ' ').replace('{', ' { ').replace('}', ' } ').split('\n')
	# split each line in words
	lines = [ [i.strip() for i in line.replace('\t', ' ').split(' ') if i] for line in lines if line ]
	# remove comments and empty lines
	lines = [i for i in lines if len(i) > 0 and not i[0].startswith('#')]
	# join commands in {} with a list of words
	aux = [w for w in lines]
	lines = []
	counter = 0
	num_lines = len(aux)
	while counter < num_lines:
		w = aux[counter]
		if w[-1] == '{':
			w = w[0:-1]
			ww = []
			counter += 1
			while counter < num_lines:
				wnext = aux[counter]
				counter += 1
				if wnext[-1] == '}':
					if len(wnext) > 1:
						ww.append(wnext[0:-1])
					break
				ww.append(wnext)
			w.append(ww)
			lines.append(w)
		else:
			lines.append(w)
			counter += 1
	return lines

# parse a node
def _parse_node(doc, words):
	n = len(words)
	id = int(words[1])
	pos = (float(words[2]), float(words[3]), float(words[4]))
	mass = None
	if len(words) > 5 and words[5] == '-mass':
		if n == 9:
			mass = (float(words[6]), float(words[7]), float(words[8]),
					0.0, 0.0, 0.0)
		elif n == 12:
			mass = (float(words[6]), float(words[7]), float(words[8]),
					float(words[9]), float(words[10]), float(words[11]))
	doc.nodes[id] = node_t(id, pos, mass)

# parse a mass command
def _parse_mass(doc, words):
	n = len(words)
	id = int(words[1])
	if n == 5:
		mass = (float(words[2]), float(words[3]), float(words[4]),
				0.0, 0.0, 0.0)
	elif n == 8:
		mass = (float(words[2]), float(words[3]), float(words[4]),
				float(words[5]), float(words[6]), float(words[7]))
	doc.nodes[id].mass = mass

# parse a limit curve
def _parse_limit_curve(doc, words):
	n = len(words)
	if (n == 15 or n == 17) and (words[1] == 'Rotation'):
		id = int(words[2])
		pos = 3
		params = []
		for i in range(5):
			params.append(int(words[pos]))
			pos += 1
		for i in range(7):
			params.append(float(words[pos]))
			pos += 1
		if n == 17:
			params.append(float(words[16]))
		else:
			params.append(None)
	try:
		doc.limit_curves[id] = limit_curve_t(id, params)
	except:
		print(words)
		raise

# parse an element
def _parse_element(doc, words):
	n = len(words)
	eltype = words[1]
	id = int(words[2])
	pos = 3
	done = True
	params = []
	if eltype in settings.BEAMS:
		type = element_t.EDGE
		nodes = [ doc.nodes[int(words[i+pos])] for i in range(2) ]
		pos += 2
		# specific parameters
		if eltype == 'elasticBeamColumn':
			params = [ float(words[i+pos]) for i in range(6) ]
			pos += 6
		elif eltype == 'ElasticTimoshenkoBeam':
			params = [ float(words[i+pos]) for i in range(8) ]
			pos += 8
		else:
			params = [ int(words[i+pos]) for i in range(2) ] # num int points / section tag
			pos += 2
		# store transformation tag here (append it at last)
		ttag = int(words[pos])
		# release pairs (default 0,0 for all)
		ry = 0
		rz = 0
		if eltype == 'elasticBeamColumn':
			try:
				where = words.index('-releasey')
				ry = int(words[where+1])
			except: pass
			try:
				where = words.index('-releasez')
				rz = int(words[where+1])
			except: pass
		params.append((ry, rz))
		# append transformation tag now
		params.append(ttag)
		
	elif eltype in settings.SHELLS:
		type = element_t.FACE
		nodes = [ doc.nodes[int(words[i+pos])] for i in range(4) ]
		pos += 4
		params = [int(words[pos])]
	
	elif eltype in settings.WALLS:
		type = element_t.FACE
		nodes = [ doc.nodes[int(words[i+pos])] for i in range(4) ]
		pos += 4
		NFib = int(words[pos])
		pos += 1
		loc = words.index('-thick')
		thick = [float(words[i+loc+1]) for i in range(NFib)]
		loc = words.index('-width')
		width = [float(words[i+loc+1]) for i in range(NFib)]
		loc = words.index('-rho')
		rho = [float(words[i+loc+1]) for i in range(NFib)]
		loc = words.index('-matConcrete')
		matConcrete = [int(words[i+loc+1]) for i in range(NFib)]
		loc = words.index('-matSteel')
		matSteel = [int(words[i+loc+1]) for i in range(NFib)]
		loc = words.index('-matShear')
		matShear = int(words[loc+1])
		try:
			CoR = float(words[words.index('-CoR')+1])
		except:
			CoR = None
		try:
			ThickMod = float(words[words.index('-ThickMod')+1])
		except:
			ThickMod = None
		try:
			Poisson = float(words[words.index('-Poisson')+1])
		except:
			Poisson = None
		try:
			Density = float(words[words.index('-Density')+1])
		except:
			Density = None
		params = [NFib, thick, width, rho, matConcrete, matSteel, matShear, CoR, ThickMod, Poisson, Density]
	
	elif eltype in settings.LINKS:
		type = element_t.LINK
		nodes = [ doc.nodes[int(words[i+pos])] for i in range(2) ]
		pos += 2
		if eltype == 'zeroLength':
			mats = []
			dirs = []
			try:
				dir_pos = words.index('-dir')
				num_mats = dir_pos - pos - 1
				for i in range(num_mats):
					mats.append(int(words[pos+i+1]))
					dirs.append(int(words[dir_pos+i+1]))
			except:
				mats = []
				dirs = []
			if '-orient' in words:
				ori_pos = words.index('-orient')
				ori = [float(words[i]) for i in range(ori_pos+1, ori_pos+7)]
			else:
				ori = [1.0, 0.0, 0.0,   0.0, 1.0, 0.0]
			params = [mats, dirs, ori]
		elif eltype == 'zeroLengthSection':
			sec = int(words[pos])
			pos += 1
			if '-orient' in words:
				ori_pos = words.index('-orient')
				ori = [float(words[i]) for i in range(ori_pos+1, ori_pos+7)]
			else:
				ori = [1.0, 0.0, 0.0,   0.0, 1.0, 0.0]
			do_rayleigh = ('-doRayleigh' in words)
			params = [sec, do_rayleigh, ori]
		elif eltype == 'twoNodeLink':
			mats = []
			dirs = []
			try:
				dir_pos = words.index('-dir')
				num_mats = dir_pos - pos - 1
				for i in range(num_mats):
					mats.append(int(words[pos+i+1]))
					dirs.append(int(words[dir_pos+i+1]))
				pos += 2*(num_mats+1)
			except:
				mats = []
				dirs = []
			nwords = len(words)
			nrem = nwords-pos
			orient = [1.0, 0.0, 0.0,   0.0, 1.0, 0.0]
			pDelta = None
			shearDist = None
			doRayleigh = False
			mass = None
			do_orient = False
			while pos < nwords:
				iw = words[pos]
				if iw == '-shearDist':
					shearDist = [float(words[i]) for i in range(pos+1, pos+3)]
					pos += 3
				elif iw == '-doRayleigh':
					doRayleigh = True
					pos += 1
				elif iw == '-pDelta':
					pDelta = [float(words[i]) for i in range(pos+1, pos+5)]
					pos += 5
				elif iw == '-mass':
					mass = float(words[pos+1])
					pos += 1
				elif iw == '-orient':
					do_orient = True
					try: # try all of them
						orient = [float(words[i]) for i in range(pos+1, pos+7)]
						pos += 7
					except:
						# try only y
						iy = [float(words[i]) for i in range(pos+1, pos+4)]
						ix = v3normalize(v3sub(nodes[1].pos,nodes[0].pos))
						orient = [ix[0], ix[1], ix[2], iy[0], iy[1], iy[2]]
						pos += 4
				else:
					pos += 1
			params = [mats, dirs, orient, [pDelta, shearDist, doRayleigh, mass, do_orient]]
	else:
		type = element_t.UNKNOWN
		done = False
	
	if done:
		if type == element_t.LINK:
			doc.links[id] = element_t(id, nodes, type, eltype, params)
		else:
			doc.elements[id] = element_t(id, nodes, type, eltype, params)

# parse geom transf
def _parse_geom_trans(doc, words):
	id = int(words[2])
	type = words[1]
	pos = 3
	vz = tuple([float(words[i+pos]) for i in range(3)])
	pos += 3
	n = len(words) - pos
	if n == 7:
		pos += 1
		offi = tuple([float(words[i+pos]) for i in range(3)])
		pos += 3
		offj = tuple([float(words[i+pos]) for i in range(3)])
		pos += 3
	else:
		offi = None
		offj = None
	doc.geom_trans[id] = geom_trans_t(id, type, vz, offi, offj)

# parse a fiber section
def _parse_fiber_sec(doc, words):
	id = int(words[2])
	gj = float(words[4])
	f = fibers_t(id)
	f.gj = gj
	for w in words[5]:
		if w[0] == 'patch':
			f.patches.append(
				fiber_patch_t(
					int(w[2]), int(w[3]), int(w[4]),
					(float(w[5]), float(w[6])),
					(float(w[7]), float(w[8])),
					(float(w[9]), float(w[10])),
					(float(w[11]), float(w[12]))
					)
				)
		elif w[0] == 'layer':
			f.layers.append(
				fiber_layer_t(
					int(w[2]), int(w[3]), float(w[4]),
					(float(w[5]), float(w[6])),
					(float(w[7]), float(w[8]))
					)
			)
		elif w[0] == 'fiber':
			f.points.append(fiber_point_t(int(w[4]), float(w[1]), float(w[2]), float(w[3])))
	doc.fibers[id] = f

# parse an aggregator
def _parse_aggregator(doc, words):
	n = len(words)
	id = int(words[2])
	a = aggregator_t(id)
	i = 3
	while(i < n):
		j = i+1
		wi = words[i]
		wj = words[j]
		if wi == '-section':
			a.section = int(wj)
		else:
			a.materials.append((wj, int(wi)))
		i += 2
	doc.aggregators[id] = a

# parse an elastic membrane plate section
def _parse_plate_section(doc, words):
	id = int(words[2])
	ep_mod = 1.0
	if len(words) > 7: ep_mod = float(words[7])
	s = plate_section_t(
		id, float(words[3]), float(words[4]), 
		float(words[5]), float(words[6]), ep_mod
		)
	doc.plate_sections[id] = s

# parse an elastic beam section
def _parse_section(doc, words):
	id = int(words[2])
	# $E $A $Iz $Iy $G $J <$alphaY $alphaZ>
	E = float(words[3])
	A = float(words[4])
	Iz = float(words[5])
	Iy = float(words[6])
	G = float(words[7])
	J = float(words[8])
	params = [E,G,A,J,Iy,Iz]
	s = section_t(id, params)
	if len(words) > 10:
		aY = float(words[9])
		aZ = float(words[10])
		s.Avy = self.A*aY
		s.Avz = self.A*aZ
	doc.sections[id] = s

# parse a uniaxial material
def _parse_material_1d(doc, words):
	name = words[1]
	id = int(words[2])
	n = len(words)
	if name == 'Elastic' or name == 'ENT':
		# $E
		params = [float(words[3])]
	elif name == 'Concrete01':
		# $fpc $epsc0 $fpcu $epsU
		params = [float(words[i]) for i in range(3, n)]
	elif name == 'Concrete02':
		# $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets
		params = [float(words[i]) for i in range(3, n)]
	elif name == 'Steel02':
		# $Fy $E $b $R0 $cR1 $cR2
		params = [float(words[i]) for i in range(3, n)]
	elif name == 'Hysteretic':
		# (18) $matTag $s1p $e1p $s2p $e2p <$s3p $e3p> $s1n $e1n $s2n $e2n <$s3n $e3n> $pinchX $pinchY $damage1 $damage2 <$beta>
		# (17) $matTag $s1p $e1p $s2p $e2p <$s3p $e3p> $s1n $e1n $s2n $e2n <$s3n $e3n> $pinchX $pinchY $damage1 $damage2
		# (13) $matTag $s1p $e1p $s2p $e2p $s1n $e1n $s2n $e2n $pinchX $pinchY $damage1 $damage2
		params = [float(words[i]) for i in range(3, n)]
	elif name == 'HystereticSM':
		# support only the backward compatible version
		nlast = 0
		for i in reversed(range(3, n)): # use reverse to get the last negative value
			value = float(words[i])
			if value < 0.0:
				nlast = i
				break
		nextra = n-nlast+1
		nvals = int((nlast-3+1)/2)
		ppos = [float(words[i]) for i in range(3, 3+nvals)]
		pneg = [float(words[i]) for i in range(3+nvals, 3+nvals*2)]
		pmis = [float(words[i]) for i in range(nlast+1, n)]
		params = [ppos, pneg, pmis]
	elif name == 'MinMax':
		# $otherTag <-min $minStrain> <-max $maxStrain>
		other = int(words[3])
		mmin = None
		mmax = None
		counter = 4
		while(counter < n):
			wi = words[counter]
			wj = words[counter+1]
			if wi == '-min':
				mmin = float(wj)
			elif wi == '-max':
				mmax = float(wj)
			counter += 2
		params = [other, mmin, mmax]
	elif name == 'Series':
		# $matTag $tag1 $tag2 ...
		# note: use float just in case the source file has a .0 ...
		params = [int(float(words[i])) for i in range(3, n)]
	elif name == 'Pinching4':
		# $matTag 
		# $ePf1 $ePd1 $ePf2 $ePd2 $ePf3 $ePd3 $ePf4 $ePd4 
		# <$eNf1 $eNd1 $eNf2 $eNd2 $eNf3 $eNd3 $eNf4 $eNd4> 
		# $rDispP $rForceP $uForceP 
		# <$rDispN $rForceN $uForceN > 
		# $gK1 $gK2 $gK3 $gK4 $gKLim $gD1 $gD2 $gD3 $gD4 $gDLim $gF1 $gF2 $gF3 $gF4 $gFLim $gE $dmgType
		nn = n - 3
		pos = 3
		eP = [float(words[i+pos]) for i in range(8)]
		pos += 8
		if nn == 39:
			eN = [float(words[i+pos]) for i in range(8)]
			pos += 8
		else:
			eN = [None]*8
		rP = [float(words[i+pos]) for i in range(3)]
		pos += 3
		if nn == 39:
			rN = [float(words[i+pos]) for i in range(3)]
			pos += 3
		else:
			rN = [None]*3
		gK = [float(words[i+pos]) for i in range(5)]
		pos += 5
		gD = [float(words[i+pos]) for i in range(5)]
		pos += 5
		gF = [float(words[i+pos]) for i in range(5)]
		pos += 5
		gE = float(words[pos])
		dType = words[pos+1]
		params = [eP, eN, rP, rN, gK, gD, gF, gE, dType]
	elif name == 'ViscousDamper':
		# $matTag
		# $K $Cd $alpha <$LGap> < $NM $RelTol $AbsTol $MaxHalf>
		nn = n - 3
		pos = 3
		p1 = [float(words[i+pos]) for i in range(3)]
		p2 = None
		p3 = None
		pos += 3
		if nn == 4:
			p2 = float(words[pos])
			pos += 1
		elif nn == 8:
			p2 = float(words[pos])
			pos += 1
			p3 = [float(words[i+pos]) for i in range(4)]
		params = [p1, p2, p3]
	elif name == 'LimitState':
		# $s1p $e1p $s2p $e2p $s3p $e3p $s1n $e1n $s2n $e2n $s3n $e3n $pinchX $pinchY $damage1 $damage2 $beta $curveTag $curveType
		nn = n - 3
		pos = 3
		params = [None]*19
		for i in range(17):
			params[i] = float(words[pos+i])
		params[17] = int(words[pos+17])
		params[18] = int(words[pos+18])
	else:
		params = []
	mat = material_1d_t(id, name, params)
	doc.materials_1d[id] = mat

# parse regions
def _parse_regions(doc, words):
	id = int(words[1])
	n = len(words)
	pos = 2
	nodes = []
	eles = []
	ijob = 0
	while pos < n:
		wi = words[pos]
		if wi == '-nodeRange':
			n1 = int(words[pos+1])
			n2 = int(words[pos+2])
			for i in range(n1, n2+1):
				nodes.append(i)
			pos += 2
		if wi == '-eleRange':
			n1 = int(words[pos+1])
			n2 = int(words[pos+2])
			for i in range(n1, n2+1):
				eles.append(i)
			pos += 2
		elif wi == '-node':
			ijob = 1
		elif wi == '-ele':
			ijob = 2
		else:
			if ijob != 0:
				try:
					theint = int(wi)
				except:
					ijob = 0
				if ijob == 1:
					nodes.append(int(wi))
				elif ijob == 2:
					eles.append(int(wi))
		pos += 1
	if len(nodes) == 0:
		nodes = None
	if len(eles) == 0:
		eles = None
	doc.regions[id] = region_t(id, the_nodes = nodes, the_elements = eles)

# parse recorder
def _parse_recorders(doc, words):
	w0 = words[0]
	if w0 == 'Node':
		key = '-node'
	elif w0 == 'Element':
		key = '-ele'
	else:
		print(w0)
		raise Exception("Invalid recorder format")
	a = words.index(key)
	b = a+1
	while True:
		try:
			int(words[b])
			b += 1
		except:
			break
	rid = 1
	if len(doc.regions) > 0:
		rid = max(doc.regions.keys())+1
	pre = words[:a]
	reg = ['region', rid] + words[a:b]
	post = words[b:]
	doc.recorders.append(recorder_t(pre, rid, post))
	_parse_regions(doc, reg)

# parse fix
def _parse_fix(doc, words):
	node = int(words[1])
	pos = 2
	n = len(words)
	dofs = 0
	for i in range(6):
		iflag = words[i+pos] == '1'
		if iflag:
			dofs |= (1 << i)
	stage = doc.stages[-1]
	if dofs in stage.fix:
		fix = stage.fix[dofs]
	else:
		fix = fix_t(dofs)
		stage.fix[dofs] = fix
	fix.nodes.append(node)

# parse diaphragm
def _parse_diaphragm(doc, words):
	dir = int(words[1])
	master = int(words[2])
	key = (dir, master)
	stage = doc.stages[-1]
	if key in stage.diaphragms:
		dia = stage.diaphragms[key]
	else:
		dia = diaphragm_t(dir, master)
		stage.diaphragms[key] = dia
	for i in range(3, len(words)):
		dia.slaves.append(int(words[i]))

# parse equal dof
def _parse_edof(doc, words):
	master = int(words[1])
	slave = int(words[2])
	pos = 3
	dofs = 0
	for i in range(pos, len(words)):
		j = int(words[i])-1
		dofs |= (1 << j)
	stage = doc.stages[-1]
	key = (master, dofs)
	if key in stage.edofs:
		edof = stage.edofs[key]
	else:
		edof = edof_t(master, dofs)
		stage.edofs[key] = edof
	edof.slaves.append(slave)

# renumber physical properties
def _renumber_props(doc, stage = 1):
	
	# here we give an offset for the different components
	# that will go into STKO's physical properties.
	# So that we have some room for inserting new ones
	# during parsing if necessary
	offset = _parser_utils.prop_offset
	
	# here we also renumber definitions (limitCurves) starting with an offset
	# because ID = 1 is taken by the Linear (default) timeSeries
	new_id = offset+1
	aux = dict(doc.limit_curves)
	doc.limit_curves = {}
	for old_id, icurve in aux.items():
		if icurve.id != new_id:
			icurve.id = new_id
			# change entities referencing this curve (they can be:
			# uniaxial materials of type (LimiState))
			for _, mat in doc.materials_1d.items():
				if mat.name == 'LimitState':
					if mat.params[-2] == old_id:
						mat.params[-2] = new_id
		# add it back
		doc.limit_curves[icurve.id] = icurve
		# increment
		new_id += 1
	
	# renumber physical properties
	new_id = 1
	#
	# 1 - renumber uniaxial materials
	aux = dict(doc.materials_1d)
	doc.materials_1d = {}
	material_1d_map = {} # old to new id map (used for mat referencing mat)
	for old_id, imat in aux.items():
		if imat.id != new_id:
			# map old to new
			material_1d_map[old_id] = new_id
			# change id
			imat.id = new_id
			# change entities referencing this mat 1d (they can be:
			# fiber sections, aggregators, zero length elements (links), or WALLS)
			for _, fiber in doc.fibers.items():
				for item in fiber.patches:
					if item.mat == old_id:
						item.mat = new_id
				for item in fiber.layers:
					if item.mat == old_id:
						item.mat = new_id
			for _, aggr in doc.aggregators.items():
				for i in range(len(aggr.materials)):
					item = aggr.materials[i]
					if item[1] == old_id:
						aggr.materials[i] = (item[0], new_id)
			for _, ele in doc.links.items():
				if ele.name == 'zeroLength':
					mats = ele.params[0]
					for i in range(len(mats)):
						if mats[i] == old_id:
							mats[i] = new_id
			for _, ele in doc.elements.items():
				if ele.name == 'MVLEM_3D':
					for mats in (ele.params[4], ele.params[5]): # Concrete and Steel
						for i in range(len(mats)):
							if mats[i] == old_id:
								mats[i] = new_id
					if ele.params[6] == old_id: # Shear
						ele.params[6] = new_id
		# add it back
		doc.materials_1d[imat.id] = imat
		# increment
		new_id += 1
	#
	# 1.1 - update materials referencing other materials!
	for _, imat in doc.materials_1d.items():
		if imat.name == 'MinMax':
			other_id = imat.params[0]
			other_id_mod = material_1d_map.get(other_id, None)
			if other_id_mod is not None:
				imat.params[0] = other_id_mod
		elif imat.name == 'Series':
			for other_loc, other_id in enumerate(imat.params):
				other_id_mod = material_1d_map.get(other_id, None)
				if other_id_mod is not None:
					imat.params[other_loc] = other_id_mod
	#
	# 2 - renumber fiber sections
	if stage==1:
		new_id = offset*1 + 1
	aux = dict(doc.fibers)
	doc.fibers = {}
	for old_id, ifib in aux.items():
		if ifib.id != new_id:
			ifib.id = new_id
			# change entities referencing this fiber (they can be:
			# aggregators, elements, zeroLengthSection
			for _, aggr in doc.aggregators.items():
				for i in range(len(aggr.materials)):
					if aggr.section == old_id:
						aggr.section = new_id
			for _, ele in doc.elements.items():
				if ele.name in tclimport.tclobjects.settings.BEAMS_WITH_SEC:
					if ele.params[1] == old_id:
						ele.params[1] = new_id
			for _, ele in doc.links.items():
				if ele.name in tclimport.tclobjects.settings.LINKS_WITH_SEC:
					if ele.params[0] == old_id:
						ele.params[0] = new_id
		# add it back
		doc.fibers[ifib.id] = ifib
		# increment
		new_id += 1
	#
	# 3 - renumber elastic beam sections
	if stage==1:
		new_id = offset*2 + 1
	#
	# 4 - renumber aggregators
	if stage==1:
		new_id = offset*3 + 1
	aux = dict(doc.aggregators)
	doc.aggregators = {}
	for old_id, iaggr in aux.items():
		if iaggr.id != new_id:
			iaggr.id = new_id
			# change entities referencing this aggregator (they can be:
			# elements, zeroLengthSection)
			for _, ele in doc.elements.items():
				if ele.name in tclimport.tclobjects.settings.BEAMS_WITH_SEC:
					# 0 = int points , 1 = sec id -> from parse_element
					if ele.params[1] == old_id:
						ele.params[1] = new_id
			for _, ele in doc.links.items():
				if ele.name in tclimport.tclobjects.settings.LINKS_WITH_SEC:
					if ele.params[0] == old_id:
						ele.params[0] = new_id
		# add it back
		doc.aggregators[iaggr.id] = iaggr
		# increment
		new_id += 1
	#
	# 5 - renumber beam properties
	if stage==1:
		new_id = offset*4 + 1
	#
	# 6 - renumber plate sections
	if stage==1:
		new_id = offset*5 + 1
	aux = dict(doc.plate_sections)
	doc.plate_sections = {}
	for old_id, isec in aux.items():
		if isec.id != new_id:
			isec.id = new_id
			# change entities referencing this plate section (they can be:
			# elements)
			for _, ele in doc.elements.items():
				if ele.name in tclimport.tclobjects.settings.SHELLS:
					if ele.params[0] == old_id:
						ele.params[0] = new_id
		# add it back
		doc.plate_sections[isec.id] = isec
		# increment
		new_id += 1

# converts elastic beam properties into elastic beam sections
def _convert_elastic_beams(doc):
	new_id = _parser_utils.prop_offset*2+1
	aux = {}
	for _, ele in doc.elements.items():
		if ele.name in tclimport.tclobjects.settings.BEAMS_ELA:
			# here params contains all the elastic section parameters (first) and the transformation tag (last)
			if ele.name == 'elasticBeamColumn':
				s = section_t(0, [ele.params[i] for i in [1,2,0,3,4,5]])
			else:
				s = section_t(0, ele.params[0:8])
			if not s in aux:
				s.id = new_id
				aux[s] = new_id
				new_id += 1
				# get the id of this new section
				sec_id = s.id
			else:
				# get the id of the existing section
				sec_id = aux[s]
			# now convert the elastic beam parameters with [sec_id, releaseCodesYZ, transformation id]
			ele.params = [sec_id, ele.params[-2], ele.params[-1]]
	for k,_ in aux.items():
		doc.sections[k.id] = k

# disp and force beam columns in STKO use the new format
def _make_beam_props(doc):
	new_id = _parser_utils.prop_offset*4+1
	aux = {}
	for _, ele in doc.elements.items():
		if ele.name in tclimport.tclobjects.settings.BEAMS_WITH_SEC:
			nump = ele.params[0]
			sec = ele.params[1]
			key = (sec, nump)
			if key in aux:
				prop_id = aux[key]
			if not key in aux:
				prop_id = new_id
				aux[key] = new_id
				new_id += 1
			# now convert the beam parameters with [sec_id, releaseCodesYZ, transformation id]
			ele.params = [prop_id, ele.params[-2], ele.params[-1]]
	for key, prop_id in aux.items():
		p = beam_prop_t(prop_id)
		p.section = key[0]
		p.nump = key[1]
		doc.beam_props[prop_id] = p

# process offsets in geom trans to pass it to sections.
def _process_offset(doc):
	# maps original items with duplicated ones with offset
	# the key is a tuple (original_id, offy, offz) -> value = new id
	# section offset
	class key_t:
		def __init__(self, id, y, z):
			self.id = id
			self.y = y
			self.z = z
		def __hash__(self):
			return hash((self.id, self.y, self.z))
		def __eq__(self, b):
			return (self.id == b.id) and (abs(self.y-b.y)<1.0e-10) and (abs(self.z-b.z)<1.0e-10)
	# search for a key that containts the given id
	def has_id_in_key(map, id):
		for key in map:
			if key.id == id:
				return True
		return False
	# maps a section (Elastic or Fiber) to the proper id/offset key
	# returns the id of the mapped section with requested offset.
	# automatically creates a new one if necessary
	def remap_section(sec_id, oy, oz, collection, map, update_off):
		key = key_t(sec_id, oy, oz)
		if key in map:
			new_sec_id = map[key]
		else:
			# get old section
			sec = collection[sec_id]
			# if there is a key with this old sec_id (but with different offsets) ...
			if has_id_in_key(map, sec_id):
				# we copy and create a new one ...
				new_sec_id = max(collection) + 1
				sec = sec.copy(new_sec_id)
				collection[new_sec_id] = sec
			else:
				# ... otherwise we modify the existing one and use it as new_sec_id
				new_sec_id = sec_id
			# update now the offset and map it
			if update_off:
				sec.offy = oy
				sec.offz = oz
			map[key] = new_sec_id
		return new_sec_id
	# one map for every document component that may have offsets
	map_sec = {}
	map_fib = {}
	map_agg = {}
	map_pro = {}
	# a small tolerance
	tol = 1.0e-10
	# process all geometric transformations
	for transf_id, transf in doc.geom_trans.items():
		oi = transf.offi
		oj = transf.offj
		if oi is None or oj is None:
			oi = (0.0, 0.0, 0.0)
			oj = (0.0, 0.0, 0.0)
		noi = v3norm(oi)
		noj = v3norm(oj)
		diff = v3sub(oi, oj)
		if v3norm(diff) > tol:
			raise Exception("Beam-end offsets must be equal")
		# find all elements with this transformation
		for _, ele in doc.elements.items():
			if ele.name in tclimport.tclobjects.settings.BEAMS:
				if ele.params[-1] == transf_id: # transformation id at the end (sec, releaseCodesYZ, tid)
					p1 = ele.nodes[0].pos
					p2 = ele.nodes[1].pos
					dx = v3normalize(v3sub(p2,p1))
					dy = v3normalize(v3cross(transf.vz, dx))
					dz = v3normalize(v3cross(dx, dy))
					ox = v3dot(dx, oi)
					if abs(ox) > tol:
						raise Exception('Axial offset not supported')
					oy = v3dot(dy, oi)
					oz = v3dot(dz, oi)
					# find its section. sections which may have offsets are:
					# - Elastic
					# - Fiber
					# - At this point, an element can have either an elastic section,
					#   or a beam property (with fiber, elastic, or aggregator)
					sec_id = ele.params[0]
					if sec_id in doc.sections:
						new_sec_id = remap_section(sec_id, oy, oz, doc.sections, map_sec, True)
						ele.params[0] = new_sec_id
					elif sec_id in doc.beam_props:
						beam_prop = doc.beam_props[sec_id]
						if beam_prop.section in doc.aggregators:
							aggregator = doc.aggregators[beam_prop.section]
							if aggregator.section is not None:
								if aggregator.section in doc.fibers:
									# remap the beam_prop (do not update offset)
									new_sec_id = remap_section(sec_id, oy, oz, doc.beam_props, map_pro, False)
									ele.params[0] = new_sec_id
									# now we need to remap the inner aggregator
									new_sec = doc.beam_props[new_sec_id]
									new_sec.section = remap_section(new_sec.section, oy, oz, doc.aggregators, map_agg, True)
									# now we need to remap the inner fiber section
									new_aggr = doc.aggregators[new_sec.section]
									new_aggr.section = remap_section(new_aggr.section, oy, oz, doc.fibers, map_fib, True)
								elif aggregator.section in doc.sections:
									# remap the beam_prop (do not update offset)
									new_sec_id = remap_section(sec_id, oy, oz, doc.beam_props, map_pro, False)
									ele.params[0] = new_sec_id
									# now we need to remap the inner aggregator
									new_sec = doc.beam_props[new_sec_id]
									new_sec.section = remap_section(new_sec.section, oy, oz, doc.aggregators, map_agg, True)
									# now we need to remap the inner fiber section
									new_aggr = doc.aggregators[new_sec.section]
									new_aggr.section = remap_section(new_aggr.section, oy, oz, doc.fibers, map_fib, True)
						elif beam_prop.section in doc.fibers:
							# remap the beam_prop (do not update offset)
							new_sec_id = remap_section(sec_id, oy, oz, doc.beam_props, map_pro, False)
							ele.params[0] = new_sec_id
							# now we need to remap the inner fiber section
							new_sec = doc.beam_props[new_sec_id]
							new_sec.section = remap_section(new_sec.section, oy, oz, doc.fibers, map_fib, True)
						elif beam_prop.section in doc.sections:
							# remap the beam_prop (do not update offset)
							new_sec_id = remap_section(sec_id, oy, oz, doc.beam_props, map_pro, False)
							ele.params[0] = new_sec_id
							# now we need to remap the inner elstic section
							new_sec = doc.beam_props[new_sec_id]
							new_sec.section = remap_section(new_sec.section, oy, oz, doc.sections, map_sec, True)

# process local axes
def _process_axes(doc):
	# a map of new local axes to be put in STKO
	# if different from default ones
	class vz_t:
		def __init__(self, v):
			self.v = v
		def __hash__(self):
			return hash(self.v)
		def __eq__(self, b):
			return v3dot(self.v, b.v) > 0.999
	locax_id = 1
	locax_map = {}
	# process all elements
	for ele_id, ele in doc.elements.items():
		if ele.name in tclimport.tclobjects.settings.BEAMS:
			transf_id = ele.params[-1] # transf id at the end
			transf = doc.geom_trans[transf_id]
			# element local X axis
			p1 = ele.nodes[0].pos
			p2 = ele.nodes[1].pos
			dx = v3normalize(v3sub(p2, p1))
			# user-defined local Z axis (not necessary orthogonal to dx)
			Vz = transf.vz
			Vy = v3normalize(v3cross(Vz, dx))
			Vz = v3normalize(v3cross(dx, Vy))
			# now that Vz (user-defined) is ortho-normal to dx...
			# compute default orientation as per STKO
			dy = ( 0.0, -1.0, 0.0)
			dz = ( 0.0,  0.0, 1.0)
			# project on curve (tangent = dx)
			if abs(v3dot(dz, dx)) > 0.999:
				# use trial dy
				dz = v3normalize(v3cross(dx, dy))
				dy = v3normalize(v3cross(dz, dx))
			else:
				# use trial dz
				dy = v3normalize(v3cross(dz, dx))
				dz = v3normalize(v3cross(dx, dy))
			# compare
			if v3dot(dz, Vz) < 0.99:
				# we need to assign a different local axes
				# let's see if one with the same orientation exists
				key = vz_t(Vz)
				if key in locax_map:
					new_id = locax_map[key].id
				else:
					new_id = locax_id
					new_locax = geom_trans_t(locax_id, None, Vz, None, None)
					locax_id += 1
					locax_map[key] = new_locax
				# update element
				ele.params[-1] = new_id
			else:
				ele.params[-1] = None
			ele.params.append(transf.type)  # transf id at the end, now transf type at the end
	# remove old geom transf
	doc.geom_trans = {}
	for _, value in locax_map.items():
		doc.geom_trans[value.id] = value
	#
	# now we do the same for link elements, whose key should be Vx-Vy
	class vxvy_t:
		def __init__(self, vx, vy):
			self.vx = vx
			self.vy = vy
		def __hash__(self):
			return hash(self.vx+self.vy)
		def __eq__(self, b):
			return (v3dot(self.vx, b.vx) > 0.999) and (v3dot(self.vy, b.vy) > 0.999)
	# clear map
	locax_map = {}
	# create a default key with global orientation
	key_global = vxvy_t((1.0,0.0,0.0), (0.0,1.0,0.0))
	# process all links
	for ele_id, ele in doc.links.items():
		ori = ele.params[2]
		vx = (ori[0], ori[1], ori[2])
		vy = (ori[3], ori[4], ori[5])
		key = vxvy_t(vx, vy)
		if key == key_global:
			ele.params[2] = None # no local axes defined (use global ones)
			continue
		# map
		if key in locax_map:
			new_id = locax_map[key].id
		else:
			new_id = locax_id
			new_locax = geom_trans_t(locax_id, None, None, None, None)
			new_locax.setxy(vx, vy)
			locax_id += 1
			locax_map[key] = new_locax
		# update element
		ele.params[2] = new_id
	# add them to doc
	for _, value in locax_map.items():
		doc.geom_trans[value.id] = value

# process link materials
def _process_link_props(doc):
	new_id = _parser_utils.prop_offset*6 + 1
	# map key = tuple(mats, dirs, trans_id) to a link_mat
	aux = {}
	for _, ele in doc.links.items():
		if ele.name == 'zeroLength' or ele.name == 'twoNodeLink':
			mats = ele.params[0]
			dirs = ele.params[1]
			trans_id = ele.params[2]
			key = (tuple(mats), tuple(dirs), trans_id)
			if key in aux:
				pid = aux[key]
			else:
				p = tclimport.tclobjects.link_material_t(new_id, mats[:], dirs[:])
				aux[new_id] = p
				pid = new_id
				new_id += 1
			if ele.name == 'twoNodeLink':
				ele.params = [pid, trans_id, ele.params[3]]
				p.type = tclimport.tclobjects.link_material_t.LINK
			else:
				ele.params = [pid, trans_id]
		elif ele.name == 'zeroLengthSection':
			# make it compatible
			ele.params = [ele.params[0], ele.params[2]]
	# add to document
	for id, p in aux.items():
		doc.link_materials[id] = p

# parse patterns
def _ele_load_error(words):
	raise Exception('Non-supported eleLoad format: {}\nThe supported format is "eleLoad -ele $tag -type -beamUniform $Wy $Wz <$Wx>"'.format(words))
def _parse_pattern(doc, words):
	id = words[2]
	values = words[4]
	p = tclimport.tclobjects.pattern_t(id)
	for item in values:
		if item[0] == 'load':
			node = int(item[1])
			load = [float(item[i+2]) for i in range(6)]
			p.add_load(load, node)
		elif item[0] == 'eleLoad':
			if item[1] != '-ele': _ele_load_error(words)
			if item[3] != '-type': _ele_load_error(words)
			if item[4] != '-beamUniform': _ele_load_error(words)
			ele = int(item[2])
			if len(words) == 8:
				load = (float(item[5]), float(item[6]), float(item[7]))
			else:
				load = (float(item[5]), float(item[6]), 0.0)
			p.add_ele_load(load, ele)
			# eleLoad -ele 38 -type -beamUniform 0.0 0.0 -2.25
	p.finalize()
	doc.stages[-1].patterns[id] = p


# main parsing function
def parse_tcl(doc, filecontents):
	# before init, include the first stage in the doc
	doc.stages.append(stage_t(1))
	# 0. get lines
	lines = _parse_get_lines(doc, filecontents)
	# 1. parse nodes
	for words in lines:
		w0 = words[0]
		if w0 == 'node':
			_parse_node(doc, words)
	# 2. parse elements
	for words in lines:
		w0 = words[0]
		if w0 == 'element':
			_parse_element(doc, words)
	#
	# 3. parse all other stuff
	# store recorders lines for later... we need them after regions
	rec_lines = []
	for words in lines:
		w0 = words[0]
		if w0 == 'geomTransf':
			_parse_geom_trans(doc, words)
		elif w0 == 'limitCurve':
			_parse_limit_curve(doc, words)
		elif w0 == 'uniaxialMaterial':
			_parse_material_1d(doc, words)
		elif w0 == 'section':
			w1 = words[1]
			if w1 == 'fiberSec' or w1 == 'Fiber':
				_parse_fiber_sec(doc, words)
			elif w1 == 'Aggregator':
				_parse_aggregator(doc, words)
			elif w1 == 'ElasticMembranePlateSection':
				_parse_plate_section(doc, words)
			elif w1 == 'Elastic':
				_parse_section(doc, words)
		elif w0 == 'region':
			_parse_regions(doc, words)
		elif w0 == 'fix':
			_parse_fix(doc, words)
		elif w0 == 'rigidDiaphragm':
			_parse_diaphragm(doc, words)
		elif w0 == 'equalDOF':
			_parse_edof(doc, words)
		elif w0 == 'pattern' and words[1] == 'Plain':
			_parse_pattern(doc, words)
		elif w0 == 'recorder':
			rec_lines.append(words[1:])
		elif w0 == 'mass':
			_parse_mass(doc, words)
	#
	# 3.1 process all recorder lines
	for words in rec_lines:
		_parse_recorders(doc, words)
	#
	# 4. now we need to do some post processing
	#
	# - renumber mat1d/sections that go into physical properties,
	#   so their IDs cannot overlap with eachother
	_renumber_props(doc, stage=1)
	#
	# - elastic beam elements in stko need an elastic section
	_convert_elastic_beams(doc)
	#
	# - disp and force beam columns in new format
	_make_beam_props(doc)
	#
	# - process offsets in geom trans to pass it to sections.
	#   section will eventually be duplicated in case of offsets
	_process_offset(doc)
	#
	# - process local axes
	_process_axes(doc)
	#
	# - process link materials
	_process_link_props(doc)