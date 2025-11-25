from PyMpc import *
from scipy.spatial import cKDTree
import numpy as np

#from opspro.parameters.ExpressionGuiTools import example_expression_line_edit
#example_expression_line_edit()

App.setColorWorkTreeByUsageFlag(False)
App.clearTerminal()
doc = App.caeDocument()

# geometry ID
geom_id = 1

# make a selection set with all isolators
sset_all_iso = doc.getSelectionSet(1)
print(f"Selection set '{sset_all_iso.name}' has {len(sset_all_iso.geometries[geom_id].edges)} edges")

# isolator file
fname = 'isolators.txt'

isol_coords = []
isol_types = []
with open(fname, 'r') as f:
	for line in f.read().splitlines():
		tok = line.split('\t')
		x, y, typ = float(tok[3]), float(tok[4]), int(tok[15])
		isol_coords.append((x, y))
		isol_types.append(typ)

isol_coords = np.array(isol_coords)
isol_types = np.array(isol_types)

# build KDTree
tree = cKDTree(isol_coords)

# tolerance in model units (e.g. meters or whatever your model uses)
tol = 0.1   # adjust as needed

# --- create selection sets ---
sset_map = {}
for key in np.unique(isol_types):
	sset = MpcSelectionSet()
	sset.name = f'ISOL_{key}'
	sset.id = doc.selectionSets.getlastkey(0) + 1
	gmap = MpcSelectionSetItem()
	gmap.wholeGeometry = False
	sset.geometries[geom_id] = gmap
	sset_map[key] = sset
	doc.addSelectionSet(sset)

geom = doc.getGeometry(geom_id)
N = len(sset_all_iso.geometries[geom_id].edges)
print(f"Total edges: {N}")

nfailed = 0

for edge_id in sset_all_iso.geometries[1].edges:
	# get first vertex position
	vertices = geom.shape.getSubshapeChildren(edge_id, MpcSubshapeType.Edge, MpcSubshapeType.Vertex)
	p = geom.shape.vertexPosition(vertices[0])
	query_point = np.array([[p.x, p.y]])

	# find closest isol point
	dist, idx = tree.query(query_point)
	if dist[0] <= tol:
		isol_type = isol_types[idx[0]]
		sset_type = sset_map[isol_type]
		sset_type.geometries[geom_id].edges.append(edge_id)
	else:
		nfailed += 1

print(f"Processed {N} edges, {nfailed} not matched (tolerance = {tol})")
print(__file__)