from PyMpc import *
from scipy.spatial import cKDTree
import numpy as np

#from opspro.parameters.ExpressionGuiTools import example_expression_line_edit
#example_expression_line_edit()

App.setColorWorkTreeByUsageFlag(False)
App.clearTerminal()
doc = App.caeDocument()

sset_all_iso = doc.getSelectionSet(1)
num_isolators = len(sset_all_iso.geometries[1].edges)
print(f"Total isolators in selection set '{sset_all_iso.name}': {num_isolators}")



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
	sset.geometries[1] = gmap
	sset_map[key] = sset
	doc.addSelectionSet(sset)

geom = doc.getGeometry(1)
N = len(sset_all_iso.geometries[1].edges)
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
		sset_type.geometries[1].edges.append(edge_id)
	else:
		nfailed += 1

print(f"Processed {N} edges, {nfailed} not matched (tolerance = {tol})")