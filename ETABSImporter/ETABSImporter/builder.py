from PyMpc import *
from ETABSImporter.graph import graph
from ETABSImporter.document import document, frame, area
from ETABSImporter.stko_interface import stko_interface
import itertools
from typing import List, Dict, Tuple

# a utility class to store the components of the ETABS model
# that will be part of a single geometry in STKO
class _geometry:
    def __init__(self):
        # the index in STKO
        self.id = 0
        # The vertices dictionary is used to store the coordinates of the points in 3D space
        self.vertices : Dict[int, Math.vec3] = {}
        # The frames dictionary is used to store the connectivity of the frame members
        self.frames : Dict[int, frame] = {}
        # The areas dictionary is used to store the connectivity of the slab members
        self.areas : Dict[int, area] = {}
        # The geometry type
        self.type : TopAbs_ShapeEnum = TopAbs_ShapeEnum.TopAbs_COMPOUND

# this class is used to build the geometry of the model
# it uses the ETABS document to build the STKO document
class builder:

    def __init__(
            self,
            etabs_doc : document,
            stko : stko_interface,
        ):
        
        # the ETABS document to import in STKO
        self.etabs_doc = etabs_doc

        # the STKO document interface
        self.stko = stko

        # the geometries
        self._geoms : List[_geometry] = []

        # process
        # todo: merge coincident nodes in the ETABS model!
        try:
            self.stko.start()
            self._build_geometry_groups()
            self._build_geometries()
        except Exception as e:
            raise
        finally:
            self.stko.stop()

    # build the graph of the model.
    # the graph is used to find the connected components of the model.
    def _build_geometry_groups(self):
        # the number of vertices
        num_vertices = len(self.etabs_doc.vertices)
        # create a mapping from the ETABS node id to a 0-based index for the graph
        vertex_map = {}
        for node_id in self.etabs_doc.vertices.keys():
            vertex_map[node_id] = len(vertex_map)
        # create the graph and compute connected vertices
        the_graph = graph(num_vertices)
        def _add_ele_to_graph(ele):
            for i,j in itertools.combinations(ele.nodes, 2):
                the_graph.add_edge(vertex_map[i], vertex_map[j])
        for _, ele in self.etabs_doc.frames.items():
            _add_ele_to_graph(ele)
        for _, ele in self.etabs_doc.areas.items():
            _add_ele_to_graph(ele)
        vertex_groups = the_graph.connected_components()
        # now compute an inverse map to return the connected components in the original ETABS node id
        inv_vertex_map = {j:i for i,j in vertex_map.items()}
        # now convert the connected components to the original ETABS node id
        for i in range(len(vertex_groups)):
            vertex_groups[i] = [inv_vertex_map[j] for j in vertex_groups[i]]
        # now group frames and areas with geometries.
        # note: if a group is made by just 1 vertex, make it in an external geometry
        geom_floating_nodes = _geometry()
        self._geoms : List[_geometry] = []
        for group in vertex_groups:
            # if the group has only one vertex, make it in an external geometry
            if len(group) == 1:
                geom_floating_nodes.vertices[group[0]] = self.etabs_doc.vertices[group[0]]
                continue
            # create a new geometry for the group
            geom = _geometry()
            for node_id in group:
                geom.vertices[node_id] = self.etabs_doc.vertices[node_id]
            # now add the frames and areas to the geometry
            for ele_id, ele in self.etabs_doc.frames.items():
                for node_id in ele.nodes:
                    if node_id in group:
                        geom.frames[ele_id] = ele
                        break
            for ele_id, ele in self.etabs_doc.areas.items():
                for node_id in ele.nodes:
                    if node_id in group:
                        geom.areas[ele_id] = ele
                        break
            # determine geometry type and add (skip empty)
            if len(geom.frames) > 0:
                if len(geom.areas) > 0:
                    geom.type = TopAbs_ShapeEnum.TopAbs_COMPOUND
                else:
                    geom.type = TopAbs_ShapeEnum.TopAbs_WIRE if len(geom.frames) > 1 else TopAbs_ShapeEnum.TopAbs_EDGE
            elif len(geom.areas) > 0:
                geom.type = TopAbs_ShapeEnum.TopAbs_SHELL if len(geom.areas) > 1 else TopAbs_ShapeEnum.TopAbs_FACE
            else:
                continue
            self._geoms.append(geom)
        # add the floating nodes to the external geometry
        if len(geom_floating_nodes.vertices) > 0:
            geom_floating_nodes.type = TopAbs_ShapeEnum.TopAbs_COMPOUND if len(self._geoms) > 0 else TopAbs_ShapeEnum.TopAbs_VERTEX
            self._geoms.append(geom_floating_nodes)
    
    # adds geometries to STKO
    def _build_geometries(self):
        geom_id = self.stko.new_geometry_id()
        for geom in self._geoms:
            # store the STKO id
            geom.id = geom_id
            geom_id += 1
            # build vertices:
            # map the ETABS node id to the STKO vertex object
            V : Dict[int, FxOccShape] = {}
            for i, v in geom.vertices.items():
                V[i] = FxOccBuilder.makeVertex(v)
            # build frames:
            # map the ETABS frame id to the STKO edge object
            # also map the tuple of vertices to the STKO a the edge object (assume no overlapping edges!)
            # (will be used to obtain boundary edges for faces)
            # keep track of edges that are not connected to faces!
            E : Dict[int, FxOccShape] = {}
            EN : Dict[Tuple[int, int], int] = {}
            EF : List[int] = []
            for i, f in geom.frames.items():
                # create the edge from existing vertices in V
                n1, n2 = f.nodes
                edge = FxOccBuilder.makeEdge(V[n1], V[n2])
                E[i] = edge
                EN[tuple(sorted((n1, n2)))] = i
            # build areas
            # map the ETABS area id to the STKO face object
            F : Dict[int, FxOccShape] = {}
            for i, a in geom.areas.items():
                # create the face from existing edges in E
                # get the outer wire (only one supported now)
                outer_edges = []
                for j in range(len(a.nodes)):
                    n1 = a.nodes[j]
                    n2 = a.nodes[(j+1)%len(a.nodes)]
                    eid = EN.get(tuple(sorted((n1, n2))), None)
                    if eid is None:
                        # if the edge is not found, it's only a boundary edge, make it here
                        edge = FxOccBuilder.makeEdge(V[n1], V[n2])
                    else:
                        # the edge is also a frame, mark it as used in face
                        edge = E[eid]
                        EF.append(eid)
                    outer_edges.append(edge)
                outer_wire = FxOccBuilder.makeWire(outer_edges)
                # create the face from the edges
                face = FxOccBuilder.makeFace(outer_wire)
                F[i] = face
            # obtain floating edges
            EF = set(EF)
            EFloat = {i:j for i,j in E.items() if i not in EF}
            # make the final shape
            shape = None
            if geom.type == TopAbs_ShapeEnum.TopAbs_VERTEX:
                shape = list(V.values())[0]
            elif geom.type == TopAbs_ShapeEnum.TopAbs_EDGE:
                shape = list(EFloat.values())[0]
            elif geom.type == TopAbs_ShapeEnum.TopAbs_WIRE:
                # order floating edges in a so that the i-th edge shares a node with a previous one
                # EFloat_ordered = []
                # EN_inv = {j:i for i,j in EN.items()} # key = edge_id (same in EFloat), value = tuple of vertices
                # EVlist_source = [EN_inv[edge_id] for edge_id in EFloat.keys()]
                # EVlist = [EVlist_source.pop(0)] # pop first element
                # while EVlist_source:
                #     for _loc, _evitem in enumerate(EVlist_source):
                #         if any(_prev in _evitem for _prev in EVlist[-1]):
                #             EVlist.append(EVlist_source.pop(_loc))
                #             break
                # shape = FxOccBuilder.makeWire(list(EFloat[EN[n1n2]] for n1n2 in EVlist))
                wire_input = FxOccShapeVector()
                for val in EFloat.values():
                    wire_input.append(val)
                shapes_in_wire = FxOccFactoryCurve.makeWire(wire_input)
                shape = shapes_in_wire[0]
            elif geom.type == TopAbs_ShapeEnum.TopAbs_FACE:
                shape = list(F.values())[0]
            elif geom.type == TopAbs_ShapeEnum.TopAbs_SHELL:
                shape = FxOccBuilder.makeShell(list(F.values()))
            else:
                if len(E) == 0 and len(F) == 0:
                    # a compound of vertices
                    shape = FxOccBuilder.makeCompound(list(V.values()))
                else:
                    # if we are here, it's a compound of edges and faces
                    # join the values of E and F, concatenating the genetors
                    shape = FxOccBuilder.makeCompound(list(itertools.chain(EFloat.values(), F.values())))
            # create the STKO geometry
            geom_name = f'Geometry_{geom_id}'
            stko_geom = MpcGeometry(geom_id, geom_name, shape)
            self.stko.add_geometry(stko_geom)
        # run a regenerate command
        self.stko.regenerate()
            