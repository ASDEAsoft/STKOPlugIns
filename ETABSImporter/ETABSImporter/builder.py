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

# maps an ETABS entity (vertex, frame, area) to a STKO geometry and subshape id
class _geometry_map_item:
    def __init__(self, geom:MpcGeometry, subshape_id:int):
        self.geom = geom
        self.subshape_id = subshape_id

# maps and ETABS entity (diaphragm) to a STKO interaction
class _interaction_map_item:
    def __init__(self, interaction:MpcInteraction):
        self.interaction = interaction

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

        # the grouped geometries (key = id in STKO)
        self._geoms : Dict[int, _geometry] = {}

        # mapping of the ETABS entities to the STKO geometries
        self._vertex_map : Dict[int, _geometry_map_item] = {}
        self._frame_map : Dict[int, _geometry_map_item] = {} 
        self._area_map : Dict[int, _geometry_map_item] = {}
        self._diaphram_map : Dict[str, _interaction_map_item] = {}

        # process
        try:
            self.stko.start()
            self._build_geometry_groups()
            self._build_geometries()
            self._build_interactions()
        except Exception as e:
            raise
        finally:
            self.stko.stop()

    # builds a graph of connected vertices
    def _make_vertex_groups(self) -> List[List[int]]:
        # the number of vertices
        num_vertices = len(self.etabs_doc.vertices)
        # create the graph and compute connected vertices
        the_graph = graph(num_vertices)
        def _add_ele_to_graph(ele):
            for i,j in itertools.combinations(ele.nodes, 2):
                the_graph.add_edge(i, j)
        for _, ele in self.etabs_doc.frames.items():
            _add_ele_to_graph(ele)
        for _, ele in self.etabs_doc.areas.items():
            _add_ele_to_graph(ele)
        # done
        return the_graph.connected_components()

    # build the graph of the model.
    # the graph is used to find the connected components of the model.
    def _build_geometry_groups(self):
        next_geom_id = self.stko.new_geometry_id()
        # vertex groups
        vertex_groups = self._make_vertex_groups()
        # now group frames and areas with geometries.
        # note: if a group is made by just 1 vertex, make it in an external geometry
        geom_floating_nodes = _geometry()
        for group in vertex_groups:
            # if the group has only one vertex, make it in an external geometry
            if len(group) == 1:
                geom_floating_nodes.vertices[group[0]] = self.etabs_doc.vertices[group[0]]
                continue
            # create a new geometry for the group
            geom = _geometry()
            geom.id = next_geom_id
            next_geom_id += 1
            # add the vertices to the geometry
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
            self._geoms[geom.id] = geom
        # add the floating nodes to the external geometry
        if len(geom_floating_nodes.vertices) > 0:
            geom_floating_nodes.id = next_geom_id
            next_geom_id += 1
            geom_floating_nodes.type = TopAbs_ShapeEnum.TopAbs_COMPOUND if len(self._geoms) > 0 else TopAbs_ShapeEnum.TopAbs_VERTEX
            self._geoms[geom_floating_nodes.id] = geom_floating_nodes
    
    # adds geometries to STKO
    def _build_geometries(self):
        for geom_id, geom in self._geoms.items():
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
                shapes_in_wire = FxOccFactoryCurve.makeWire(list(EFloat.values()))
                if len(shapes_in_wire) != 1:
                    raise Exception('Wire creation failed. Found more than one shape.')
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
                    # joint floating edges in 1 or multiple wires
                    float_wires = FxOccFactoryCurve.makeWire(list(EFloat.values()))
                    float_faces = FxOccFactorySurfaces.makeShell(list(F.values()), self.etabs_doc.tolerance, False, True)
                    shape = FxOccBuilder.makeCompound(list(itertools.chain(float_wires, float_faces)))
            # create the STKO geometry
            geom_name = f'Geometry_{geom.id}'
            stko_geom = MpcGeometry(geom.id, geom_name, shape)
            self.stko.add_geometry(stko_geom)
            self._make_map(stko_geom)

    # maps an ETABS _geometry to an STKO MpcGeometry
    def _make_map(self, geom : MpcGeometry):
        import time
        time_start = time.time()
        # map vertices
        tol = self.etabs_doc.tolerance
        def _vertex_key(v:Math.vec3) -> Tuple[int,int,int]:
            return (int(v.x/tol), int(v.y/tol), int(v.z/tol))
        unique_vertices : Dict[Tuple[int,int,int], int] = {}
        for i in range(geom.shape.getNumberOfSubshapes(MpcSubshapeType.Vertex)):
            key = _vertex_key(geom.shape.vertexPosition(i))
            if key not in unique_vertices:
                unique_vertices[key] = i
        # map edges
        def _edge_key(n1:int, n2:int) -> Tuple[int,int]:
            return tuple(sorted((n1, n2)))
        unique_edges : Dict[Tuple[int,int], int] = {}
        for i in range(geom.shape.getNumberOfSubshapes(MpcSubshapeType.Edge)):
            edge_vertices = geom.shape.getSubshapeChildren(i, MpcSubshapeType.Edge, MpcSubshapeType.Vertex)
            p1 = geom.shape.vertexPosition(edge_vertices[0])
            p2 = geom.shape.vertexPosition(edge_vertices[1])
            key = _edge_key(unique_vertices[_vertex_key(p1)], unique_vertices[_vertex_key(p2)])
            if key not in unique_edges:
                unique_edges[key] = i
        # map faces
        def _face_key(ekeys : List[Tuple[int, int]]) -> List[Tuple[int,int]]:
            return tuple(sorted(ekeys))
        unique_faces : Dict[List[Tuple[int,int]], int] = {}
        for i in range(geom.shape.getNumberOfSubshapes(MpcSubshapeType.Face)):
            face_edges = geom.shape.getSubshapeChildren(i, MpcSubshapeType.Face, MpcSubshapeType.Edge)
            ekeys : List[Tuple[int,int]] = []
            for j in range(len(face_edges)):
                edge_vertices = geom.shape.getSubshapeChildren(face_edges[j], MpcSubshapeType.Edge, MpcSubshapeType.Vertex)
                p1 = geom.shape.vertexPosition(edge_vertices[0])
                p2 = geom.shape.vertexPosition(edge_vertices[1])
                ekeys.append(_edge_key(unique_vertices[_vertex_key(p1)], unique_vertices[_vertex_key(p2)]))
            unique_faces[_face_key(ekeys)] = i
        # now we can map the ETABS entities to the STKO geometry/sub-geometry ids
        # using the ETABS entities in this geometry group
        etabs_geom = self._geoms[geom.id]
        # link the vertices to the STKO geometry
        for i, v in etabs_geom.vertices.items():
            key = _vertex_key(v)
            subshape_id = unique_vertices.get(key, None)
            if subshape_id is None:
                raise Exception(f'Vertex {i} not found in geometry {geom.id}')
            self._vertex_map[i] = _geometry_map_item(geom, subshape_id)
        # link the frames to the STKO geometry
        for i, f in etabs_geom.frames.items():
            # get the edge id from the map
            edge_vertices = f.nodes
            p1 = etabs_geom.vertices[edge_vertices[0]]
            p2 = etabs_geom.vertices[edge_vertices[1]]
            key = _edge_key(unique_vertices[_vertex_key(p1)], unique_vertices[_vertex_key(p2)])
            subshape_id = unique_edges.get(key, None)
            if subshape_id is None:
                raise Exception(f'Frame {i} not found in geometry {geom.id}')
            self._frame_map[i] = _geometry_map_item(geom, subshape_id)
        # link the areas to the STKO geometry
        for i, a in etabs_geom.areas.items():
            # get the face id from the map
            face_vertices = a.nodes
            ekeys : List[Tuple[int,int]] = []
            for j in range(len(face_vertices)):
                n1 = a.nodes[j]
                n2 = a.nodes[(j+1)%len(a.nodes)]
                p1 = etabs_geom.vertices[n1]
                p2 = etabs_geom.vertices[n2]
                key = _edge_key(unique_vertices[_vertex_key(p1)], unique_vertices[_vertex_key(p2)])
                ekeys.append(key)
            subshape_id = unique_faces.get(_face_key(ekeys), None)
            if subshape_id is None:
                raise Exception(f'Area {i} not found in geometry {geom.id}')
            self._area_map[i] = _geometry_map_item(geom, subshape_id)
        
    # adds interactions to STKO
    def _build_interactions(self):
        next_id = self.stko.new_interaction_id()
        # process diaphragms
        for name, items in self.etabs_doc.diaphragms.items():
            # get data
            retained_id = items[-1]
            constrained_ids = items[:-1]
            # make a new interaction
            interaction = MpcInteraction(next_id, f'Diaphragm {name}')
            next_id += 1
            # retained item
            retained_data = self._vertex_map.get(retained_id, None)
            if retained_data is None:
                raise Exception(f'Diagram {name} retained vertex {retained_id} not found in geometry')
            interaction.items.masters.append(MpcInteractionItem(retained_data.geom, MpcSubshapeType.Vertex, retained_data.subshape_id))
            # constrained items
            for i in constrained_ids:
                constrained_data = self._vertex_map.get(i, None)
                if constrained_data is None:
                    raise Exception(f'Diagram {name} constrained vertex {i} not found in geometry')
                interaction.items.slaves.append(MpcInteractionItem(constrained_data.geom, MpcSubshapeType.Vertex, constrained_data.subshape_id))
            # add the interaction to the document
            self.stko.add_interaction(interaction)
            # map the diaphragm to the interaction
            self._diaphram_map[name] = _interaction_map_item(interaction)
