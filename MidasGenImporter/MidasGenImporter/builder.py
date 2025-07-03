from PyMpc import *
from MidasGenImporter.graph import graph
from MidasGenImporter.document import *
from MidasGenImporter.stko_interface import stko_interface
from MidasGenImporter.units_utils import unit_system
from MidasGenImporter.localaxes_utils import MIDASLocalAxesConvention
import itertools
from typing import List, Dict, Tuple, DefaultDict, Union
import math
from collections import defaultdict

class _globals:
    fix_labels = ('Ux', 'Uy', 'Uz', 'Rx', 'Ry', 'Rz/3D')

# a utility class to store the components of the MIDAS model
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

# maps an MIDAS entity (vertex, frame, area) to a STKO geometry and subshape id
class _geometry_map_item:
    def __init__(self, geom:MpcGeometry, subshape_id:int):
        self.geom = geom
        self.subshape_id = subshape_id

# maps and ETABS entity (diaphragm) to a STKO interaction
class _interaction_map_item:
    def __init__(self, interaction:MpcInteraction):
        self.interaction = interaction

# this class is used to build the geometry of the model
# it uses the MIDAS document to build the STKO document
class builder:

    def __init__(
            self,
            midas_doc : document,
            stko : stko_interface,
        ):
        
        # the MIDAS document to import in STKO
        self.midas_doc = midas_doc

        # the STKO document interface
        self.stko = stko

        # the grouped geometries (key = id in STKO)
        self._geoms : Dict[int, _geometry] = {}
        self._diaphram_retained_nodes_geom : MpcGeometry = None # the geometry for the retained nodes of the diaphragms

        # mapping of the MIDAS entities to the STKO entities
        self._vertex_map : Dict[int, _geometry_map_item] = {}
        self._frame_map : Dict[int, _geometry_map_item] = {} 
        self._area_map : Dict[int, _geometry_map_item] = {}
        self._diaphram_map : Dict[str, _interaction_map_item] = {}
        self._elastic_material_map : Dict[str, Tuple[MpcProperty, MpcProperty]] = {} # value = tuple (uniaxial material, NDMaterial)
        self._area_property_map : Dict[str, Union[MpcProperty, MpcElementProperty]] = {}
        self._frame_section_map : Dict[str, MpcProperty] = {}

        # keep track of stko indices for different entities
        # the default linear time series
        self._linear_time_series_id : int = 0
        # all sp constraints
        self._sp_ids : List[int] = []
        # all mp constraints
        self._mp_ids : List[int] = []
        # pattern-name to force condition ids
        self._pattern_load_ids : DefaultDict[str, List[int]] = defaultdict(list)
        self._pattern_eleload_ids : DefaultDict[str, List[int]] = defaultdict(list)

        # process
        try:
            self.stko.start()
            self._build_geometry_groups()
            self._build_geometries()
            self._build_interactions()
            self._build_local_axes()
            self._build_definitions()
            self._build_frame_sections()
            self._build_area_sections()
            self._build_conditions_restraints()
            self._build_conditions_diaphragms()
            self._build_conditions_masses()
            self._build_load_cases()
            self._build_finalization()
        except Exception as e:
            raise
        finally:
            self.stko.stop()

    # builds a graph of connected vertices
    def _make_vertex_groups(self) -> List[List[int]]:
        # the number of vertices
        num_vertices = len(self.midas_doc.vertices)
        # create the graph and compute connected vertices
        the_graph = graph(num_vertices)
        def _add_ele_to_graph(ele):
            for i,j in itertools.combinations(ele.nodes, 2):
                the_graph.add_edge(i, j)
        for _, ele in self.midas_doc.frames.items():
            _add_ele_to_graph(ele)
        for _, ele in self.midas_doc.areas.items():
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
                geom_floating_nodes.vertices[group[0]] = self.midas_doc.vertices[group[0]]
                continue
            # create a new geometry for the group
            geom = _geometry()
            geom.id = next_geom_id
            next_geom_id += 1
            # add the vertices to the geometry
            for node_id in group:
                geom.vertices[node_id] = self.midas_doc.vertices[node_id]
            # now add the frames and areas to the geometry
            for ele_id, ele in self.midas_doc.frames.items():
                for node_id in ele.nodes:
                    if node_id in group:
                        geom.frames[ele_id] = ele
                        break
            for ele_id, ele in self.midas_doc.areas.items():
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
            # map the MIDAS node id to the STKO vertex object
            V : Dict[int, FxOccShape] = {}
            for i, v in geom.vertices.items():
                V[i] = FxOccBuilder.makeVertex(v)
            # build frames:
            # map the MIDAS frame id to the STKO edge object
            # also map the tuple of vertices to the STKO a the edge object (assume no overlapping edges!)
            # (will be used to obtain boundary edges for faces)
            # keep track of edges that are not connected to faces!
            E : Dict[int, FxOccShape] = {} # edges from frame members
            EN : Dict[Tuple[int, int], int] = {} # vertex_map to E
            for i, f in geom.frames.items():
                # create the edge from existing vertices in V
                n1, n2 = f.nodes
                edge = FxOccBuilder.makeEdge(V[n1], V[n2])
                E[i] = edge
                EN[tuple(sorted((n1, n2)))] = i
            # build areas
            # map the MIDAS area id to the STKO face object
            EF : List[int] = [] # list of indices in E, that are also in faces
            EFN : Dict[Tuple[int, int], FxOccShape] = {} # maps a edge_key to an edge object generated for areas (not share with frames)
            F : Dict[int, FxOccShape] = {}
            FG_E : Dict[Tuple[int, int], int] = {} # maps and edge_key to a 0-based index for edges in faces (for face graphs)
            FG_E_F : DefaultDict[int, List[int]] = DefaultDict(list) # maps a face id to the edges in the face (edges indices in FE_E)
            for i, a in geom.areas.items():
                # create the face from existing edges in E
                # get the outer wire (only one supported now)
                outer_edges = []
                for j in range(len(a.nodes)):
                    n1 = a.nodes[j]
                    n2 = a.nodes[(j+1)%len(a.nodes)]
                    edge_key = tuple(sorted((n1, n2)))
                    # let's see if the edge is already in the list of edges from frames
                    eid = EN.get(edge_key, None)
                    if eid is None:
                        # not from frames, create a new edge or get existing one from other faces
                        edge = EFN.get(edge_key, None)
                        if edge is None:
                            # create a new edge from the vertices
                            edge = FxOccBuilder.makeEdge(V[n1], V[n2])
                            EFN[edge_key] = edge
                    else:
                        # the edge is also a frame, mark it as used in face
                        edge = E[eid]
                        EF.append(eid)
                    outer_edges.append(edge)
                    # for face graph...
                    fg_edge_index = FG_E.get(edge_key, None)
                    if fg_edge_index is None:
                        fg_edge_index = len(FG_E)
                        FG_E[edge_key] = fg_edge_index
                    FG_E_F[i].append(fg_edge_index)
                # create the face from the edges
                outer_wire = FxOccBuilder.makeWire(outer_edges)
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
                    # join floating edges in 1 or multiple wires
                    float_wires = FxOccFactoryCurve.makeWire(list(EFloat.values()))
                    # join floating faces in 1 or multiple shells
                    # float_faces = FxOccFactorySurfaces.makeShell(list(F.values()), self.etabs_doc.tolerance, False, True)
                    #float_face = FxOccBuilder.makeShell(list(F.values()))
                    fg_size = len(FG_E)
                    fg_graph = graph(fg_size)
                    for face_id, face_edges in FG_E_F.items():
                        for i,j in itertools.combinations(face_edges, 2):
                            fg_graph.add_edge(i, j)
                    fg_edge_groups = fg_graph.connected_components()
                    float_faces = []
                    for fg_edge_group in fg_edge_groups:
                        faces_for_shell = []
                        for face_id, face_edges in FG_E_F.items():
                            # if at least one edge is in the group, add the face
                            for face_edge in face_edges:
                                if face_edge in fg_edge_group:
                                    faces_for_shell.append(F[face_id])
                                    break
                        if len(faces_for_shell) == 0:
                            raise Exception('No faces found for shell')
                        if len(faces_for_shell) == 1:
                            float_faces.append(faces_for_shell[0])
                        else:
                            float_faces.append(FxOccBuilder.makeShell(faces_for_shell))
                    # make a compound
                    shape = FxOccBuilder.makeCompound(list(itertools.chain(float_wires, float_faces)))
            # create the STKO geometry
            geom_name = f'Geometry_{geom.id}'
            stko_geom = MpcGeometry(geom.id, geom_name, shape)
            self.stko.add_geometry(stko_geom)
            self._make_map(stko_geom)

    # maps an MIDAS _geometry to an STKO MpcGeometry (called by _build_geometries)
    def _make_map(self, geom : MpcGeometry):
        import time
        time_start = time.time()
        # map vertices
        tol = self.midas_doc.tolerance
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
        # now we can map the MIDAS entities to the STKO geometry/sub-geometry ids
        # using the MIDAS entities in this geometry group
        midas_geom = self._geoms[geom.id]
        # link the vertices to the STKO geometry
        for i, v in midas_geom.vertices.items():
            key = _vertex_key(v)
            subshape_id = unique_vertices.get(key, None)
            if subshape_id is None:
                raise Exception(f'Vertex {i} not found in geometry {geom.id}')
            self._vertex_map[i] = _geometry_map_item(geom, subshape_id)
        # link the frames to the STKO geometry
        for i, f in midas_geom.frames.items():
            # get the edge id from the map
            edge_vertices = f.nodes
            p1 = midas_geom.vertices[edge_vertices[0]]
            p2 = midas_geom.vertices[edge_vertices[1]]
            key = _edge_key(unique_vertices[_vertex_key(p1)], unique_vertices[_vertex_key(p2)])
            subshape_id = unique_edges.get(key, None)
            if subshape_id is None:
                raise Exception(f'Frame {i} not found in geometry {geom.id}')
            self._frame_map[i] = _geometry_map_item(geom, subshape_id)
        # link the areas to the STKO geometry
        for i, a in midas_geom.areas.items():
            # get the face id from the map
            face_vertices = a.nodes
            ekeys : List[Tuple[int,int]] = []
            for j in range(len(face_vertices)):
                n1 = a.nodes[j]
                n2 = a.nodes[(j+1)%len(a.nodes)]
                p1 = midas_geom.vertices[n1]
                p2 = midas_geom.vertices[n2]
                key = _edge_key(unique_vertices[_vertex_key(p1)], unique_vertices[_vertex_key(p2)])
                ekeys.append(key)
            subshape_id = unique_faces.get(_face_key(ekeys), None)
            if subshape_id is None:
                raise Exception(f'Area {i} not found in geometry {geom.id}')
            self._area_map[i] = _geometry_map_item(geom, subshape_id)
        
    # adds interactions to STKO
    def _build_interactions(self):
        # first, build a new compound geometry to hold the retained vertices of the diaphragms
        # created on-the-fly here
        all_retained_nodes = [] # for each diap, contains a vertex for the retained node
        all_retained_pos : Dict[str, Math.vec3] = {} # for each diap, contains the position of the retained node
        all_retained_ids : Dict[str, int] = {} # for each diap, contains the id of the retained vertex in the compound geometry
        all_constrained_ids : Dict[str, List[int]] = {} # for each diap, contains the ids of the constrained vertices
        for name, diap in self.midas_doc.diaphragms.items():
            # elevaltion of the diaphragm
            diap.Z
            # find all vertices with the same elevation excluding the ones in the released vertices
            slave_ids = [i for i,v in self.midas_doc.vertices.items() if abs(v.z - diap.Z) < self.midas_doc.tolerance and i not in self.midas_doc.diaphragm_released_vertices]
            all_constrained_ids[name] = slave_ids
            # find the retained vertex position as the average of the constrained vertices
            if len(slave_ids) == 0:
                raise Exception(f'Diaphragm {name} has no constrained vertices at elevation {diap.Z}')
            center = Math.vec3(0.0, 0.0, 0.0)
            for i in slave_ids:
                v = self._vertex_map.get(i, None)
                if v is None:
                    raise Exception(f'Diaphragm {name} constrained vertex {i} not found in geometry')
                center += v.geom.shape.vertexPosition(v.subshape_id)
            center /= len(slave_ids)
            all_retained_pos[name] = center
            all_retained_nodes.append(FxOccBuilder.makeVertex(center))
        # create the compound geometry for the retained nodes
        retained_geom = MpcGeometry(self.stko.new_geometry_id(), 'Diaphragm Retained Nodes', FxOccBuilder.makeCompound(all_retained_nodes))
        self.stko.add_geometry(retained_geom)
        self._diaphram_retained_nodes_geom = retained_geom
        # re-obtain the master vertex ids in the compound geometry
        for name, diap_retained_pos in all_retained_pos.items():
            retained_vertex_index = None
            for i in range(retained_geom.shape.getNumberOfSubshapes(MpcSubshapeType.Vertex)):
                v = retained_geom.shape.vertexPosition(i)
                if (v - diap_retained_pos).norm() < self.midas_doc.tolerance:
                    retained_vertex_index = i
                    break
            if retained_vertex_index is None:
                raise Exception(f'Diaphragm {diap.name} retained vertex not found in geometry')
            all_retained_ids[name] = retained_vertex_index
        # next interaction id in STKO
        next_id = self.stko.new_interaction_id()
        # process diaphragms
        for name, diap in self.midas_doc.diaphragms.items():
            # elevaltion of the diaphragm
            diap.Z
            # find all vertices with the same elevation excluding the ones in the released vertices
            constrained_ids = all_constrained_ids[name]
            # find the master vertex id
            retained_id = all_retained_ids[name]
            # make a new interaction
            interaction = MpcInteraction(next_id, f'Diaphragm {name}')
            next_id += 1
            # retained item
            interaction.items.masters.append(MpcInteractionItem(retained_geom, MpcSubshapeType.Vertex, retained_id))
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

    # builds the local axes and assign them to frame and area geometries in STKO to match
    # either the default or the user defined local axes in MIDAS
    def _build_local_axes(self):
        '''
        Frames:
        - STKO defaults:
            - vertical elements: local_y ~= -global_y (0,-1,0)
            - others: local_z ~= global_z (0,0,1)
        - MIDAS defaults:
            - vertical elements: local_y ~= global_x (1,0,0)
                                 local_z ~= global_x (1,0,0)
            - others: local_y ~= global_z (0,0,1)
        Areas:
        - STKO and MIDAS defaults are the same
        '''
        # tolerance for normalized axis components
        tol = 1.0e-3

        # utilities for reversing subgeoms
        def _check_reverse_subshape(stype:MpcSubshapeType, midas_id:int, axis:Math.vec3):
            # compare to geom
            stko_dz = Math.vec3()
            stko_dx = Math.vec3()
            stko_dy = Math.vec3()
            if stype == MpcSubshapeType.Edge:
                edge_data = self._frame_map.get(midas_id, None)
                if edge_data is None:
                    raise Exception(f'Frame {midas_id} not found in geometry')
                edge_data.geom.getLocalAxesOnEdge(edge_data.subshape_id, 0.0, stko_dx, stko_dy, stko_dz)
                if axis.dot(stko_dx) < 0.0:
                    #print(f'Edge {midas_id} has inconsistent orientation with the local axes. Reversing it...')
                    if edge_data.geom.id == 1 and edge_data.subshape_id == 2077:
                        print('reversing edge 2077 in geometry 1 from midas_id {} ({})'.format(midas_id, stype))
                    edge_data.geom.shape.toggleEdgeReversedFlag(edge_data.subshape_id)
            elif stype == MpcSubshapeType.Face:
                area_data = self._area_map.get(midas_id, None)
                if area_data is None:
                    raise Exception(f'Area {midas_id} not found in geometry')
                area_data.geom.getLocalAxesOnFace(area_data.subshape_id, 0.0, 0.0, stko_dx, stko_dy, stko_dz)
                if axis.dot(stko_dz) < 0.0:
                    #print(f'Area {midas_id} has inconsistent orientation with the local axes. Reversing it...')
                    area_data.geom.shape.toggleFaceReversedFlag(area_data.subshape_id)
            else:
                raise Exception(f'Unsupported subshape type {stype}')

        # first, we need to reverse sub-edges and sub-faces in STKO if their director
        # does not match the director from the source frame or area in MIDAS.
        for midas_id, elem in self.midas_doc.frames.items():
            dx = MIDASLocalAxesConvention.get_director([self.midas_doc.vertices[i] for i in elem.nodes])
            _check_reverse_subshape(MpcSubshapeType.Edge, midas_id, dx)
        for midas_id, elem in self.midas_doc.areas.items():
            dz = MIDASLocalAxesConvention.get_director([self.midas_doc.vertices[i] for i in elem.nodes])
            _check_reverse_subshape(MpcSubshapeType.Face, midas_id, dz)

        # next local axes id in STKO
        next_locax_id = self.stko.new_local_axes_id()
        # maps a tuple of (x,y,z) to the local axis id
        # the key is a tuple of integers obtained as the rounded values of the axis components
        locax_id_map : Dict[Tuple[Tuple[int,int,int], Tuple[int,int,int], Tuple[int,int,int]], int] = {}
        def _make_key(x:Math.vec3, y:Math.vec3, z:Math.vec3) -> Tuple[Tuple[int,int,int], Tuple[int,int,int], Tuple[int,int,int]]:
            return (
                tuple(int(round(v/tol)) for v in x), 
                tuple(int(round(v/tol)) for v in y), 
                tuple(int(round(v/tol)) for v in z))
        # maps the local axis id to the local axis object in STKO
        locax_map : Dict[int, MpcLocalAxes] = {}
        # maps a frame to the local axis id
        frame_to_locax_id_map : Dict[int, int] = {}
        # maps an area to the local axis id
        area_to_locax_id_map : Dict[int, int] = {}

        # define local axes for frames as per MIDAS
        stko_dx = Math.vec3()
        stko_dy = Math.vec3()
        stko_dz = Math.vec3()
        for midas_id, elem in self.midas_doc.frames.items():
            # get local axes from MIDAS defaults
            T_midas = MIDASLocalAxesConvention.get_local_axes([self.midas_doc.vertices[i] for i in elem.nodes])
            # rotate about the local x if user angle is defined
            dx = T_midas.col(0)
            dy = T_midas.col(1)
            dz = T_midas.col(2)
            if abs(elem.angle) > 1.0e-10:
                q = Math.quaternion.fromAxisAngle(dx, elem.angle/180.0*math.pi)
                dy = q.rotate(dy).normalized()
                dz = q.rotate(dz).normalized()
            # get the default from STKO
            area_data = self._frame_map[midas_id]
            area_data.geom.getLocalAxesOnEdge(area_data.subshape_id, 0.0, stko_dx, stko_dy, stko_dz)
            # check similarity
            if stko_dx.dot(dx) < 0.99:
                raise Exception(f'Frame {midas_id} has inconsistent local axes with the MIDAS defaults. Please check the frame orientation in MIDAS.')
            # now we are sure the main axis is aligned with the frame
            # next step: if the local y (or z) is not aligned with the one from STKO's defaults,
            # we need to make a local axis object in STKO
            if dy.dot(stko_dy) < 0.99:
                # map the frame to the local axis
                locax_key = _make_key(dx, dy, dz)
                locax_id = locax_id_map.get(locax_key, None)
                if locax_id is None:
                    locax_id = next_locax_id
                    locax_id_map[locax_key] = next_locax_id
                    next_locax_id += 1
                frame_to_locax_id_map[midas_id] = locax_id

        # define local axes for areas as per MIDAS
        for midas_id, elem in self.midas_doc.areas.items():
            # get local axes from MIDAS defaults
            T_midas = MIDASLocalAxesConvention.get_local_axes([self.midas_doc.vertices[i] for i in elem.nodes])
            # rotate about the local z if user angle is defined
            dx = T_midas.col(0)
            dy = T_midas.col(1)
            dz = T_midas.col(2)
            if abs(elem.angle) > 1.0e-10:
                q = Math.quaternion.fromAxisAngle(dz, elem.angle/180.0*math.pi)
                dy = q.rotate(dy).normalized()
                dx = q.rotate(dx).normalized()
            # get the default from STKO
            area_data = self._area_map[midas_id]
            area_data.geom.getLocalAxesOnFace(area_data.subshape_id, 0.0, 0.0, stko_dx, stko_dy, stko_dz)
            # check similarity
            if stko_dz.dot(dz) < 0.99:
                raise Exception(f'Area {midas_id} has inconsistent local axes with the MIDAS defaults. Please check the frame orientation in MIDAS.')
            # now we are sure the main axis is aligned with the area normal
            # next step: if the local y (or zx) is not aligned with the one from STKO's defaults,
            # we need to make a local axis object in STKO
            if dy.dot(stko_dy) < 0.99:
                # map the frame to the local axis
                locax_key = _make_key(dx, dy, dz)
                locax_id = locax_id_map.get(locax_key, None)
                if locax_id is None:
                    locax_id = next_locax_id
                    locax_id_map[locax_key] = next_locax_id
                    next_locax_id += 1
                area_to_locax_id_map[midas_id] = locax_id
        
        # now generate local axes in STKO
        for (x,y,_), id in locax_id_map.items():
            # create the local axis object
            p1 = self.midas_doc.bbox.minPoint - (self.midas_doc.bbox.maxPoint - self.midas_doc.bbox.minPoint).norm()*0.1*Math.vec3(1,1,1)
            p2 = p1 + Math.vec3(float(x[0])*tol, float(x[1])*tol, float(x[2])*tol)
            p3 = p1 + Math.vec3(float(y[0])*tol, float(y[1])*tol, float(y[2])*tol)
            locax = MpcLocalAxes(id, f'Local Axes {id}', MpcLocalAxesType.Rectangular, p1, p2, p3)
            # add it to the document
            self.stko.add_local_axes(locax)
            # map it
            locax_map[id] = locax
        
        # assign to frames
        for frame_id, locax_id in frame_to_locax_id_map.items():
            # get the geometry and subshape id
            frame_data = self._frame_map.get(frame_id, None)
            if frame_data is None:
                raise Exception(f'Frame {frame_id} not found in geometry')
            # set the local axes id
            frame_data.geom.assign(locax_map[locax_id], frame_data.subshape_id, MpcSubshapeType.Edge)
        
        # assign to areas
        for area_id, locax_id in area_to_locax_id_map.items():
            # get the geometry and subshape id
            area_data = self._area_map.get(area_id, None)
            if area_data is None:
                raise Exception(f'Area {area_id} not found in geometry')
            # set the local axes id
            area_data.geom.assign(locax_map[locax_id], area_data.subshape_id, MpcSubshapeType.Face)

    # build definitions
    def _build_definitions(self):
        # 1. The default linear time series
        definition = MpcDefinition()
        definition.id = self.stko.new_definition_id()
        definition.name = 'Linear Time Series'
        # define xobject
        meta = self.stko.doc.metaDataDefinition('timeSeries.Linear')
        xobj = MpcXObject.createInstanceOf(meta)
        definition.XObject = xobj
        # add the definition to the document
        self.stko.add_definition(definition)
        definition.commitXObjectChanges()
        # track the time series id
        self._linear_time_series_id = definition.id

    # builds the frame sections in STKO
    def _build_frame_sections(self):
        '''
        In midas elements have 1) a material for E and v, 2) a section and 3) section modifiers.
        For each frame, we have to compute combinations of material + section + modifiers
        '''

        # create 1 element property for all frames
        ele_prop = MpcElementProperty()
        ele_prop.id = self.stko.new_element_property_id()
        ele_prop.name = 'Elastic Beam Element'
        meta = self.stko.doc.metaDataElementProperty('beam_column_elements.elasticBeamColumn')
        xobj = MpcXObject.createInstanceOf(meta)
        ele_prop.XObject = xobj
        self.stko.add_element_property(ele_prop)
        ele_prop.commitXObjectChanges()

        # map frame to the index in the section modifier list
        map_frame_mod : Dict[int, int] = {}
        for mod_id, mod in enumerate(self.midas_doc.section_scale_factors):
            for ele_id in mod.elements:
                map_frame_mod[ele_id] = mod_id
        print_mod = len(map_frame_mod) > 1

        # for each frame that is not a link, we need to create a key = tuple(material, section, modifier)
        map_sec_ele : DefaultDict[Tuple[int, int, int], List[int]] = defaultdict(list)
        for ele_id, ele in self.midas_doc.frames.items():
            # skip links
            if ele.link is not None:
                continue 
            # get material and section
            mat = ele.mat
            sec = ele.sec
            mod = map_frame_mod.get(ele_id, -1)
            key = (mat, sec, mod)
            # add the frame to the map
            map_sec_ele[key].append(ele_id)

        # now we can build the sections
        for (mat_id, sec_id, mod_id), ele_ids in map_sec_ele.items():
            # get the material
            mat = self.midas_doc.elastic_materials.get(mat_id, None)
            if mat is None:
                raise Exception(f'Material {mat_id} not found for frame section with material {mat_id}, section {sec_id}, modifier {mod_id}')
            # get the section
            sec = self.midas_doc.sections.get(sec_id, None)
            if sec is None:
                raise Exception(f'Section {sec_id} not found for frame section with material {mat_id}, section {sec_id}, modifier {mod_id}')
            # get the modifier
            mod = self.midas_doc.section_scale_factors[mod_id] if mod_id >= 0 else None
            # generate the frame section
            prop = MpcProperty()
            prop.id = self.stko.new_physical_property_id()
            if print_mod:
                prop.name = f'Frame Section {mat.name} - {sec.name} - Modifier {mod_id if mod else "None"}'
            else:
                prop.name = f'Frame Section {mat.name} - {sec.name}'
            meta = self.stko.doc.metaDataPhysicalProperty('sections.Elastic')
            xobj = MpcXObject.createInstanceOf(meta)
            # common properties
            xobj.getAttribute('E').quantityScalar.value = mat.E
            xobj.getAttribute('G').quantityScalar.value = mat.E/(2*(1+mat.poiss))
            if mod is not None:
                xobj.getAttribute('A_modifier').real = max(1.0e-3, mod.A)
                xobj.getAttribute('Asy_modifier').real = max(1.0e-3, mod.Asy)
                xobj.getAttribute('Asz_modifier').real = max(1.0e-3, mod.Asz)
                xobj.getAttribute('Iyy_modifier').real = max(1.0e-3, mod.Iyy)
                xobj.getAttribute('Izz_modifier').real = max(1.0e-3, mod.Izz)
                xobj.getAttribute('J_modifier').real = max(1.0e-3, mod.Ixx)
            xobj.getAttribute('Y/section_offset').quantityScalar.value = sec.offset_y
            xobj.getAttribute('Z/section_offset').quantityScalar.value = sec.offset_z
            xobj.getAttribute('Shear Deformable').boolean = True # make it shear
            # define the section
            if sec.type == section.shape_type.SB:
                H, B = sec.shape_info[:2]
                stko_section = MpcBeamSection(
                    MpcBeamSectionShapeType.Rectangular,
                    'section_Box', 'user', self.midas_doc.units[1].lower(),
                    [H, B]
                )
            elif sec.type == section.shape_type.L:
                H, B, tw, tf = sec.shape_info[:4]
                stko_section = MpcBeamSection(
                    MpcBeamSectionShapeType.LU,
                    'section_L_U', 'user', self.midas_doc.units[1].lower(),
                    [H, B, tw, tf]
                )
            elif sec.type == section.shape_type.T:
                H, B, tw, tf = sec.shape_info[:4]
                stko_section = MpcBeamSection(
                    MpcBeamSectionShapeType.T,
                    'section_T', 'user', self.midas_doc.units[1].lower(),
                    [H, B, tw, tf]
                )
            else:
                raise Exception(f'Unsupported section type {sec.type} for frame section with material {mat_id}, section {sec_id}, modifier {mod_id}')
            # set the section
            xobj.getAttribute('Section').customObject = stko_section
            # set the xobject
            prop.XObject = xobj
            self.stko.add_physical_property(prop)
            prop.commitXObjectChanges()

            # now we can assign this section to all target frames
            geom_sub_map : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
            for ele_id in ele_ids:
                # get the geometry and subshape id
                frame_data = self._frame_map.get(ele_id, None)
                if frame_data is None:
                    raise Exception(f'Frame {ele_id} not found in geometry')
                geom_sub_map[frame_data.geom].append(frame_data.subshape_id)
            for geom, subshape_ids in geom_sub_map.items():
                # assign the section to the frame
                geom.assign(prop, subshape_ids, MpcSubshapeType.Edge)
                # assign also the element property
                geom.assign(ele_prop, subshape_ids, MpcSubshapeType.Edge)

    # builds the area sections in STKO
    def _build_area_sections(self):

        # create 1 element property for all areas
        ele_prop_q4 = MpcElementProperty()
        ele_prop_q4.id = self.stko.new_element_property_id()
        ele_prop_q4.name = 'Elastic Shell Element'
        meta = self.stko.doc.metaDataElementProperty('shell.ASDShellQ4')
        xobj = MpcXObject.createInstanceOf(meta)
        ele_prop_q4.XObject = xobj
        self.stko.add_element_property(ele_prop_q4)
        ele_prop_q4.commitXObjectChanges()

        # create 1 element property for all areas
        ele_prop_t3 = MpcElementProperty()
        ele_prop_t3.id = self.stko.new_element_property_id()
        ele_prop_t3.name = 'Elastic Shell Element'
        meta = self.stko.doc.metaDataElementProperty('shell.ASDShellT3')
        xobj = MpcXObject.createInstanceOf(meta)
        ele_prop_t3.XObject = xobj
        self.stko.add_element_property(ele_prop_t3)
        ele_prop_t3.commitXObjectChanges()

        # map areas to the index in the thickness modifier list
        map_area_mod : Dict[int, int] = {}
        for mod_id, mod in enumerate(self.midas_doc.thickness_scale_factors):
            for ele_id in mod.elements:
                map_area_mod[ele_id] = mod_id
        print_mod = len(map_area_mod) > 1
        
        # for each area we need to create a key = tuple(material, thickness, modifier)
        map_thick_ele : DefaultDict[Tuple[int, int, int], List[int]] = defaultdict(list)
        for ele_id, ele in self.midas_doc.areas.items():
            # get material and section
            mat = ele.mat
            sec = ele.sec
            mod = map_area_mod.get(ele_id, -1)
            key = (mat, sec, mod)
            # add the frame to the map
            map_thick_ele[key].append(ele_id)

        # now we can build the sections
        for (mat_id, sec_id, mod_id), ele_ids in map_thick_ele.items():
            # get the material
            mat = self.midas_doc.elastic_materials.get(mat_id, None)
            if mat is None:
                raise Exception(f'Material {mat_id} not found for area section with material {mat_id}, section {sec_id}, modifier {mod_id}')
            # get the section
            sec = self.midas_doc.thicknesses.get(sec_id, None)
            if sec is None:
                raise Exception(f'Thickness {sec_id} not found for area section with material {mat_id}, section {sec_id}, modifier {mod_id}')
            if sec.offset != 0.0:
                raise Exception(f'Thickness {sec_id} has an offset {sec.offset} which is not supported in STKO. Please remove the offset from the thickness in MIDAS.')
            # get the modifier
            mod = self.midas_doc.thickness_scale_factors[mod_id] if mod_id >= 0 else None
            # now we can compute the stiffnesses of the plate and create a layered shell if there is an offset
            thickness = sec.in_thick
            oop_thickness = sec.out_thick
            if mod is not None:
                thickness *= mod.ip_mod
                oop_thickness *= mod.oop_mod
            oop_scale = oop_thickness / thickness
            # generate the elastic shell section
            prop = MpcProperty()
            prop.id = self.stko.new_physical_property_id()
            if print_mod:
                prop.name = f'Area Section {mat.name} - {sec.name} - Modifier {mod_id if mod else "None"}'
            else:
                prop.name = f'Area Section {mat.name} - {sec.name}'
            meta = self.stko.doc.metaDataPhysicalProperty('sections.ElasticMembranePlateSection')
            xobj = MpcXObject.createInstanceOf(meta)
            # common properties
            xobj.getAttribute('E').quantityScalar.value = mat.E
            xobj.getAttribute('nu').real = mat.poiss
            xobj.getAttribute('h').quantityScalar.value = thickness
            xobj.getAttribute('Ep_mod').real = oop_scale
            # set the xobject
            prop.XObject = xobj
            self.stko.add_physical_property(prop)
            prop.commitXObjectChanges()

            # now we can assign this section to all target frames
            geom_sub_map : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
            geom_sub_map_t3 : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
            geom_sub_map_q4 : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
            for ele_id in ele_ids:
                # get the geometry and subshape id
                area_data = self._area_map.get(ele_id, None)
                if area_data is None:
                    raise Exception(f'Area {ele_id} not found in geometry')
                geom_sub_map[area_data.geom].append(area_data.subshape_id) # for physical property
                if len(self.midas_doc.areas[ele_id].nodes) == 3:
                    geom_sub_map_t3[area_data.geom].append(area_data.subshape_id) # for t3 elements
                else:
                    geom_sub_map_q4[area_data.geom].append(area_data.subshape_id) # for q4 elements
            # assign the section to the area
            for geom, subshape_ids in geom_sub_map.items():
                geom.assign(prop, subshape_ids, MpcSubshapeType.Face)
            # assign also the element property
            for geom, subshape_ids in geom_sub_map_t3.items():
                geom.assign(ele_prop_t3, subshape_ids, MpcSubshapeType.Face)
            for geom, subshape_ids in geom_sub_map_q4.items():
                geom.assign(ele_prop_q4, subshape_ids, MpcSubshapeType.Face)

    # builds the conditions for restraints in STKO
    def _build_conditions_restraints(self):
        # utils
        def _make_xobj_gen():
            # define xobject
            meta = self.stko.doc.metaDataCondition('Constraints.sp.fix')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('Dimension').string = '3D'
            xobj.getAttribute('ModelType').string = 'U-R (Displacement+Rotation)'
            xobj.getAttribute('2D').boolean = False
            xobj.getAttribute('3D').boolean = True
            xobj.getAttribute('U (Displacement)').boolean = False
            xobj.getAttribute('U-P (Displacement+Pressure)').boolean = False
            xobj.getAttribute('U-R (Displacement+Rotation)').boolean = True
            return xobj
        # group restraints by type
        # key = tuple(int * 6), value = list of vertex ids
        restraints : DefaultDict[Tuple[int,int,int,int,int,int], List[int]] = DefaultDict(list)
        for i, r in self.midas_doc.constraints.items():
            restraints[r].append(i)
        # process each group of restraints
        for rtype, asn_nodes in restraints.items():
            # create a new condition
            condition = MpcCondition()
            condition.id = self.stko.new_condition_id()
            condition.name = f'Fix {rtype}'
            # define xobject
            xobj = _make_xobj_gen()
            for dof_flag, dof_label in zip(rtype, _globals.fix_labels):
                xobj.getAttribute(dof_label).boolean = dof_flag == 1
            condition.XObject = xobj
            # add the assigned nodes to the condition
            for i in asn_nodes:
                # get the geometry and subshape id
                vertex_data = self._vertex_map.get(i, None)
                if vertex_data is None:
                    raise Exception(f'Vertex {i} not found in geometry')
                sset = MpcConditionIndexedSubSet()
                sset.vertices.append(vertex_data.subshape_id)
                condition.assignTo(vertex_data.geom, sset)
            # add the condition to the document
            self.stko.add_condition(condition)
            condition.commitXObjectChanges()
            # track the condition id
            self._sp_ids.append(condition.id)
        # fix master nodes for rigid diaphragm
        if self._diaphram_retained_nodes_geom is not None:
            # create a new condition for the rigid diaphragm retained nodes
            condition = MpcCondition()
            condition.id = self.stko.new_condition_id()
            condition.name = 'Fix Rigid Diaphragm Retained Nodes'
            # define xobject
            meta = self.stko.doc.metaDataCondition('Constraints.mp.rigidDiaphragm')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('perpDirn').integer = 3
            # define xobject
            xobj = _make_xobj_gen()
            for dof_flag, dof_label in zip((0,0,1,1,1,0), _globals.fix_labels):
                xobj.getAttribute(dof_label).boolean = dof_flag == 1
            condition.XObject = xobj
            # add the assigned nodes to the condition
            for i in range(self._diaphram_retained_nodes_geom.shape.getNumberOfSubshapes(MpcSubshapeType.Vertex)):
                sset = MpcConditionIndexedSubSet()
                sset.vertices.append(i)
                condition.assignTo(self._diaphram_retained_nodes_geom, sset)
            # add the condition to the document
            self.stko.add_condition(condition)
            condition.commitXObjectChanges()
            # track the condition id
            self._mp_ids.append(condition.id)

    # builds the conditions for rigid diaphragms in STKO
    def _build_conditions_diaphragms(self):
        # build 1 condition for all diaphragms
        condition = MpcCondition()
        condition.id = self.stko.new_condition_id()
        condition.name = f'Diaphragm'
        # define xobject
        meta = self.stko.doc.metaDataCondition('Constraints.mp.rigidDiaphragm')
        xobj = MpcXObject.createInstanceOf(meta)
        xobj.getAttribute('perpDirn').integer = 3
        condition.XObject = xobj
        for name, _ in self.midas_doc.diaphragms.items():
            # add the assigned interactions to the condition
            interaction = self._diaphram_map.get(name, None)
            if interaction is None:
                raise Exception(f'Interaction for Diaphragm {name} not found in STKO')
            condition.assignTo(interaction.interaction)
         # add the condition to the document
        self.stko.add_condition(condition)
        condition.commitXObjectChanges()
        # track the condition id
        self._mp_ids.append(condition.id)

    # builds the conditions for joint masses in STKO
    def _build_conditions_masses(self):
        for (mx, my, mz), node_ids in self.midas_doc.masses.items():
            # create a new condition
            condition = MpcCondition()
            condition.id = self.stko.new_condition_id()
            condition.name = 'Joint Mass {:.2g} {:.2g} {:.2g}'.format(mx, my, mz)
            # define xobject
            meta = self.stko.doc.metaDataCondition('Mass.NodeMass')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('mass').quantityVector3.value = Math.vec3(mx, my, mz)
            condition.XObject = xobj
            # add the assigned nodes to the condition
            # group nodes by geometry
            geom_node_map : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
            for i in node_ids:
                # get the geometry and subshape id
                vertex_data = self._vertex_map.get(i, None)
                if vertex_data is None:
                    raise Exception(f'Vertex {i} not found in geometry')
                geom_node_map[vertex_data.geom].append(vertex_data.subshape_id)
            for geom, vertices in geom_node_map.items():
                sset = MpcConditionIndexedSubSet()
                for iv in vertices:
                    sset.vertices.append(iv)
                condition.assignTo(geom, sset)
            # add the condition to the document
            self.stko.add_condition(condition)
            condition.commitXObjectChanges()

    # builds the self-weight condition in STKO
    def _build_self_weight(self, lc:load_case):
        if lc.self_weight is None:
            # no self-weight defined, skip
            return
        # create a self-weight condition
        g = self.midas_doc.struct_type.grav
        for mass, node_ids in self.midas_doc.masses_dens.items():
            sw = Math.vec3(0.0, 0.0, -g * mass)
            # make condition only if self-weight is not zero
            if sw.norm() < 1.0e-10:
                continue
            condition = MpcCondition()
            condition.id = self.stko.new_condition_id()
            condition.name = f'Self Weight ({sw.z:.2g})'
            # define xobject
            meta = self.stko.doc.metaDataCondition('Loads.Force.NodeForce')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('F').quantityVector3.value = sw
            condition.XObject = xobj
            # add the assigned nodes to the condition
            # group nodes by geometry
            geom_node_map : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
            for i in node_ids:
                # get the geometry and subshape id
                vertex_data = self._vertex_map.get(i, None)
                if vertex_data is None:
                    raise Exception(f'Vertex {i} not found in geometry')
                geom_node_map[vertex_data.geom].append(vertex_data.subshape_id)
            for geom, vertices in geom_node_map.items():
                sset = MpcConditionIndexedSubSet()
                for iv in vertices:
                    sset.vertices.append(iv)
                condition.assignTo(geom, sset)
            # add the condition to the document
            self.stko.add_condition(condition)
            condition.commitXObjectChanges()
            # track the load id for the load pattern
            self._pattern_load_ids[lc.name].append(condition.id)
                
    # builds the nodal loads in STKO
    def _build_nodal_loads(self, lc:load_case):
        for node_load, node_ids in lc.nodal_loads.items():
            F = Math.vec3(*node_load.value[:3])
            M = Math.vec3(*node_load.value[3:])
            for load_vec, meta_name, label in zip((F, M), ('NodeForce', 'NodeCouple'), ('F', 'M')):
                # skip null
                if load_vec.norm() < 1.0e-10:
                    continue
                # create a new condition
                condition = MpcCondition()
                condition.id = self.stko.new_condition_id()
                condition.name = f'{meta_name} {node_load}'
                # define xobject
                meta = self.stko.doc.metaDataCondition(f'Loads.Force.{meta_name}')
                xobj = MpcXObject.createInstanceOf(meta)
                xobj.getAttribute(label).quantityVector3.value = load_vec
                condition.XObject = xobj
                # add the assigned nodes to the condition
                geom_node_map : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
                for i in node_ids:
                    # get the geometry and subshape id
                    vertex_data = self._vertex_map.get(i, None)
                    if vertex_data is None:
                        raise Exception(f'Vertex {i} not found in geometry')
                    geom_node_map[vertex_data.geom].append(vertex_data.subshape_id)
                for geom, vertices in geom_node_map.items():
                    sset = MpcConditionIndexedSubSet()
                    for iv in vertices:
                        sset.vertices.append(iv)
                    condition.assignTo(geom, sset)
                # add the condition to the document
                self.stko.add_condition(condition)
                condition.commitXObjectChanges()
                # track the load id for the load pattern
                self._pattern_load_ids[lc.name].append(condition.id)

    # builds the beam loads in STKO
    def _build_beam_loads(self, lc:load_case):
        for ele_load, ele_ids in lc.beam_loads.items():
            # create a new condition
            condition = MpcCondition()
            condition.id = self.stko.new_condition_id()
            condition.name = f'Beam Load ({ele_load.value[0]:.2g} {ele_load.value[1]:.2g} {ele_load.value[2]:.2g})'
            # define xobject
            meta = self.stko.doc.metaDataCondition('Loads.eleLoad.eleLoad_beamUniform')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('Dimension').string = '3D'
            xobj.getAttribute('2D').boolean = False
            xobj.getAttribute('3D').boolean = True
            xobj.getAttribute('use_Wx').boolean = True
            xobj.getAttribute('Wx').real = ele_load.value[0]
            xobj.getAttribute('Wy').real = ele_load.value[1]
            xobj.getAttribute('Wz').real = ele_load.value[2]
            if ele_load.local:
                xobj.getAttribute('Orientation').string = 'Local'
                xobj.getAttribute('Global').boolean = False
            else:
                xobj.getAttribute('Orientation').string = 'Global'
                xobj.getAttribute('Global').boolean = True
            condition.XObject = xobj
            # add the assigned elements to the condition
            geom_ele_map : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
            for i in ele_ids:
                # get the geometry and subshape id
                frame_data = self._frame_map.get(i, None)
                if frame_data is None:
                    raise Exception(f'Frame {i} not found in geometry')
                geom_ele_map[frame_data.geom].append(frame_data.subshape_id)
            for geom, edges in geom_ele_map.items():
                sset = MpcConditionIndexedSubSet()
                for ie in edges:
                    sset.edges.append(ie)
                condition.assignTo(geom, sset)
            # add the condition to the document
            self.stko.add_condition(condition)
            condition.commitXObjectChanges()
            # track the load id for the load pattern
            self._pattern_eleload_ids[lc.name].append(condition.id)

    # builds the pressure loads in STKO
    def _build_pressure_loads(self, lc:load_case):
        for press, ele_ids in lc.pressure_loads.items():
            # create a new condition
            condition = MpcCondition()
            condition.id = self.stko.new_condition_id()
            condition.name = f'Pressure Load ({press.value[0]:.2g}, {press.value[1]:.2g}, {press.value[2]:.2g})'
            # define xobject
            meta = self.stko.doc.metaDataCondition('Loads.Force.FaceForce')
            xobj = MpcXObject.createInstanceOf(meta)
            if press.local:
                xobj.getAttribute('Orientation').string = 'Local'
                xobj.getAttribute('Global').boolean = False
            else:
                xobj.getAttribute('Orientation').string = 'Global'
                xobj.getAttribute('Global').boolean = True
            xobj.getAttribute('F').quantityVector3.value = Math.vec3(*press.value[:3])
            condition.XObject = xobj
            # add the assigned elements to the condition
            geom_ele_map : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
            for i in ele_ids:
                # get the geometry and subshape id
                area_data = self._area_map.get(i, None)
                if area_data is None:
                    raise Exception(f'Area {i} not found in geometry')
                geom_ele_map[area_data.geom].append(area_data.subshape_id)
            for geom, faces in geom_ele_map.items():
                sset = MpcConditionIndexedSubSet()
                for ie in faces:
                    sset.faces.append(ie)
                condition.assignTo(geom, sset)
            # add the condition to the document
            self.stko.add_condition(condition)
            condition.commitXObjectChanges()
            # track the load id for the load pattern
            self._pattern_load_ids[lc.name].append(condition.id)

    # builds the floor loads in STKO
    def _build_floor_loads(self, lc:load_case):
        # build floor load as NodeForce (not Moments)
        for floor_load, node_ids in lc.floor_loads.items():
            # create a new condition
            condition = MpcCondition()
            condition.id = self.stko.new_condition_id()
            condition.name = f'Floor Load ({floor_load.value[0]:.2g}, {floor_load.value[1]:.2g}, {floor_load.value[2]:.2g})'
            # define xobject
            meta = self.stko.doc.metaDataCondition('Loads.Force.NodeForce')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('F').quantityVector3.value = Math.vec3(*floor_load.value)
            condition.XObject = xobj
            # add the assigned nodes to the condition
            geom_node_map : DefaultDict[MpcGeometry, List[int]] = defaultdict(list)
            for i in node_ids:
                # get the geometry and subshape id
                vertex_data = self._vertex_map.get(i, None)
                if vertex_data is None:
                    raise Exception(f'Vertex {i} not found in geometry')
                geom_node_map[vertex_data.geom].append(vertex_data.subshape_id)
            for geom, vertices in geom_node_map.items():
                sset = MpcConditionIndexedSubSet()
                for iv in vertices:
                    sset.vertices.append(iv)
                condition.assignTo(geom, sset)
            # add the condition to the document
            self.stko.add_condition(condition)
            condition.commitXObjectChanges()
            # track the load id for the load pattern
            self._pattern_load_ids[lc.name].append(condition.id)

    # builds the load cases in STKO
    def _build_load_cases(self):
        for lc_name, lc in self.midas_doc.load_cases.items():
            self._build_self_weight(lc)
            self._build_nodal_loads(lc)
            self._build_beam_loads(lc)
            self._build_pressure_loads(lc)
            self._build_floor_loads(lc)
        # now build the load patterns
        for lc_name, _ in self.midas_doc.load_cases.items():
            # create a new Analysis Step
            step = MpcAnalysisStep()
            step.id = self.stko.new_analysis_step_id()
            step.name = f'Load Pattern {lc_name}'
            # define xobject
            meta = self.stko.doc.metaDataAnalysisStep('Patterns.addPattern.loadPattern')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('tsTag').index = self._linear_time_series_id
            # add conditions
            # - load
            load_ids = self._pattern_load_ids[lc_name]
            if len(load_ids) > 0:
                xobj.getAttribute('load').indexVector = load_ids
            # - eleLoad
            ele_load_ids = self._pattern_eleload_ids[lc_name]
            if len(ele_load_ids) > 0:
                xobj.getAttribute('eleLoad').indexVector = ele_load_ids
            # scale
            xobj.getAttribute('-fact').boolean = True
            xobj.getAttribute('cFactor').real = 1.0
            # set xobj
            step.XObject = xobj
            # add the analysis step to the document
            self.stko.add_analysis_step(step)
            step.commitXObjectChanges()

    # do some extra stuff to finalize the import
    def _build_finalization(self):
        self.stko.assign_custom_colors()
        # TODO: keep dialog open...
        # when the dialog is closed, we can ask the user to save the document