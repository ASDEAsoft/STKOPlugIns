from PyMpc import *
from MidasGenImporter.graph import graph
from MidasGenImporter.document import *
from MidasGenImporter.stko_interface import stko_interface
from MidasGenImporter.units_utils import unit_system
import itertools
from typing import List, Dict, Tuple, DefaultDict, Union
import math
from collections import defaultdict

class _globals:
    fix_labels = ('Ux', 'Uy', 'Uz', 'Rx', 'Ry', 'Rz/3D')
    uniform_excitation_dir_map = {
        'U1' : 'dx',
        'U2' : 'dy',
        'U3' : 'dz',
        'R1' : 'Rx',
        'R2' : 'Ry',
        'R3' : 'Rz',
    }

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

        # mapping of the ETABS entities to the STKO entities
        self._vertex_map : Dict[int, _geometry_map_item] = {}
        self._frame_map : Dict[int, _geometry_map_item] = {} 
        self._area_map : Dict[int, _geometry_map_item] = {}
        self._diaphram_map : Dict[str, _interaction_map_item] = {}
        self._elastic_material_map : Dict[str, Tuple[MpcProperty, MpcProperty]] = {} # value = tuple (uniaxial material, NDMaterial)
        self._area_property_map : Dict[str, Union[MpcProperty, MpcElementProperty]] = {}
        self._frame_section_map : Dict[str, MpcProperty] = {}
        self._frame_nonlinear_hinge_map : Dict[str, MpcProperty] = {}
        self._frame_nonlinear_stiff_mat : MpcProperty = None

        # keep track of stko indices for different entities
        # the default linear time series
        self._linear_time_series_id : int = 0
        # all sp constraints
        self._sp_ids : List[int] = []
        # all mp constraints
        self._mp_ids : List[int] = []
        # pattern-name to force condition ids
        self._pattern_load_ids : DefaultDict[str, List[int]] = defaultdict(list)
        # time history name to condition ids
        self._th_ids : Dict[str, int] = {}

        # process
        try:
            self.stko.start()
            self._build_geometry_groups()
            self._build_geometries()
            self._build_interactions()
            self._build_local_axes()
            self._build_definitions()
            self._build_elastic_materials()
            self._build_area_sections()
            self._assign_area_sections()
            self._build_frame_sections()
            self._build_frame_nonlinear_hinges()
            self._assign_frame_sections()
            self._build_conditions_restraints()
            self._build_conditions_diaphragms()
            self._build_conditions_joint_loads()
            self._build_conditions_joint_masses()
            self._build_mesh()
            self._build_analysis_step_recorder()
            self._build_analysis_step_constraint_pattern()
            self._build_analysis_step_load_patterns()
            self._build_analysis('Gravity Analysis', 'Static', load_const=True, wipe=False, duration=1.0, num_incr=10)
            self._build_eigen('Post Gravity')
            done, duration, num_incr = self._build_uniform_excitation()
            if done:
                self._build_monitors()
                self._build_analysis('Dynamic Analysis', 'Transient', load_const=True, wipe=True, duration=duration, num_incr=num_incr)
            self._build_finalization()
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
            E : Dict[int, FxOccShape] = {} # edges from frame members
            EN : Dict[Tuple[int, int], int] = {} # vertex_map to E
            for i, f in geom.frames.items():
                # create the edge from existing vertices in V
                n1, n2 = f.nodes
                edge = FxOccBuilder.makeEdge(V[n1], V[n2])
                E[i] = edge
                EN[tuple(sorted((n1, n2)))] = i
            # build areas
            # map the ETABS area id to the STKO face object
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

    # maps an ETABS _geometry to an STKO MpcGeometry (called by _build_geometries)
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

    # builds the local axes and assign them to frame and area geometries in STKO to match
    # either the default or the user defined local axes in ETABS
    def _build_local_axes(self):
        '''
        Frames:
        - STKO defaults:
            - vertical elements: local_y ~= -global_y (0,-1,0)
            - others: local_z ~= global_z (0,0,1)
        - ETABS defaults:
            - vertical elements: local_y ~= global_x (1,0,0)
            - others: local_y ~= global_z (0,0,1)
        Areas:
        - STKO and ETABS defaults are the same
        '''
        # tolerance for normalized axis components
        tol = 1.0e-3

        # utilities for reversing subgeoms
        def _check_reverse_subshape(stype:MpcSubshapeType, etabs_id:int, axis:Math.vec3):
            # compare to geom
            stko_dz = Math.vec3()
            stko_dx = Math.vec3()
            stko_dy = Math.vec3()
            if stype == MpcSubshapeType.Edge:
                edge_data = self._frame_map.get(etabs_id, None)
                if edge_data is None:
                    raise Exception(f'Frame {etabs_id} not found in geometry')
                edge_data.geom.getLocalAxesOnEdge(edge_data.subshape_id, 0.0, stko_dx, stko_dy, stko_dz)
                if axis.dot(stko_dx) < 0.0:
                    print(f'Edge {etabs_id} has inconsistent orientation with the local axes. Reversing it...')
                    edge_data.geom.shape.toggleEdgeReversedFlag(edge_data.subshape_id)
            elif stype == MpcSubshapeType.Face:
                area_data = self._area_map.get(etabs_id, None)
                if area_data is None:
                    raise Exception(f'Area {etabs_id} not found in geometry')
                area_data.geom.getLocalAxesOnFace(area_data.subshape_id, 0.0, 0.0, stko_dx, stko_dy, stko_dz)
                if axis.dot(stko_dz) < 0.0:
                    print(f'Area {etabs_id} has inconsistent orientation with the local axes. Reversing it...')
                    area_data.geom.shape.toggleFaceReversedFlag(area_data.subshape_id)
            else:
                raise Exception(f'Unsupported subshape type {stype}')

        # next local axes id in STKO
        next_locax_id = self.stko.new_local_axes_id()
        # maps a tuple of (x,y,z) to the local axis id
        # the key is a tuple of integers obtained as the rounded values of the axis components
        locax_id_map : Dict[Tuple[Tuple[int,int,int], Tuple[int,int,int], Tuple[int,int,int]], int] = {}
        def _make_key(x:Math.vec3, y:Math.vec3, z:Math.vec3) -> Tuple[Tuple[int,int,int], Tuple[int,int,int], Tuple[int,int,int]]:
            return (tuple(int(v/tol) for v in x), tuple(int(v/tol) for v in y), tuple(int(v/tol) for v in z))
        # maps the local axis id to the local axis object in STKO
        locax_map : Dict[int, MpcLocalAxes] = {}
        # maps a frame to the local axis id
        frame_to_locax_id_map : Dict[int, int] = {}
        # maps an area to the local axis id
        area_to_locax_id_map : Dict[int, int] = {}

        # define local axes for frames as per ETABS
        for i, f in self.etabs_doc.frames.items():
            p1 = self.etabs_doc.vertices[f.nodes[0]]
            p2 = self.etabs_doc.vertices[f.nodes[1]]
            # first axis always alligned to the frame
            dx = (p2 - p1).normalized()
            # obtain trial dy as per ETABS convention
            dy = Math.vec3(1,0,0) if abs(dx.z) >= 1.0-tol else Math.vec3(0,0,1)
            # compute dz by right hand rule
            dz = dx.cross(dy).normalized()
            # make dy orhtogonal to dx and dz
            dy = dz.cross(dx).normalized()
            # check if the subshape is reversed in STKO
            _check_reverse_subshape(MpcSubshapeType.Edge, i, dx)
            # rotate if user angle is defined
            if abs(f.angle) > 1.0e-10:
                q = Math.quaternion.fromAxisAngle(dx, f.angle/180.0*math.pi)
                dy = q.rotate(dy)
                dz = q.rotate(dz)
                # generate ad hoc local axis reference using ETABS local axis
                GX, GY, GZ = dx, dy, dz
            else:
                # generate local axis reference using default ETABS local axis
                if(abs(dx.z) >= 1.0-tol):
                    # vertical element
                    # define a rectangular local axis object with dy = X+ (opposite for dx = -Z?)
                    GY = Math.vec3(1.0, 0.0, 0.0)
                    GZ = Math.vec3(0.0, 0.0, 1.0)
                    GX = Math.vec3(0.0,-1.0, 0.0)
                else:
                    # GY should point to the global Z axis 
                    GY = Math.vec3(0.0, 0.0, 1.0)
                    if abs(dx.x) < tol:
                        # aligned with global Y axis (enforce local Y = global Y+)
                        if dx.y > 0.0:
                            GZ = Math.vec3(1.0, 0.0, 0.0)
                            GX = Math.vec3(0.0, 1.0, 0.0)
                        else:
                            GZ = Math.vec3(-1.0, 0.0, 0.0)
                            GX = Math.vec3(0.0, -1.0, 0.0)
                    else:
                        # aligned with X or skew in XY plane projection
                        # (enforce local Y = global Y+)
                        if dx.x > 0.0:
                            GZ = Math.vec3(0.0,-1.0, 0.0)
                            GX = Math.vec3(1.0, 0.0, 0.0)
                        else:
                            GZ = Math.vec3(0.0, 1.0, 0.0)
                            GX = Math.vec3(-1.0, 0.0, 0.0)
            # map it
            locax_key = _make_key(GX, GY, GZ)
            locax_id = locax_id_map.get(locax_key, None)
            if locax_id is None:
                locax_id = next_locax_id
                locax_id_map[locax_key] = next_locax_id
                next_locax_id += 1
            frame_to_locax_id_map[i] = locax_id
        
        # define local axes for areas as per ETABS
        for i, a in self.etabs_doc.areas.items():
            # find the normal to the area
            # compute trial local x
            p1 = self.etabs_doc.vertices[a.nodes[0]]
            p2 = self.etabs_doc.vertices[a.nodes[1]]
            dx = (p2 - p1).normalized()
            # try the first point not aligned with dx
            # it may be a general polygon (aligned sides) but assumed as a plane,
            # so there should be at least 3 points that form a plane
            dz = None
            for j in range(2, len(a.nodes)):
                p3 = self.etabs_doc.vertices[a.nodes[j]]
                dy = (p3 - p1).normalized()
                _dz = dx.cross(dy).normalized()
                if _dz.norm() > 0.0:
                    dz = _dz
                    break
            if dz is None:
                raise Exception(f'Cannot find plane in Area {i}')
            # set default
            if abs(dz.z) >= 1.0-tol:
                dy = Math.vec3(0,1,0)
            else:
                dy = Math.vec3(0,0,1)
            # make dx orthogonal to dy and dz
            dx = dy.cross(dz).normalized()
            # check if the subshape is reversed in STKO
            _check_reverse_subshape(MpcSubshapeType.Face, i, dz)
            # rotate if user define angle is defined
            if abs(a.angle) > 1.0e-10:
                q = Math.quaternion.fromAxisAngle(dz, a.angle/180.0*math.pi)
                dx = q.rotate(dx)
                dy = q.rotate(dy)
                # map it (only if angle is defined, because the default should be the same as in STKO)
                locax_key = _make_key(dx, dy, dz)
                locax_id = locax_id_map.get(locax_key, None)
                if locax_id is None:
                    locax_id = next_locax_id
                    locax_id_map[locax_key] = next_locax_id
                    next_locax_id += 1
                area_to_locax_id_map[i] = locax_id
        
        # now generate local axes in STKO
        for (x,y,_), id in locax_id_map.items():
            # create the local axis object
            p1 = self.etabs_doc.bbox.minPoint - (self.etabs_doc.bbox.maxPoint - self.etabs_doc.bbox.minPoint).norm()*0.1*Math.vec3(1,1,1)
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

    # builds uniaxial or nd materials in STKO (TODO: make generic material in STKO professional)
    # from the ETABS material properties
    def _build_elastic_materials(self):
        for name, mat in self.etabs_doc.elastic_materials.items():
            # generate the uniaxial material
            prop = MpcProperty()
            prop.id = self.stko.new_physical_property_id()
            prop.name = f'Elastic Material {name} (1D)'
            meta = self.stko.doc.metaDataPhysicalProperty('materials.uniaxial.Elastic')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('E').quantityScalar.value = mat.E1
            prop.XObject = xobj
            self.stko.add_physical_property(prop)
            prop.commitXObjectChanges()
            P1 = prop
            # generate the nd material
            prop = MpcProperty()
            prop.id = self.stko.new_physical_property_id()
            prop.name = f'Elastic Material {name} (3D)'
            meta = self.stko.doc.metaDataPhysicalProperty('materials.nD.ElasticIsotropic')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('E').quantityScalar.value = mat.E1
            xobj.getAttribute('v').real = mat.U12
            xobj.getAttribute('use_rho').boolean = True
            xobj.getAttribute('rho').quantityScalar.value = mat.rho
            prop.XObject = xobj
            self.stko.add_physical_property(prop)
            prop.commitXObjectChanges()
            P2 = prop
            # add the material to the map
            self._elastic_material_map[name] = (P1, P2)

    # builds the area sections in STKO
    def _build_area_sections(self):

        # keeps track of shearK material for MVLEM
        _shearK_map : Dict[float, MpcProperty] = {}
        def _get_shearK_material(shearK:float) -> MpcProperty:
            if shearK in _shearK_map:
                return _shearK_map[shearK]
            # create a new shearK material
            prop = MpcProperty()
            prop.id = self.stko.new_physical_property_id()
            prop.name = f'MVLEM Shear {shearK}'
            meta = self.stko.doc.metaDataPhysicalProperty('materials.uniaxial.Elastic')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('E').quantityScalar.value = shearK
            prop.XObject = xobj
            self.stko.add_physical_property(prop)
            prop.commitXObjectChanges()
            _shearK_map[shearK] = prop
            return prop

        for name, section in self.etabs_doc.thicknesses.items():
            # obtain the material
            mat = self.etabs_doc.elastic_materials.get(section.material, None)
            # generate the area section
            if section.conversion_info is not None:
                if section.conversion_info['type'] == 'MVLEM':
                    # do nothing, it requires just an Element Property
                    prop = MpcElementProperty()
                    prop.id = self.stko.new_element_property_id()
                    prop.name = '{} Property {}'.format('Wall' if section.is_wall else 'Slab', name)
                    meta = self.stko.doc.metaDataElementProperty('beam_column_elements.MVLEM_3D')
                    xobj = MpcXObject.createInstanceOf(meta)
                    xobj.getAttribute('m').integer = section.conversion_info['nfibers']
                    xobj.getAttribute('Thicknesses').quantityVector.value = section.conversion_info['-thick']
                    xobj.getAttribute('Widths').quantityVector.value = section.conversion_info['-width']
                    xobj.getAttribute('Reinforcing_ratios').quantityVector.value = section.conversion_info['-rho']
                    xobj.getAttribute('Concrete_tags').indexVector = [self._nonlinear_material_map[i].id for i in section.conversion_info['-matConcrete']]
                    xobj.getAttribute('Steel_tags').indexVector = [self._nonlinear_material_map[i].id for i in section.conversion_info['-matSteel']]
                    xobj.getAttribute('Shear_tag').index = _get_shearK_material(section.conversion_info['-ShearK'][0]).id
                    prop.XObject = xobj
                    self.stko.add_element_property(prop)
                    prop.commitXObjectChanges()
                    # add the section to the map
                    self._area_property_map[name] = prop
            else:
                prop = MpcProperty()
                prop.id = self.stko.new_physical_property_id()
                prop.name = '{} Property {}'.format('Wall' if section.is_wall else 'Slab', name)
                meta = self.stko.doc.metaDataPhysicalProperty('sections.ElasticMembranePlateSection')
                xobj = MpcXObject.createInstanceOf(meta)
                xobj.getAttribute('E').quantityScalar.value = mat.E1 * section.Fmod if mat else 0.0
                xobj.getAttribute('nu').real = mat.U12 if mat else 0.0
                xobj.getAttribute('h').quantityScalar.value = section.thickness
                xobj.getAttribute('rho').quantityScalar.value = mat.rho if mat else 0.0
                xobj.getAttribute('Ep_mod').real = section.Mmod / section.Fmod
                prop.XObject = xobj
                self.stko.add_physical_property(prop)
                prop.commitXObjectChanges()
                # add the section to the map
                self._area_property_map[name] = prop
    
    # assignes the area sections to the areas in STKO
    def _assign_area_sections(self):
        asd_shell = None
        for area_mat_name, area_ids in self.etabs_doc.area_sections_assignment.items():
            # get the phsyical property id
            prop = self._area_property_map.get(area_mat_name, None)
            if prop is None:
                raise Exception(f'Area material {area_mat_name} not found in STKO')
            # assign the property to the areas
            for area_id in area_ids:
                # get the geometry and subshape id
                area_data = self._area_map.get(area_id, None)
                if area_data is None:
                    raise Exception(f'Area {area_id} not found in geometry')
                # assign the property to the geometry
                area_data.geom.assign(prop, area_data.subshape_id, MpcSubshapeType.Face)
                # if prop is an MpcProperty, create on the fly a shell element property
                if type(prop) == MpcProperty:
                    # create a shell element property
                    if asd_shell is None:
                        asd_shell = MpcElementProperty()
                        asd_shell.id = self.stko.new_element_property_id()
                        asd_shell.name = f'Shell Element Property {area_mat_name}'
                        meta = self.stko.doc.metaDataElementProperty('shell.ASDShellQ4')
                        xobj = MpcXObject.createInstanceOf(meta)
                        asd_shell.XObject = xobj
                        self.stko.add_element_property(asd_shell)
                        asd_shell.commitXObjectChanges()
                    # assign the element property to the area geometry
                    area_data.geom.assign(asd_shell, area_data.subshape_id, MpcSubshapeType.Face)

    # builds the frame sections in STKO
    def _build_frame_sections(self):
        for name, section in self.etabs_doc.sections.items():
            # obtain the material
            mat = self.etabs_doc.elastic_materials.get(section.material, None)
            # generate the frame section
            prop = MpcProperty()
            prop.id = self.stko.new_physical_property_id()
            prop.name = f'Frame Section {name}'
            meta = self.stko.doc.metaDataPhysicalProperty('sections.Elastic')
            xobj = MpcXObject.createInstanceOf(meta)
            # common properties
            xobj.getAttribute('E').quantityScalar.value = mat.E1 if mat else 0.0
            xobj.getAttribute('G').quantityScalar.value = mat.G12 if mat else 0.0
            xobj.getAttribute('A_modifier').real = max(1.0e-3, section.AMod)
            xobj.getAttribute('Asy_modifier').real = max(1.0e-3, section.AsyMod)
            xobj.getAttribute('Asz_modifier').real = max(1.0e-3, section.AszMod)
            xobj.getAttribute('Iyy_modifier').real = max(1.0e-3, section.IyyMod)
            xobj.getAttribute('Izz_modifier').real = max(1.0e-3, section.IzzMod)
            xobj.getAttribute('J_modifier').real = max(1.0e-3, section.JMod)
            xobj.getAttribute('Y/section_offset').quantityScalar.value = section.Oy
            xobj.getAttribute('Z/section_offset').quantityScalar.value = section.Oz
            xobj.getAttribute('Shear Deformable').boolean = True # make it shear deformable
            # define the section
            if section.shape == section.shape_type.rectangle:
                LY, LZ = section.shape_info[:]
                stko_section = MpcBeamSection(
                    MpcBeamSectionShapeType.Rectangular,
                    'section_Box', 'user', self.etabs_doc.units[1],
                    [LZ, LY],
                )
            else:
                stko_section = MpcBeamSection(
                    MpcBeamSectionShapeType.Custom,
                    'section_Custom', 'user', self.etabs_doc.units[1],
                    [section.A, section.Iyy, section.Izz, section.J, section.Sy, section.Sz],
                )
            # set the section
            xobj.getAttribute('Section').customObject = stko_section
            # set the xobject
            prop.XObject = xobj
            self.stko.add_physical_property(prop)
            prop.commitXObjectChanges()
            # add the section to the map
            self._frame_section_map[name] = prop

    # build the frame nonlinear hinges in STKO
    def _build_frame_nonlinear_hinges(self):
        for name, hinge in self.etabs_doc.frame_nonlinear_hinges.items():
            # generate the property
            prop = MpcProperty()
            prop.id = self.stko.new_physical_property_id()
            prop.name = f'Frame Hinge {name}'
            meta = self.stko.doc.metaDataPhysicalProperty('materials.uniaxial.HystereticSM')
            xobj = MpcXObject.createInstanceOf(meta)
            # common properties
            if hinge.D[0] == 0.0:
                hinge.D[0] = hinge.F[0] / self.etabs_doc.penalty_hinges
            xobj.getAttribute('ep').quantityVector.value = hinge.D
            xobj.getAttribute('sp').quantityVector.value = hinge.F
            # set the xobject
            prop.XObject = xobj
            self.stko.add_physical_property(prop)
            prop.commitXObjectChanges()
            # add the hinge to the map
            self._frame_nonlinear_hinge_map[name] = prop
        # now we need an elastic material for the Vz direction of the hinge
        prop = MpcProperty()
        prop.id = self.stko.new_physical_property_id()
        prop.name = 'Frame Hinge Elastic Vz Material'
        meta = self.stko.doc.metaDataPhysicalProperty('materials.uniaxial.Elastic')
        xobj = MpcXObject.createInstanceOf(meta)
        xobj.getAttribute('E').quantityScalar.value = self.etabs_doc.penalty_hinges
        prop.XObject = xobj
        self.stko.add_physical_property(prop)
        prop.commitXObjectChanges()
        self._frame_nonlinear_stiff_mat = prop

    # assignes the frame sections and hinges to the frames in STKO
    def _assign_frame_sections(self):
        '''
        Note: if a frame has only elastic properties, we can make an elastic beam element 
        and assign the elastic section directly to the frame.
        '''
        # create a map of tuples (int, int), where the first is the frame section and the second the hinge (that can be none)
        combined_map : Dict[Tuple[int, int], MpcProperty] = {}
        # the elastic beam element and the hinged beam (to be created only once)
        beam_elements : List[MpcElementProperty] = [None, None]
        def _get_elastic_beam_elem() -> MpcElementProperty:
            if beam_elements[0] is None:
                # create the elastic beam element
                prop = MpcElementProperty()
                prop.id = self.stko.new_element_property_id()
                prop.name = 'Elastic Beam Element'
                meta = self.stko.doc.metaDataElementProperty('beam_column_elements.elasticBeamColumn')
                xobj = MpcXObject.createInstanceOf(meta)
                if self.etabs_doc.kinematics == 'P-Delta':
                    xobj.getAttribute('transfType').string = 'PDelta'
                elif self.etabs_doc.kinematics == 'Large Displacements':
                    xobj.getAttribute('transfType').string = 'Corotational'
                prop.XObject = xobj
                self.stko.add_element_property(prop)
                prop.commitXObjectChanges()
                beam_elements[0] = prop
            return beam_elements[0]
        def _get_hinged_beam_elem() -> MpcElementProperty:
            if beam_elements[1] is None:
                # create the hinged beam element
                prop = MpcElementProperty()
                prop.id = self.stko.new_element_property_id()
                prop.name = 'Hinged Beam Element'
                meta = self.stko.doc.metaDataElementProperty('special_purpose.BeamWithShearHinge')
                xobj = MpcXObject.createInstanceOf(meta)
                xobj.getAttribute('Beam Element').index = _get_elastic_beam_elem().id
                prop.XObject = xobj
                self.stko.add_element_property(prop)
                prop.commitXObjectChanges()
                beam_elements[1] = prop
            return beam_elements[1]
        # process each frame element
        for frame_id, frame in self.etabs_doc.frames.items():
            # get the elastic section (skip if not found)
            section_name = self.etabs_doc.frame_sections_assignment_inverse.get(frame_id, None)
            if section_name is None:
                continue
            section_prop = self._frame_section_map.get(section_name, None)
            if section_prop is None:
                raise Exception(f'Frame section {section_name} not found in STKO')
            # get the hinge if available
            hinge_name = self.etabs_doc.frame_nonlinear_hinges_assignment_inverse.get(frame_id, None)
            # geometry and subshape id
            frame_data = self._frame_map.get(frame_id, None)
            if frame_data is None:
                raise Exception(f'Frame {frame_id} not found in geometry')
            # if the hinge is none, use the elastic section
            if hinge_name is None:
                # assign the section to the geometry
                frame_data.geom.assign(section_prop, frame_data.subshape_id, MpcSubshapeType.Edge)
                frame_data.geom.assign(_get_elastic_beam_elem(), frame_data.subshape_id, MpcSubshapeType.Edge)
            else:
                # get the hinge property or create a new one
                hinge_prop = self._frame_nonlinear_hinge_map.get(hinge_name, None)
                if hinge_prop is None:
                    raise Exception(f'Frame hinge {hinge_name} not found in STKO')
                combined_key = (section_prop.id, hinge_prop.id)
                combined_prop = combined_map.get(combined_key, None)
                if combined_prop is None:
                    # make the combined physical property
                    prop = MpcProperty()
                    prop.id = self.stko.new_physical_property_id()
                    prop.name = f'Frame Section {section_name} + Hinge {hinge_name}'
                    meta = self.stko.doc.metaDataPhysicalProperty('special_purpose.BeamWithShearHingeProperty')
                    xobj = MpcXObject.createInstanceOf(meta)
                    xobj.getAttribute('Beam Property').index = section_prop.id
                    xobj.getAttribute('Vy Material').index = hinge_prop.id
                    xobj.getAttribute('Vz Material').index = self._frame_nonlinear_stiff_mat.id
                    xobj.getAttribute('K').real = self.etabs_doc.penalty_hinges
                    prop.XObject = xobj
                    self.stko.add_physical_property(prop)
                    prop.commitXObjectChanges()
                    combined_map[combined_key] = prop
                    combined_prop = prop
                # assign the combined property to the geometry
                frame_data.geom.assign(combined_prop, frame_data.subshape_id, MpcSubshapeType.Edge)
                # assign the hinged beam element to the geometry
                frame_data.geom.assign(_get_hinged_beam_elem(), frame_data.subshape_id, MpcSubshapeType.Edge)

    # builds the conditions for restraints in STKO
    def _build_conditions_restraints(self):
        # group restraints by type
        # key = tuple(int * 6), value = list of vertex ids
        restraints : DefaultDict[Tuple[int,int,int,int,int,int], List[int]] = DefaultDict(list)
        for i, r in self.etabs_doc.constraints.items():
            restraints[r].append(i)
        # process each group of restraints
        for rtype, asn_nodes in restraints.items():
            # create a new condition
            condition = MpcCondition()
            condition.id = self.stko.new_condition_id()
            condition.name = f'Restraint {rtype}'
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
    
    # builds the conditions for rigid diaphragms in STKO
    def _build_conditions_diaphragms(self):
        # group restraints by type
        # key = tuple(int * 6), value = list of vertex ids
        for name, _ in self.etabs_doc.diaphragms.items():
            # create a new condition
            condition = MpcCondition()
            condition.id = self.stko.new_condition_id()
            condition.name = f'Diaphragm {name}'
            # define xobject
            meta = self.stko.doc.metaDataCondition('Constraints.mp.rigidDiaphragm')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('perpDirn').integer = 3 # default to X-Y plane?
            condition.XObject = xobj
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

    # builds the conditions for joint loads in STKO
    def _build_conditions_joint_loads(self):
        # group loads by load pattern in the joint_load object
        # key = load pattern name, value = dict of node_id : joint loads
        pattern_loads : DefaultDict[str, Dict[int, nodal_load]] = defaultdict(dict)
        for node_id, load in self.etabs_doc.joint_loads.items():
            # add the joint load to the group
            pattern_loads[load.load_pattern][node_id] = load
        # process each pattern load list and group them by load value
        for pattern_name, node_load in pattern_loads.items():
            load_node_map : DefaultDict[Tuple[float, float, float, float, float, float], List[int]] = DefaultDict(list)
            for node_id, load in node_load.items():
                load_node_map[load.value].append(node_id)
            # process each load value and create a condition for it
            tracked_loads = self._pattern_load_ids[pattern_name]
            for load_value, node_ids in load_node_map.items():
                # split load value into force and moment components
                F = Math.vec3(*load_value[:3])
                M = Math.vec3(*load_value[3:])
                for load_vec, meta_name, label in zip((F, M), ('NodeForce', 'NodeCouple'), ('F', 'M')):
                    # skip null
                    if load_vec.norm() < 1.0e-10:
                        continue
                    # create a new condition
                    condition = MpcCondition()
                    condition.id = self.stko.new_condition_id()
                    condition.name = f'{meta_name} {load_value}'
                    # define xobject
                    meta = self.stko.doc.metaDataCondition(f'Loads.Force.{meta_name}')
                    xobj = MpcXObject.createInstanceOf(meta)
                    xobj.getAttribute(label).quantityVector3.value = load_vec
                    condition.XObject = xobj
                    # add the assigned nodes to the condition
                    for i in node_ids:
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
                    # track the pattern name to force condition id mapping
                    tracked_loads.append(condition.id)

    # builds the conditions for joint masses in STKO
    def _build_conditions_joint_masses(self):
        # process each joint mass list and group them by value 
        # (group nodes with same mass, but before split translational and rotational)
        tmass_node_map : DefaultDict[Tuple[float, float, float], List[int]] = defaultdict(list)
        rmass_node_map : DefaultDict[Tuple[float, float, float], List[int]] = defaultdict(list) 
        for node_id, mass in self.etabs_doc.joint_masses.items():
            if mass.mass_xy > 0.0 or mass.mass_z > 0.0:
                tmass = (mass.mass_xy, mass.mass_xy, mass.mass_z)
                tmass_node_map[tmass].append(node_id)
            if mass.mmi_x > 0.0 or mass.mmi_y > 0.0 or mass.mmi_z > 0.0:
                rmass = (mass.mmi_x, mass.mmi_y, mass.mmi_z)
                rmass_node_map[rmass].append(node_id)
        # process each mass value and create a condition for it
        for mass_mapp, meta_name, label in zip((tmass_node_map, rmass_node_map), ('NodeMass', 'NodeRotationalMass'), ('mass', 'mass')):
            # meta
            meta = self.stko.doc.metaDataCondition(f'Mass.{meta_name}')
            for mass, node_ids in mass_mapp.items():
                # create a new condition
                condition = MpcCondition()
                condition.id = self.stko.new_condition_id()
                condition.name = f'{meta_name} {mass}'
                # define xobject
                xobj = MpcXObject.createInstanceOf(meta)
                xobj.getAttribute(label).quantityVector3.value = Math.vec3(*mass)
                condition.XObject = xobj
                # add the assigned nodes to the condition
                for i in node_ids:
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

        # 2. time histories in etabs as path time series definitions in stko
        for name, th in self.etabs_doc.th_functions.items():
            # create a new definition
            definition = MpcDefinition()
            definition.id = self.stko.new_definition_id()
            definition.name = f'Time History {name}'
            # define xobject
            meta = self.stko.doc.metaDataDefinition('timeSeries.Path')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('list_of_values').quantityVector.value = th.values
            xobj.getAttribute('dt').real = th.dt
            definition.XObject = xobj
            # add the definition to the document
            self.stko.add_definition(definition)
            definition.commitXObjectChanges()
            # map
            self._th_ids[name] = definition.id

    # mesh the model
    def _build_mesh(self):
        # get mesh controls
        mc = self.stko.doc.meshControls
        for _, geom_control in mc.geometryControls.items():
            # make all faces quad/structured
            for face_control in geom_control.faceControls:
                face_control.topology = MpcMeshAlgoTopology.QuadHexa
                face_control.algorithm = MpcMeshAlgo.Structured
        # by default STKO uses 1 division per edge, and we are fine with that
        # mesh the geometry
        App.runCommand('BuildMesh')
        # By default partition it with 1 partition: we want to use OpenSeesMP for MUMPS 
        # but at the moment we need only 1 partition for the Eigen command!
        MpcMeshPartitioner(self.stko.doc.mesh, 1).run()

    # build analysis step : recorder
    def _build_analysis_step_recorder(self):
        # create a new analysis step
        step = MpcAnalysisStep()
        step.id = self.stko.new_analysis_step_id()
        step.name = 'Recorder'
        # define xobject
        meta = self.stko.doc.metaDataAnalysisStep('Recorders.MPCORecorder')
        xobj = MpcXObject.createInstanceOf(meta)
        xobj.getAttribute('name').string = f'{self.etabs_doc.name}.mpco'
        xobj.getAttribute('displacement').boolean = True
        xobj.getAttribute('rotation').boolean = True
        xobj.getAttribute('reactionForce').boolean = True
        xobj.getAttribute('reactionMoment').boolean = True
        xobj.getAttribute('modesOfVibration').boolean = True
        xobj.getAttribute('localForce').boolean = True
        xobj.getAttribute('section.force').boolean = True
        # TODO: add more options
        step.XObject = xobj
        # add the analysis step to the document
        self.stko.add_analysis_step(step)
        step.commitXObjectChanges()

    # build analysis step : constraint pattern
    # here we just create 1 pattern with all sp and mp constraints
    def _build_analysis_step_constraint_pattern(self):
        # create a new analysis step
        step = MpcAnalysisStep()
        step.id = self.stko.new_analysis_step_id()
        step.name = 'Constraint Pattern'
        # define xobject
        meta = self.stko.doc.metaDataAnalysisStep('Patterns.addPattern.constraintPattern')
        xobj = MpcXObject.createInstanceOf(meta)
        xobj.getAttribute('sp').indexVector = self._sp_ids
        xobj.getAttribute('mp').indexVector = self._mp_ids
        step.XObject = xobj
        # add the analysis step to the document
        self.stko.add_analysis_step(step)
        step.commitXObjectChanges()
    
    # build an analysis step : load patterns
    def _build_analysis_step_load_patterns(self):
        # process all patterns in the static load cases
        for lc_name, lc in self.etabs_doc.load_cases_static.items():
            for name, pattern_factor in lc.load_patterns:
                pattern = self.etabs_doc.load_patterns.get(name, None)
                if pattern is None:
                    raise Exception(f'Load pattern {name} not found in ETABS document')
                # create a new analysis step
                step = MpcAnalysisStep()
                step.id = self.stko.new_analysis_step_id()
                step.name = f'Load Pattern {name}'
                # define xobject
                meta = self.stko.doc.metaDataAnalysisStep('Patterns.addPattern.loadPattern')
                xobj = MpcXObject.createInstanceOf(meta)
                xobj.getAttribute('tsTag').index = self._linear_time_series_id
                # add conditions
                # - load
                load_ids = self._pattern_load_ids[name]
                if len(load_ids) > 0:
                    xobj.getAttribute('load').indexVector = load_ids
                # scale
                xobj.getAttribute('-fact').boolean = True
                xobj.getAttribute('cFactor').real = pattern_factor
                # set xobj
                step.XObject = xobj
                # add the analysis step to the document
                self.stko.add_analysis_step(step)
                step.commitXObjectChanges()
    
    # build an analysis step
    def _build_analysis(self, name : str, atype : str, load_const : bool = True, wipe : bool = True, duration : float = 1.0, num_incr : int = 10):
        # create a new analysis step
        step = MpcAnalysisStep()
        step.id = self.stko.new_analysis_step_id()
        step.name = name
        # define xobject
        meta = self.stko.doc.metaDataAnalysisStep('Analyses.AnalysesCommand')
        xobj = MpcXObject.createInstanceOf(meta)
        xobj.getAttribute('analysisType').string = atype
        xobj.getAttribute('constraints').string = 'Auto'
        xobj.getAttribute('numbererType').string = 'Parallel Reverse Cuthill-McKee Numberer'
        xobj.getAttribute('system').string = 'Mumps'
        xobj.getAttribute('algorithm').string = 'Krylov-Newton'
        xobj.getAttribute('-maxDim/KrylovNewton').boolean = True
        xobj.getAttribute('maxDim/KrylovNewton').integer = 40
        xobj.getAttribute('testCommand').string = 'Norm Displacement Increment Test'
        xobj.getAttribute('tol/NormDispIncr').real = unit_system.L(1.0e-5, 'm', self.etabs_doc.units[1]) 
        xobj.getAttribute('iter/NormDispIncr').integer = 20
        xobj.getAttribute('Time Step Type').string = 'Adaptive Time Step'
        xobj.getAttribute('numIncr').integer = num_incr
        xobj.getAttribute('loadConst').boolean = load_const
        if load_const:
            xobj.getAttribute('use pseudoTime').boolean = True
            xobj.getAttribute('pseudoTime').real = 0.0
        xobj.getAttribute('wipeAnalysis').boolean = wipe
        if atype == 'Static':
            xobj.getAttribute('staticIntegrators').string = 'Load Control'
            xobj.getAttribute('duration').real = duration
        elif atype == 'Transient':
            xobj.getAttribute('transientIntegrators').string = 'Newmark Method'
            xobj.getAttribute('gamma').real = 0.5
            xobj.getAttribute('beta').real = 0.25
            xobj.getAttribute('duration/transient').real = duration
        else:
            raise Exception(f'Unknown analysis type {atype}')
        # set xobj
        step.XObject = xobj
        # add the analysis step to the document
        self.stko.add_analysis_step(step)
        step.commitXObjectChanges()

    # build an eigenvalue analysis step
    def _build_eigen(self, name : str):
        # define num_modes based on the mode4_ratio in the load case dynamic
        num_modes = 6
        for _, lc in self.etabs_doc.load_cases_dynamic.items():
            if lc.mode4_ratio > num_modes:
                num_modes = lc.mode4_ratio
        # create a new analysis step
        step = MpcAnalysisStep()
        step.id = self.stko.new_analysis_step_id()
        step.name = f'Eigen Analysis - {name}'
        # define xobject
        meta = self.stko.doc.metaDataAnalysisStep('Analyses.eigen')
        xobj = MpcXObject.createInstanceOf(meta)
        xobj.getAttribute('numEigenvalues').integer = num_modes
        # set xobj
        step.XObject = xobj
        # add the analysis step to the document
        self.stko.add_analysis_step(step)
        step.commitXObjectChanges()
        # do the same with the modalProperties
        step = MpcAnalysisStep()
        step.id = self.stko.new_analysis_step_id()
        step.name = f'Modal Properties - {name}'
        # define xobject
        meta = self.stko.doc.metaDataAnalysisStep('Analyses.modalProperties')
        xobj = MpcXObject.createInstanceOf(meta)
        xobj.getAttribute('-print').boolean = True
        xobj.getAttribute('-file').boolean = True
        xobj.getAttribute('reportFileName').string = f'{step.name}.modal'
        # set xobj
        step.XObject = xobj
        # add the analysis step to the document
        self.stko.add_analysis_step(step)
        step.commitXObjectChanges()

    # build the damping step
    def _build_damping(self, lc):
        # create a new analysis step
        step = MpcAnalysisStep()
        step.id = self.stko.new_analysis_step_id()
        step.name = f'Damping - {lc.name}'
        # define xobject
        # now we implement it as a tcl custom script to mimick the ETABS behavior
        meta = self.stko.doc.metaDataAnalysisStep('Misc_commands.customCommand')
        xobj = MpcXObject.createInstanceOf(meta)
        if lc.pro_by == 'Direct':
            xobj.getAttribute('TCLscript').string = f'''
# Output results
set CM {lc.mass_coeff} ;# Rayleigh damping coefficient for mass proportional damping
set CK {lc.stiff_coeff} ;# Rayleigh damping coefficient for stiffness proportional damping
puts "Direct Rayleigh Damping..."
puts "CM (alpha): $CM"
puts "CK (beta):  $CK"
# define OpenSees command
# rayleigh $alphaM $betaK $betaKinit $betaKcomm
rayleigh $CM 0.0 $CK 0.0
'''
        else:
            xobj.getAttribute('TCLscript').string = f'''
# User inputs
set Mode4Ratio {lc.mode4_ratio} ;# Index of reference mode (1-based)
set ProTimeVal1 {lc.pro_time_val1 if lc.pro_by == 'Frequency Ratio' else 1.0/lc.pro_time_val1} ;# Scale factor for F1
set ProTimeVal2 {lc.pro_time_val2 if lc.pro_by == 'Frequency Ratio' else 1.0/lc.pro_time_val2} ;# Scale factor for F2
set ProDamping1 {lc.pro_damping1} ;# Damping ratio at F1
set ProDamping2 {lc.pro_damping2} ;# Damping ratio at F2

puts "Defining Rayleigh Damping..."
puts "    Mode4Ratio : $Mode4Ratio "
puts "   ProTimeVal1 : $ProTimeVal1"
puts "   ProTimeVal2 : $ProTimeVal2"
puts "   ProDamping1 : $ProDamping1"
puts "   ProDamping2 : $ProDamping2"

# Find first .modal file in current directory
set modalFile ""
foreach f [glob *.modal] {{
    set modalFile $f
    break
}}

if {{$modalFile eq ""}} {{
    puts "No .modal file found in current directory."
    exit
}}
puts "   Using modal file: $modalFile"

# Open .modal file
set in [open $modalFile "r"]
set lines [split [read $in] "\n"]
close $in

# Parse frequencies
set freqs {{}}
set inEigen 0
foreach line $lines {{
    if {{[string match "*EIGENVALUE ANALYSIS*" $line]}} {{
        set inEigen 1
    }} elseif {{$inEigen && [regexp {{^\s*\d+\s+[\d.eE+-]+\s+[\d.eE+-]+\s+([\d.eE+-]+)}} $line match freq]}} {{
        lappend freqs $freq
    }} elseif {{$inEigen && [string match "*TOTAL MASS*" $line]}} {{
        break
    }}
}}

set nFreqs [llength $freqs]
if {{$Mode4Ratio > $nFreqs}} {{
    puts "Error: Mode4Ratio ($Mode4Ratio) exceeds number of modes found ($nFreqs)."
    exit
}}

# Get reference frequency
set FR [lindex $freqs [expr {{$Mode4Ratio - 1}}]]

# Compute F1 and F2
set F1 [expr {{$ProTimeVal1 * $FR}}]
set F2 [expr {{$ProTimeVal2 * $FR}}]

# Convert to angular frequencies
set w1 [expr {{2.0 * 3.141592653589793 * $F1}}]
set w2 [expr {{2.0 * 3.141592653589793 * $F2}}]

# Damping ratios
set xi1 $ProDamping1
set xi2 $ProDamping2

# Compute Rayleigh coefficients
set numerator_alpha [expr {{$w1 * $w2 * ($w2 * $xi1 - $w1 * $xi2)}}]
set denominator [expr {{$w2*$w2 - $w1*$w1}}]
set CM [expr {{2.0 * $numerator_alpha / $denominator}}]

set numerator_beta [expr {{$w2 * $xi2 - $w1 * $xi1}}]
set CK [expr {{2.0 * $numerator_beta / $denominator}}]

# Output results
puts "Reference Frequency (FR): $FR Hz"
puts "F1: $F1 Hz"
puts "F2: $F2 Hz"
puts "CM (alpha): $CM"
puts "CK (beta):  $CK"

# define OpenSees command
# rayleigh $alphaM $betaK $betaKinit $betaKcomm
rayleigh $CM 0.0 $CK 0.0
'''
        # set xobj
        step.XObject = xobj
        # add the analysis step to the document
        self.stko.add_analysis_step(step)
        step.commitXObjectChanges()

    # builds default monitors
    def _build_monitors(self):
        # find a reference node
        # try the top-most node in the rigid diaphragms
        control_geom_info = None
        zmax = -float('inf')
        for name, items in self.etabs_doc.diaphragms.items():
            if len(items) > 0:
                retained_id = items[-1]
                retained_z = self.etabs_doc.vertices[retained_id].z
                if retained_z > zmax:
                    zmax = retained_z
                    # get the geometry and subshape id
                    retained_data = self._vertex_map.get(retained_id, None)
                    if retained_data is None:
                        raise Exception(f'Diagram {name} retained vertex {retained_id} not found in geometry')
                    control_geom_info = retained_data
        # if not found, exit
        if control_geom_info is None:
            return
        # get all information about the reaction nodes at restraints
        reaction_geom_info : List[_geometry_map_item] = []
        for i, _ in self.etabs_doc.constraints.items():
            # get the geometry and subshape id
            vertex_data = self._vertex_map.get(i, None)
            if vertex_data is None:
                raise Exception(f'Restraint {i} not found in geometry')
            # add to the list
            reaction_geom_info.append(vertex_data)
        # create 2 selection sets, one for the control node and one for the reaction nodes
        # control node selection set
        control_sset = MpcSelectionSet()
        control_sset.id = self.stko.new_selection_set_id()
        control_sset.name = 'Control Node'
        sset_item = MpcSelectionSetItem()
        sset_item.wholeGeometry = False
        sset_item.vertices.append(control_geom_info.subshape_id)
        control_sset.geometries[control_geom_info.geom.id] = sset_item
        self.stko.add_selection_set(control_sset)
        # reaction nodes selection set
        reaction_sset = MpcSelectionSet()
        reaction_sset.id = self.stko.new_selection_set_id()
        reaction_sset.name = 'Reaction Nodes'
        for vertex_data in reaction_geom_info:
            if vertex_data.geom.id in reaction_sset.geometries:
                # if the geometry is already in the selection set, just append the vertex
                sset_item = reaction_sset.geometries[vertex_data.geom.id]
            else:
                # if the geometry is not in the selection set, create a new item
                sset_item = MpcSelectionSetItem()
                sset_item.wholeGeometry = False
                reaction_sset.geometries[vertex_data.geom.id] = sset_item
            sset_item.vertices.append(vertex_data.subshape_id)
        self.stko.add_selection_set(reaction_sset)
        # now we can build the monitors
        # 1. Force-Displacement-X
        # 2. Force-Displacement-Y
        # 3. Acceleration-Time-X
        # 4. Acceleration-Time-Y
        for mon_type, mon_name in enumerate(('Force-Displacement', 'Acceleration-Time')):
            for dir_label in ('X', 'Y'):
                # create a new monitor
                monitor = MpcAnalysisStep()
                monitor.id = self.stko.new_analysis_step_id()
                monitor.name = f'{mon_name}-{dir_label}'
                # define xobject
                meta = self.stko.doc.metaDataAnalysisStep('Misc_commands.monitor')
                xobj = MpcXObject.createInstanceOf(meta)
                xobj.getAttribute('Monitor Plot').boolean = True
                xobj.getAttribute('Use Custom Name').boolean = True
                xobj.getAttribute('Custom Name').string = f'{mon_name}-{dir_label}'
                if mon_type == 0:
                    # plot-x (displacement)
                    xobj.getAttribute('Type/X').string = 'Results X Axis Plot'
                    xobj.getAttribute('Result/X').string = 'Displacement'
                    xobj.getAttribute('Component/X').string = dir_label
                    xobj.getAttribute('Selection Set/X').index = control_sset.id
                    xobj.getAttribute('Operation/X').string = 'Average'
                    # plot-y (-reaction force)
                    xobj.getAttribute('Type/Y').string = 'Results Y Axis Plot'
                    xobj.getAttribute('Result/Y').string = 'Reaction Force'
                    xobj.getAttribute('Component/Y').string = dir_label
                    xobj.getAttribute('Selection Set/Y').index = reaction_sset.id
                    xobj.getAttribute('Operation/Y').string = 'Sum'
                    xobj.getAttribute('ScaleFactor/Y').real = -1.0
                    # labels
                    xobj.getAttribute('XLabelAppend').string = f'[{self.etabs_doc.units[1]}]'
                    xobj.getAttribute('YLabelAppend').string = f'[{self.etabs_doc.units[0]}]'
                else:
                    # plot-x (time)
                    xobj.getAttribute('Type/X').string = 'Pseudo Time'
                    # plot-y (acceleration)
                    xobj.getAttribute('Type/Y').string = 'Results Y Axis Plot'
                    xobj.getAttribute('Result/Y').string = 'Acceleration'
                    xobj.getAttribute('Component/Y').string = dir_label
                    xobj.getAttribute('Selection Set/Y').index = control_sset.id
                    xobj.getAttribute('Operation/Y').string = 'Average'
                    # labels
                    xobj.getAttribute('XLabelAppend').string = 's'
                    xobj.getAttribute('YLabelAppend').string = f'[{self.etabs_doc.units[1]}/s^2]'
                # set xobj
                monitor.XObject = xobj
                # add the analysis step to the document
                self.stko.add_analysis_step(monitor)
                monitor.commitXObjectChanges()

    # build uniform excitation.
    # returns a tuple (bool, float, int)
    # where the first is True if the uniform excitation was built,
    # the second is the total duration of the uniform excitation,
    # and the third is the number of increments.
    def _build_uniform_excitation(self) -> Tuple[bool, float, int]:
        '''
        TODO: to be decided.
        # now we ask the user what to choose for the uniform excitation.
        '''
        # quick return if no time histories are defined
        if len(self.etabs_doc.load_cases_dynamic) == 0:
            return False, 0.0, 0
        # if there are more than 1 time histories defined,
        # gather information and let the user choose which one to use.
        index = 0
        if len(self.etabs_doc.load_cases_dynamic) > 1:
            # gather information as a list of string.
            # each one has:
            # name, duration, AMax (one for each comp), ProBy
            # we will use this to let the user choose which one to use.
            info = []
            for _, lc in self.etabs_doc.load_cases_dynamic.items():
                dt = lc.step_size
                duration = dt * lc.num_steps
                lc_info = f'{lc.name} - Step: {dt:.3f}s, Duration: {duration:.3f}s, ProBy: {lc.pro_by}'
                for function, direction, scale in lc.functions:
                    # get the time history
                    th = self.etabs_doc.th_functions.get(function, None)
                    if th is None:
                        raise Exception(f'Time history {function} not found in ETABS document')
                    # get the max value
                    amax = max(abs(v) for v in th.values) * scale
                    # add to info
                    lc_info += f', {direction}: max(A): {amax:.3f}'
                info.append(lc_info)
            # ask the user to choose which one to use
            index = self.stko.select_from_list(
                'Please select the dynamic load case to use for the analysis.',
                info,
                default_index=0,
            )
            if index < 0 or index >= len(self.etabs_doc.load_cases_dynamic):
                # user cancelled
                return False, 0.0, 0
        # get the selected load case
        lc = list(self.etabs_doc.load_cases_dynamic.values())[index]
        # get duration and increments for the result
        dt = lc.step_size
        duration = dt * lc.num_steps
        retval = (True, duration, lc.num_steps)
        # create the uniform excitations (one for each item in lc.functions)
        for function, direction, scale in lc.functions:
            # get the time history
            th = self.etabs_doc.th_functions.get(function, None)
            if th is None:
                raise Exception(f'Time history {function} not found in ETABS document')
            # create a analysis step
            step = MpcAnalysisStep()
            step.id = self.stko.new_analysis_step_id()
            step.name = f'({lc.name} - Uniform Excitation {function} - {direction}'
            # define xobject
            meta = self.stko.doc.metaDataAnalysisStep('Patterns.addPattern.UniformExcitation')
            xobj = MpcXObject.createInstanceOf(meta)
            xobj.getAttribute('tsTag').index = self._th_ids[function]
            xobj.getAttribute('direction').string = _globals.uniform_excitation_dir_map[direction]
            xobj.getAttribute('-fact').boolean = True
            xobj.getAttribute('cFactor').real = scale
            # set xobj
            step.XObject = xobj
            # add the analysis step to the document
            self.stko.add_analysis_step(step)
            step.commitXObjectChanges()
        # build the damping step
        self._build_damping(lc)
        # return the result
        return retval

    # do some extra stuff to finalize the import
    def _build_finalization(self):
        self.stko.assign_custom_colors()
        # TODO: keep dialog open...
        # when the dialog is closed, we can ask the user to save the document