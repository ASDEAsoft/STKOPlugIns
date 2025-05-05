from PyMpc import *
from ETABSImporter.graph import graph
from ETABSImporter.document import document, frame, area
from ETABSImporter.stko_interface import stko_interface
import itertools
from typing import List, Dict, Tuple, DefaultDict
import math

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
            self._build_local_axes()
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
        
        
