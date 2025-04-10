from typing import List, Dict, Tuple
from collections import defaultdict
from shapely.geometry import Polygon, LineString, MultiLineString
from shapely.ops import split, unary_union, linemerge
import numpy as np
from PyMpc import *

# The frame class is used to store the connectivity of a frame member
class frame:
    def __init__(self, nodes:List[int]):
        if len(nodes) != 2:
            raise ValueError('Frame must have exactly 2 nodes')
        self.nodes = nodes
    def __str__(self):
        return '({:8},{:8})'.format(*self.nodes)
    def __repr__(self):
        return self.__str__()

# The area class is used to store the connectivity of a slab member
class area:
    def __init__(self, nodes:List[int]):
        if len(nodes) != 4:
            raise ValueError('Area must have exactly 4 nodes')
        self.nodes = nodes
    def __str__(self):
        return '({:8},{:8},{:8},{:8})'.format(*self.nodes)
    def __repr__(self):
        return self.__str__()

# The document class is used to store the model data
class document:

    # The document class is used to store the model data
    def __init__(self):
        # The vertices dictionary is used to store the coordinates of the points in 3D space
        self.vertices : Dict[int, Math.vec3] = {}
        # The frames dictionary is used to store the connectivity of the frame members
        self.frames : Dict[int, frame] = {}
        # The areas dictionary is used to store the connectivity of the slab members
        self.areas : Dict[int, area] = {}
        # diaphragm dictionary is used to store the rigid diaphragm members, where the key is the name
        self.diaphragms : Dict[str, List[int]] = {}

    # return the string representation of the document
    def __str__(self):
        from io import StringIO
        f = StringIO()
        f.write('Vertices\n')
        for i,v in self.vertices.items():
            f.write(f"{i:8} : ({v.x:8.3g}, {v.y:8.3g},{v.z:8.3g}'\n")
        f.write('Frames\n')
        for i,v in self.frames.items():
            f.write(f"{i:8} : {v}\n")
        f.write('Areas\n')
        for i,v in self.areas.items():
            f.write(f"{i:8} : {v}\n")
        return f.getvalue()
    
    # return the string representation of the document
    def __repr__(self):
        return self.__str__()
    
    # method to process the parsed document, to be called
    # after the document is created
    def process(self):
        # first renumber vertices in 0-based indexing
        ...





    # The splits area members that intersect with frame memebers
    def _split_area_members(self):
        # tolerance
        tol = self._tolerance()
        # process each area.
        # maps the id of the are to remove, with a dictionary of the new area members
        to_add : Dict[int, Dict[int, area]] = {}
        new_vertex_id = max(self.vertices.keys()) + 1
        new_area_id = max(self.areas.keys()) + 1
        for area_id, the_area in self.areas.items():
            print('Processing area {}'.format(area_id))
            # compute the area info
            T, planar, shape, area_z = self._area_info(the_area)
            # skip non-planar areas
            if not planar:
                continue
            # iterate over the frame members, if they are in the area plane,
            # collect them into a splitter geometry, after that do the split
            test_lines = []
            for _, the_frame in self.frames.items():
                # get the coordinates of the line in local coordinates
                pl1 = T * self.vertices[the_frame.nodes[0]].as_vec3()
                pl2 = T * self.vertices[the_frame.nodes[1]].as_vec3()
                # skip if either of the line points is above/below the area plane
                if abs(pl1[2] - area_z) > tol or abs(pl2[2] - area_z) > tol:
                    continue
                # create a line from the line member
                the_line = FxOccBuilder.makeEdge(FxOccBuilder.makeVertex(pl1), FxOccBuilder.makeVertex(pl2))
                test_lines.append(the_line)
            # merge them into a single geometry
            if len(test_lines) == 0: 
                continue
            splitter = FxOccBuilder.makeCompound(test_lines)
            # try the boolean operation
            merged = FxOccBuilder.merge([shape, splitter], tol)
            nfaces = merged.getNumberOfSubshapes(MpcSubshapeType.Face)
            if nfaces > 1:
                # it has been splitted!
                TT = T.transpose() # local to global
                new_areas = {}
                for subid in range(nfaces):
                    subface = merged.getSubshape(subid, MpcSubshapeType.Face)
                    nvert = subface.getNumberOfSubshapes(MpcSubshapeType.Vertex)
                    if nvert != 4:
                        print('Warning: area member has been split into a face with {} vertices'.format(nvert))
                        continue
                    # convert points in global coords and make the new area member
                    new_area_indices = []
                    for vid in range(nvert):
                        gp = TT * subface.vertexPosition(vid)
                        new_vertex = vertex(gp.x, gp.y, gp.z)
                        self.vertices[new_vertex_id] = new_vertex
                        new_area_indices.append(new_vertex_id)
                        new_vertex_id += 1
                    new_areas[new_area_id] = area(new_area_indices)
                    new_area_id += 1
                # add the new areas to the list of areas to add
                to_add[area_id] = new_areas
        # now, remove old and add new areas
        for area_id, new_areas in to_add.items():
            # remove old area
            del self.areas[area_id]
            # add new areas
            for new_area_id, new_area in new_areas.items():
                self.areas[new_area_id] = new_area

    def _process(self):
        # we want to merge any coincident vertices in the model (using a tolerance)
        # but I walso want to map old vertices' indices to new vertices indices
        tol = self._tolerance()

        def _make_vertex_key(v):
            # create a key for the vertex based on its coordinates and the tolerance
            return (round(v.x/tol), round(v.y/tol), round(v.z/tol))

        class _vertex_value:
            # this class is used to store the vertex value and its index
            def __init__(self):
                self.v = None
                self.i = None
        
        self.vertices = {
            1 : vertex(0,0,0),
            2 : vertex(1+1e-6,0,0),
            3 : vertex(0,1,0),
            4 : vertex(1,0,0),
        }
        self.frames = {}
        self.areas = {}

        unique_vertices = defaultdict(_vertex_value)
        vertex_map : Dict[int, int] = {}
        for i,v in self.vertices.items():
            # create a key for the vertex
            key = _make_vertex_key(v)
            # get existing or create new vertex value
            val = unique_vertices[key]
            # if val.v is none, set it to the vertex
            if val.v is None:
                val.v = v
                val.i = i
            # map old it to new (unique) id
            vertex_map[i] = val.i
        
        print(vertex_map)

    def plot(self, vertices:bool=True, frames:bool=True, areas:bool=True):
        from shapely.geometry import Polygon, LineString
        from shapely.ops import split
        import plotly
        import plotly.graph_objects as go
        
        # plot this document using plotly
        fig = plotly.graph_objects.Figure()

        # settings
        lighting_effects = dict(
            ambient=0.4, diffuse=0.5, roughness = 0.9, specular=0.6, fresnel=0.1, 
            facenormalsepsilon=1.0e-12, vertexnormalsepsilon=1.0e-12)
        light_position = dict(x=2000, y=2000, z=10000)

        # plot the vertices
        if vertices:
            # plot the vertices as points
            x = []
            y = []
            z = []
            for i,v in self.vertices.items():
                x.append(v.x)
                y.append(v.y)
                z.append(v.z)
            fig.add_trace(plotly.graph_objects.Scatter3d(x=x, y=y, z=z, mode='markers', name='Vertices'))

        # plot the frames
        if frames:
            Ex = []
            Ey = []
            Ez = []
            for i,f in self.frames.items():
                v1 = self.vertices[f.nodes[0]]
                v2 = self.vertices[f.nodes[1]]
                Ex.extend([v1.x, v2.x, None])
                Ey.extend([v1.y, v2.y, None])
                Ez.extend([v1.z, v2.z, None])
            fig.add_trace(go.Scatter3d(x=Ex, y=Ey, z=Ez,
                               mode='lines', 
                               line = dict(color = 'rgb(0,0,0)', width=1.0),
                               name='Frames'))

        # plot the areas as triangulation using go.Mesh3d
        if areas:
            x = []
            y = []
            z = []
            i = []
            j = []
            k = []
            nmap : Dict[int, int] = {}
            for _,a in self.areas.items():
                an = []
                for nid in a.nodes:
                    loc = nmap.get(nid, -1)
                    if loc < 0:
                        loc = len(x)
                        nmap[nid] = loc
                    npos = self.vertices[nid]
                    x.append(npos.x)
                    y.append(npos.y)
                    z.append(npos.z)
                    an.append(loc)
                # triangulate the area into two triangles
                t1 = [an[0], an[1], an[2]]
                t2 = [an[0], an[2], an[3]]
                for t in (t1, t2):
                    i.append(t[0])
                    j.append(t[1])
                    k.append(t[2])
            # add the triangulation to the plot
            fig.add_trace(go.Mesh3d(
                x=x, y=y, z=z, i=i, j=j, k=k, 
                color = 'rgb(0,150,255)',
                lighting = lighting_effects,
                lightposition = light_position,
                name='Areas'))
            
        # set the layout
        fig.update_layout(title='ETABS Model', scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'))
        # show the plot
        fig.show()
    
    '''
    Private methods
    '''

    # computes the tolerance based on the size of the bounding box of the model
    def _tolerance(self):
        # compute the bounding box
        min_x = min([v.x for v in self.vertices.values()])
        max_x = max([v.x for v in self.vertices.values()])
        min_y = min([v.y for v in self.vertices.values()])
        max_y = max([v.y for v in self.vertices.values()])
        min_z = min([v.z for v in self.vertices.values()])
        max_z = max([v.z for v in self.vertices.values()])
        # compute the size of the bounding box
        size_x = max_x - min_x
        size_y = max_y - min_y
        size_z = max_z - min_z
        # compute the tolerance based on the size of the bounding box
        tolerance = max(size_x, size_y, size_z) * 1e-6
        return tolerance
    
    # computes area info
    def _area_info(self, a : area) -> Tuple[np.ndarray, bool, FxOccShape]:
        p1 = self.vertices[a.nodes[0]].as_vec3()
        p2 = self.vertices[a.nodes[1]].as_vec3()
        p3 = self.vertices[a.nodes[2]].as_vec3()
        p4 = self.vertices[a.nodes[3]].as_vec3()
        # directors of first sub-triangle
        dx = p2 - p1
        dy = p3 - p1
        dz = dx.cross(dy)
        dy = dz.cross(dx)
        dz.normalize()
        dy.normalize()
        dx.normalize()
        # normal vector of second sub-triangle
        n2 = (p3-p1).cross(p4-p1).normalized()
        # check if planar
        planar = dz.dot(n2) > 0.9999
        T = Math.mat3(dx, dy, dz)
        T.transposeInPlace()
        shape = None
        if planar:
            pl = [T*p for p in (p1,p2,p3,p4)]
            v = [FxOccBuilder.makeVertex(p) for p in pl]
            e = [FxOccBuilder.makeEdge(v[i],v[j]) for i,j in ((0,1),(1,2),(2,3),(3,0))]
            w = FxOccBuilder.makeWire(e)
            shape = FxOccBuilder.makeFace(w)
        # return the transpose of the orientation matrix (global to local) and the planarity
        return (T, planar, shape, pl[0].z)