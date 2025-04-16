from typing import List, Dict, Tuple, DefaultDict
from collections import defaultdict
from PyMpc import *

# The frame class is used to store the connectivity of a frame member
class frame:
    def __init__(self, nodes:List[int], angle:float=0.0):
        if len(nodes) != 2:
            raise ValueError('Frame must have exactly 2 nodes')
        self.nodes = nodes
        self.angle = angle
    def __str__(self):
        return '({:8},{:8}), A={:.3f}°'.format(*self.nodes, self.angle)
    def __repr__(self):
        return self.__str__()

# The area class is used to store the connectivity of a slab member
class area:
    def __init__(self, nodes:List[int], angle:float=0.0):
        if len(nodes) != 4:
            raise ValueError('Area must have exactly 4 nodes')
        self.nodes = nodes
        self.angle = angle
    def __str__(self):
        return '({:8},{:8},{:8},{:8}), A = {:.3f}°'.format(*self.nodes, self.angle)
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
        # computed tolerance
        self.tolerance = 1.0e-6

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
        # first renumber vertices in 0-based indexing while mergin them
        self.tolerance = self._tolerance()
        self._merge_vertices()

    # compute tolerance
    def _tolerance(self):
        # compute tolerance based on the average size of the model
        bbox = FxBndBox()
        for v in self.vertices.values():
            bbox.add(v)
        dl = bbox.maxPoint - bbox.minPoint
        avg = (dl.x + dl.y + dl.z) / 3.0
        return 1.0e-6*avg
        
    # merge coincident vertices
    def _merge_vertices(self):
        # convert vertices to tuples of int using tol as tolerance
        # to avoid floating point errors
        def _unique_key(v:Math.vec3, tol:float) -> Tuple[int,int,int]:
            return (int(v.x/tol), int(v.y/tol), int(v.z/tol))
        class _unique_val:
            def __init__(self):
                self.id = -1
                self.vertex : Math.vec3 = None
        # create a dictionary to store unique 0-based indices for merged vertices
        unique_vertices : DefaultDict[Tuple[int,int,int], _unique_val] = defaultdict(_unique_val)
        # maps the old indices to the new unique ones
        old_to_new : Dict[int,int] = {} 
        unique_count = 0
        for i,v in self.vertices.items():
            unique_vertex = unique_vertices[_unique_key(v, self.tolerance)]
            if unique_vertex.vertex is None:
                unique_vertex.vertex = v
                unique_vertex.id = unique_count
                unique_count += 1
            old_to_new[i] = unique_vertex.id
        # create a new dictionary to store the unique vertices
        num_merged = len(self.vertices) - len(unique_vertices)
        # remap the vertices to the new unique ones, in all components
        # vertices
        self.vertices = {uv.id:uv.vertex for uv in unique_vertices.values()}
        # frames
        for i,f in self.frames.items():
            f.nodes = [old_to_new[n] for n in f.nodes]
        # areas
        for i,a in self.areas.items():
            a.nodes = [old_to_new[n] for n in a.nodes]
        # diaphragms
        for k,v in self.diaphragms.items():
            self.diaphragms[k] = [old_to_new[n] for n in v]



    