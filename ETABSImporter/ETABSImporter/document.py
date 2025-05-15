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

# The elastic material class is used to store the properties of an elastic material
class elastic_material:
    def __init__(self, name:str, dens_type:str, sw:float, rho:float, E1:float, G12:float, U12:float, A1:float):
        self.name = name
        self.dens_type = dens_type
        self.sw = sw
        self.rho = rho
        self.E1 = E1
        self.G12 = G12
        self.U12 = U12
        self.A1 = A1
    def __str__(self):
        return '{} {} {} {} {} {} {} {}'.format(self.name, self.dens_type, self.sw, self.rho, self.E1, self.G12, self.U12, self.A1)
    def __repr__(self):
        return self.__str__()

# The area material class is used to store the properties of an area material
class area_material:
    def __init__(self, name:str, type:str, material:str, thickness:float, Fmod:float, Mmod:float, is_wall:bool=False):
        self.name = name
        self.type = type
        self.material = material
        self.thickness = thickness
        self.Fmod = Fmod
        self.Mmod = Mmod
        self.is_wall = is_wall
    def __str__(self):
        return '{} {} {} {} {} {} {}'.format(self.name, self.type, self.material, self.thickness, self.Fmod, self.Mmod, 'Wall' if self.is_wall else 'Slab')
    def __repr__(self):
        return self.__str__()

# a load pattern in etabs (Name,IsAuto,Type,SelfWtMult)
class load_pattern:
    def __init__(self, name:str, is_auto:bool, type:str, self_wt_mult:float):
        self.name = name
        self.is_auto = is_auto
        self.type = type
        self.self_wt_mult = self_wt_mult
    def __str__(self):
        return '{}: {} {} {}'.format(self.name, self.is_auto, self.type, self.self_wt_mult)
    def __repr__(self):
        return self.__str__()

# a joint load in etabs (load_pattern,value)
# where value is a tuple of 6 floats (Fx,Fy,Fz,Mx,My,Mz)
class joint_load:
    def __init__(self, load_pattern:str, value:Tuple[float,float,float,float,float,float]):
        self.load_pattern = load_pattern
        self.value = value
    def __str__(self):
        return '{} {}'.format(self.load_pattern, self.value)
    def __repr__(self):
        return self.__str__()

# a joint concentrated mass in etabs
# NodeID,MassXY,MassZ,MMIX,MMIY,MMIZ
class joint_mass:
    def __init__(self, mass_xy:float, mass_z:float, mmi_x:float, mmi_y:float, mmi_z:float):
        self.mass_xy = mass_xy
        self.mass_z = mass_z
        self.mmi_x = mmi_x
        self.mmi_y = mmi_y
        self.mmi_z = mmi_z
    def __str__(self):
        return '{} {} {} {} {}'.format(self.mass_xy, self.mass_z, self.mmi_x, self.mmi_y, self.mmi_z)
    def __repr__(self):
        return self.__str__()

# time history function in etabs
class th_function:
    def __init__(self, name:str, dt:float, values:List[float]):
        self.name = name
        self.dt = dt
        self.values = values
    def __str__(self):
        return '{} {} {}'.format(self.name, self.dt, self.values)
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
        # The elastic materials dictionary is used to store the properties of the elastic materials
        self.elastic_materials : Dict[str, elastic_material] = {}
        # The area materials dictionary is used to store the properties of the area materials
        self.area_materials : Dict[str, area_material] = {}
        # diaphragm dictionary is used to store the rigid diaphragm members, where the key is the name
        self.diaphragms : Dict[str, List[int]] = {}
        # restraints (key = vertex id, value = list of restraint ids 1 or 0 for 6 DOFs)
        self.restraints : Dict[int, Tuple[int,int,int,int,int,int]] = {}
        # load patterns (key = load pattern name, value = load pattern object)
        self.load_patterns : Dict[str, load_pattern] = {}
        # joint loads (key = vertex id, value = joint load object)
        self.joint_loads : Dict[int, joint_load] = {}
        # joint masses (key = vertex id, value = joint mass object)
        self.joint_masses : Dict[int, joint_mass] = {}
        # time history functions (key = function name, value = function object)
        self.th_functions : Dict[str, th_function] = {}
        # computed tolerance
        self.bbox = FxBndBox()
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
        self._tolerance()
        self._merge_vertices()

    # compute tolerance
    def _tolerance(self):
        # compute tolerance based on the average size of the model
        self.bbox = FxBndBox()
        for v in self.vertices.values():
            self.bbox.add(v)
        dl = self.bbox.maxPoint - self.bbox.minPoint
        avg = (dl.x + dl.y + dl.z) / 3.0
        self.tolerance = 1.0e-6*avg
        
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
        # restraints (node id is in the key)
        _restraints = {i:j for i,j in self.restraints.items()}
        self.restraints = {}
        for k,v in _restraints.items():
            self.restraints[old_to_new[k]] = v
        # joint loads (node id is in the key)
        _joint_loads = {i:j for i,j in self.joint_loads.items()}
        self.joint_loads = {}
        for k,v in _joint_loads.items():
            self.joint_loads[old_to_new[k]] = v
        # joint masses (node id is in the key)
        _joint_masses = {i:j for i,j in self.joint_masses.items()}
        self.joint_masses = {}
        for k,v in _joint_masses.items():
            self.joint_masses[old_to_new[k]] = v



    