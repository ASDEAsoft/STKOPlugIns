from typing import List, Dict, Tuple, DefaultDict, Any
from collections import defaultdict
import math
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

# The nonlinear material class is used to store the properties of a nonlinear material
class nonlinear_material:
    def __init__(self, name:str, mat_type:str, pos_env:List[Tuple[float,float]], neg_env:List[Tuple[float,float]]):
        self.name = name
        self.mat_type = mat_type
        self.pos_env = pos_env
        self.neg_env = neg_env
    def __str__(self):
        return '{} {} {} {}'.format(self.name, self.mat_type, self.pos_env, self.neg_env)
    def __repr__(self):
        return self.__str__()

# The area section class is used to store the properties of an area material
class area_section:
    def __init__(self, name:str, type:str, material:str, thickness:float, Fmod:float, Mmod:float, is_wall:bool=False):
        self.name = name
        self.type = type
        self.material = material
        self.thickness = thickness
        self.Fmod = Fmod
        self.Mmod = Mmod
        self.is_wall = is_wall
        self.conversion_info : Dict[str, Any] = None # conversion info for the area section, if needed
    def __str__(self):
        return '{} {} {} {} {} {} {}'.format(self.name, self.type, self.material, self.thickness, self.Fmod, self.Mmod, 'Wall' if self.is_wall else 'Slab')
    def __repr__(self):
        return self.__str__()

# The frame section class is used to store the properties of a frame material
class frame_section:
    class shape_type:
        generic = 0
        rectangle = 1
    def __init__(self, name:str, shape:shape_type, material:str,
                 A:float, Iyy:float, Izz:float, J:float, Sy:float, Sz:float,
                 Oy:float, Oz:float,
                 AMod:float, AsyMod:float, AszMod:float, IyyMod:float, IzzMod:float, JMod:float,
                 shape_info:List[float]=None):
        self.name = name
        self.shape = shape
        self.material = material
        self.A = A # cross-sectional area
        self.Iyy = Iyy # moment of inertia about y-axis
        self.Izz = Izz # moment of inertia about z-axis
        self.J = J # torsional constant
        self.Sy = Sy # shear correction factor in the y-axis
        self.Sz = Sz # shear correction factor in the z-axis
        self.Oy = Oy # offset in the y-axis
        self.Oz = Oz # offset in the z-axis
        self.AMod = AMod # modifier for the cross-sectional area
        self.AsyMod = AsyMod # modifier for shear area in the y-axis
        self.AszMod = AszMod # modifier for shear area in the z-axis
        self.IyyMod = IyyMod # modifier for moment of inertia about y-axis
        self.IzzMod = IzzMod # modifier for moment of inertia about z-axis
        self.JMod = JMod # modifier for torsional constant
        self.shape_info = shape_info # additional information about the shape (depending on the shape type)
    def __str__(self):
        return f'{self.name} {self.shape} {self.material} {self.A} {self.Iyy} {self.Izz} {self.J} {self.Sy} {self.Sz} {self.Oy} {self.Oz} {self.AMod} {self.AsyMod} {self.AszMod} {self.IyyMod} {self.IzzMod} {self.JMod}'
    def __repr__(self):
        return self.__str__()

# The frame nonlinear hinge class is used to store the properties of a frame nonlinear hinge
class frame_nonlinear_hinge:
    def __init__(self, name:str, D:List[float], F:List[float]):
        self.name = name
        self.D = D
        self.F = F
        N = len(self.D)
        if N != len(self.F):
            raise ValueError('D and F must have the same number of values')
        if N < 2 or N > 7:
            raise ValueError('D and F must have between 2 and 7 values')
    def __str__(self):
        return '{} {} {}'.format(self.name, self.D, self.F)
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

# the static load case
class load_case_static:
    def __init__(self, name:str):
        self.name = name
        self.load_patterns : List[Tuple[str, float]] = [] # list of tuples (load_pattern_name, multiplier)
    def __str__(self):
        return '{}: {}'.format(self.name, self.load_patterns)
    def __repr__(self):
        return self.__str__()

# The document class is used to store the model data
class document:

    # The document class is used to store the model data
    def __init__(self, name:str='ETABS Model'):
        # name
        self.name : str = name
        # units for length, force, temperature
        self.units : Tuple[str,str,str] = ('m', 'N', 'C')
        # The vertices dictionary is used to store the coordinates of the points in 3D space
        self.vertices : Dict[int, Math.vec3] = {}
        # The frames dictionary is used to store the connectivity of the frame members
        self.frames : Dict[int, frame] = {}
        # The areas dictionary is used to store the connectivity of the slab members
        self.areas : Dict[int, area] = {}
        # The elastic materials dictionary is used to store the properties of the elastic materials
        self.elastic_materials : Dict[str, elastic_material] = {}
        # The nonlinear materials dictionary is used to store the properties of the nonlinear materials
        self.nonlinear_materials : Dict[str, nonlinear_material] = {}
        # The area section dictionary is used to store the properties of the area materials
        self.area_sections : Dict[str, area_section] = {}
        # The area section assignment dictionary (key = area section name, value = list of area ids in ETABS)
        self.area_sections_assignment : DefaultDict[str, List[int]] = defaultdict(list)
        # The inverse of the area section assignment dictionary (key = area id, value = area section name)
        self.area_sections_assignment_inverse : Dict[int, str] = {}
        # The frame section dictionary is used to store the properties of the frame materials
        self.frame_sections : Dict[str, frame_section] = {}
        # The frame section assignment dictionary (key = frame section name, value = list of frame ids in ETABS)
        self.frame_sections_assignment : DefaultDict[str, List[int]] = defaultdict(list)
        # The inverse of the frame section assignment dictionary (key = frame id, value = frame section name)
        self.frame_sections_assignment_inverse : Dict[int, str] = {}
        # The frame nonlinear hinge dictionary is used to store the properties of the frame nonlinear hinges
        self.frame_nonlinear_hinges : Dict[str, frame_nonlinear_hinge] = {}
        # TODO: this works only for: 1 hinge per frame, assumed at reldist = 0.5 and converted to BeamWithShearHinge!
        # The frame nonlinear hinge assignment dictionary (key = frame nonlinear hinge name, value = list of frame ids in ETABS)
        self.frame_nonlinear_hinges_assignment : DefaultDict[str, List[int]] = defaultdict(list)
        # The inverse of the frame nonlinear hinge assignment dictionary (key = frame id, value = frame nonlinear hinge name)
        self.frame_nonlinear_hinges_assignment_inverse : Dict[int, str] = {}
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
        # static load cases (key = load case name, value = load case object)
        self.load_cases_static : Dict[str, load_case_static] = {}
        # computed tolerance
        self.bbox = FxBndBox()
        self.tolerance = 1.0e-6
        # penalty value for the model
        self.penalty_hinges = 1.0e12
        # the kinematic type: Linear or P-Delta (TODO: ask kristijan)
        self.kinematics : str = ''

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
        self._compute_penalty()

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

    # compute the penalty for the model based on the stiffness of the elements
    def _compute_penalty(self):
        # penalty for hinge frame elements
        # to be used in hinge properties with rigid-plastic behavior
        Kmax = 0.0
        for frame_id, frame in self.frames.items():
            section_name = self.frame_sections_assignment_inverse.get(frame_id, None)
            if section_name is None:
                continue
            section = self.frame_sections[section_name]
            material = self.elastic_materials.get(section.material, None)
            if material is None:
                continue
            E = material.E1
            G = material.G12
            A = section.A
            Iyy = section.Iyy
            Izz = section.Izz
            p1 = self.vertices[frame.nodes[0]]
            p2 = self.vertices[frame.nodes[1]]
            L = (p2-p1).norm()
            Kaxial = E*A/L # the axial stiffness that relates the axial force to the axial displacement
            Kv2 = 12*E*Iyy/L**3 # the deflection stiffness that relates the shear force to the deflection
            Kv3 = 12*E*Izz/L**3 # the deflection stiffness that relates the shear force to the deflection
            Km2 = 4*E*Iyy/L # the moment stiffness that relates the moment to the rotation
            Km3 = 4*E*Izz/L # the moment stiffness that relates the moment to the rotation
            Kmax = max(Kmax, Kaxial, Kv2, Kv3, Km2, Km3)
        # penalty for hinge frame elements
        self.penalty_hinges = 10.0 ** (math.ceil(math.log10(Kmax)) + 2)
        if self.penalty_hinges == 0.0:
            self.penalty_hinges = 1.0e12