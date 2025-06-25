from typing import List, Dict, Tuple, DefaultDict, Any
from collections import defaultdict
import math
from PyMpc import *

# some global infoprmation about the structure
class structure_type:
    def __init__(self, self_weight:bool, grav:float):
        self.self_weight = self_weight
        self.grav = grav
    def __str__(self):
        return 'self_weight={}, grav={}'.format(self.self_weight, self.grav)
    def __repr__(self):
        return self.__str__()

# The frame class is used to store the connectivity of a frame member
class frame:
    def __init__(self, mat:int, sec:int, nodes:List[int], angle:float=0.0):
        if len(nodes) != 2:
            raise ValueError('Frame must have exactly 2 nodes')
        self.mat = mat
        self.sec = sec
        self.nodes = nodes
        self.angle = angle
    def __str__(self):
        return '({:8},{:8}), M = {:8}, S = {:8}, A = {:.3f}°'.format(*self.nodes, self.mat, self.sec, self.angle)
    def __repr__(self):
        return self.__str__()

# The area class is used to store the connectivity of a slab member
class area:
    def __init__(self, mat:int, sec:int, nodes:List[int], angle:float=0.0):
        if len(nodes) != 4 and len(nodes) != 3:
            raise ValueError('Area must have exactly 3 or 4 nodes')
        self.mat = mat
        self.sec = sec
        self.nodes = nodes
        self.angle = angle
    def __str__(self):
        return '({}), A = {:.3f}°'.format(self.nodes, self.angle)
    def __repr__(self):
        return self.__str__()

# The group class is used to store the groups of nodes and elements
class group:
    def __init__(self, name:str, nodes:List[int], elements:List[int]):
        self.name = name
        self.nodes = nodes
        self.elements = elements
    def __str__(self):
        return 'Group: {}, Nodes: {}, Elements: {}'.format(self.name, self.nodes, self.elements)
    def __repr__(self):
        return self.__str__()

# The elastic material class is used to store the properties of an elastic material
class elastic_material:
    def __init__(self, name:str, E:float, poiss:float, rho:float):
        self.name = name
        self.E = E  # Young's modulus
        self.poiss = poiss  # Poisson's ratio
        self.rho = rho  # density
    def __str__(self):
        return '{} {} {} {}'.format(self.name, self.E, self.poiss, self.rho)
    def __repr__(self):
        return self.__str__()

# The section class is used to store the properties of a frame cross-section
class section:
    class shape_type:
        SB = 1
        L = 2
        T = 3
    def __init__(self, name:str, type:int, shape_info:List[float], offset_y:float, offset_z:float):
        self.name = name
        self.type = type  # type of the section (SB, L, T)
        self.shape_info = shape_info  # list of floats with the shape information
        self.offset_y = offset_y  # offset in the y direction
        self.offset_z = offset_z  # offset in the z direction

# The thickness class is used to store the properties of an area material
class thickness:
    def __init__(self, name:str, in_thick:float, out_thick:float, offset:float, offset_type:int):
        self.name = name
        self.in_thick = in_thick  # thickness of the area in the in-plane
        self.out_thick = out_thick  # thickness of the area in the out-of-plane
        self.offset = offset  # offset of the area in the out-of-plane
        self.offset_type = offset_type  # type of the offset
    def __str__(self):
        return '{} {} {} {} {}'.format(self.name, self.in_thick, self.out_thick, self.offset, self.offset_type)
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

# the dynamic load case
class load_case_dynamic:
    def __init__(self, name:str, load_type:str, num_steps:int, step_size:float,
                 pro_by:str, mass_coeff:float, stiff_coeff:float,
                 pro_time_val1:float, pro_damping1:float,
                 pro_time_val2:float, pro_damping2:float,
                 mode4_ratio:int):
        self.name = name
        self.functions : List[Tuple[str, str, float]] = [] # list of tuples (function_name, direction (U1, U2....), multiplier)
        self.load_type = load_type
        self.num_steps = num_steps
        self.step_size = step_size
        self.pro_by = pro_by
        self.mass_coeff = mass_coeff
        self.stiff_coeff = stiff_coeff
        self.pro_time_val1 = pro_time_val1
        self.pro_damping1 = pro_damping1
        self.pro_time_val2 = pro_time_val2
        self.pro_damping2 = pro_damping2
        self.mode4_ratio = mode4_ratio
    def __str__(self):
        return ('Name: {}, Load Type: {}, Functions: {}, '
                'Num Steps: {}, Step Size: {}, Pro By: {}, '
                'Mass Coeff: {}, Stiff Coeff: {}, Pro Time Val1: {}, Pro Damping1: {}, '
                'Pro Time Val2: {}, Pro Damping2: {}, Mode4 Ratio: {}').format(
            self.name, self.load_type, self.functions,
            self.num_steps, self.step_size, self.pro_by,
            self.mass_coeff, self.stiff_coeff, self.pro_time_val1, self.pro_damping1,
            self.pro_time_val2, self.pro_damping2, self.mode4_ratio)
    def __repr__(self):
        return self.__str__()

# The document class is used to store the model data
class document:

    # The document class is used to store the model data
    def __init__(self, name:str='Midas Gen Model'):
        # name
        self.name : str = name
        # units for force, length, temperature
        self.units : Tuple[str,str,str] = ('N', 'M', 'C')
        # the structure type (self_weight, grav)
        self.struct_type : structure_type = structure_type(False, 0.0)
        # The vertices dictionary is used to store the coordinates of the points in 3D space
        self.vertices : Dict[int, Math.vec3] = {}
        # The frames dictionary is used to store the connectivity of the frame members
        self.frames : Dict[int, frame] = {}
        # The areas dictionary is used to store the connectivity of the slab members
        self.areas : Dict[int, area] = {}
        # groups of nodes and elements
        self.groups : Dict[str, group] = {} 
        # The elastic materials dictionary is used to store the properties of the elastic materials
        self.elastic_materials : Dict[str, elastic_material] = {}
        # The frame section dictionary is used to store the properties of the frame materials
        self.sections : Dict[str, section] = {}
        # The thickness dictionary is used to store the properties of the area cross sections
        self.thicknesses : Dict[str, thickness] = {}


        # The area section assignment dictionary (key = area section name, value = list of area ids in ETABS)
        self.area_sections_assignment : DefaultDict[str, List[int]] = defaultdict(list)
        # The inverse of the area section assignment dictionary (key = area id, value = area section name)
        self.area_sections_assignment_inverse : Dict[int, str] = {}
        
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
        # dynamic load cases (key = load case name, value = load case object)
        self.load_cases_dynamic : Dict[str, load_case_dynamic] = {}
        # computed tolerance
        self.bbox = FxBndBox()
        self.tolerance = 1.0e-6
        # penalty value for the model
        self.penalty_hinges = 1.0e12
        # the kinematic type: (None, P-Delta, Large Displacements)
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
        self.tolerance = 1.0e-12*avg
        
    # merge coincident vertices
    def _merge_vertices(self):
        # store the initial number of vertices
        num_vertices = len(self.vertices)
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
        # get the final number of vertices
        num_vertices_final = len(self.vertices)
        print(f'Merged {num_merged} vertices, from {num_vertices} to {num_vertices_final} within tolerance {self.tolerance:.3g}')

    # compute the penalty for the model based on the stiffness of the elements
    def _compute_penalty(self):
        # penalty for hinge frame elements
        # to be used in hinge properties with rigid-plastic behavior
        Kmax = 0.0
        for frame_id, frame in self.frames.items():
            section_name = self.frame_sections_assignment_inverse.get(frame_id, None)
            if section_name is None:
                continue
            section = self.sections[section_name]
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
        if Kmax == 0.0:
            Kmax = 1.0
        self.penalty_hinges = 10.0 ** (math.ceil(math.log10(Kmax)) + 2)
        if self.penalty_hinges == 0.0:
            self.penalty_hinges = 1.0e12