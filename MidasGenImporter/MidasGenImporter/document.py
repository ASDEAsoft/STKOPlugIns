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

# the link_data class is used to store the properties of a link element
class link_data:
    def __init__(self, 
                 id:int,
                 rigid_flags:Tuple[bool, bool, bool, bool, bool, bool],
                 stiff:Tuple[float, float, float, float, float, float],
                 shear:bool, dry:float, drz:float):
        self.id = id  # ID of the link element in MGT
        self.rigid_flags = rigid_flags # flags for the 6 DOFs (translation in x, y, z and rotation around x, y, z)
        self.stiff = stiff  # stiffness in the 6 DOFs (translation in x, y, z and rotation around x, y, z)
        self.shear = shear  # shear stiffness flag?
        self.dry = dry  # damping in the y direction?
        self.drz = drz  # damping in the z direction?
    def __str__(self):
        return 'Rigid Flags: {}, Stiffness: {}, Shear: {}, Dry: {}, Drz: {}'.format(
            self.rigid_flags, self.stiff, self.shear, self.dry, self.drz)
    def __repr__(self):
        return self.__str__()

# The frame class is used to store the connectivity of a frame member or for a 2-node elastic link
class frame:
    def __init__(self, mat:int, sec:int, nodes:List[int], angle:float=0.0, link:link_data=None):
        if len(nodes) != 2:
            raise ValueError('Frame must have exactly 2 nodes')
        self.mat = mat
        self.sec = sec
        self.nodes = nodes
        self.angle = angle
        self.link = link  # link data, if this is a 2-node elastic link
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
    def __init__(self, name:str, type:int, shape_info:List[float], offset_y:float, offset_z:float, area:float):
        self.name = name
        self.type = type  # type of the section (SB, L, T)
        self.shape_info = shape_info  # list of floats with the shape information
        self.offset_y = offset_y  # offset in the y direction
        self.offset_z = offset_z  # offset in the z direction
        self.area = area  # area of the section
    def __str__(self):
        return '{} {} {} {} {} {}'.format(self.name, self.type, self.shape_info, self.offset_y, self.offset_z, self.area)
    def __repr__(self):
        return self.__str__()

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

# The section scale factors class is used to store the scale factors for a section
# AREA_SF, ASY_SF, ASZ_SF, IXX_SF, IYY_SF, IZZ_SF, WGT_SF, GROUP, iPart
class section_scale_factors:
    def __init__(self, A:float, Asy:float, Asz:float,
                 Ixx:float, Iyy:float, Izz:float,
                 elements:List[int]=None):
        self.A = A # Area scale factor
        self.Asy = Asy # Asy scale factor
        self.Asz = Asz # Asz scale factor   
        self.Ixx = Ixx # Ixx scale factor
        self.Iyy = Iyy # Iyy scale factor
        self.Izz = Izz # Izz scale factor
        self.elements = elements if elements is not None else []  # list of elements that use this scale factor
    def __str__(self):
        return 'A: {}, Asy: {}, Asz: {}, Ixx: {}, Iyy: {}, Izz: {}'.format(
            self.A, self.Asy, self.Asz, self.Ixx, self.Iyy, self.Izz)
    def __repr__(self):
        return self.__str__()

# The thickness scale factors class is used to store the scale factors for a thickness
class thickness_scale_factors:
    def __init__(self, ip_mod:float, oop_mod:float, elements:List[int]=None):
        self.ip_mod = ip_mod  # in-plane modification factor
        self.oop_mod = oop_mod # out-of-plane modification factor
        self.elements = elements if elements is not None else []  # list of elements that use this scale factor
    def __str__(self):  
        return  'IP Mod: {}, OOP Mod: {}'.format(self.ip_mod, self.oop_mod)
    def __repr__(self):
        return self.__str__()

# The diaphragm class is used to store the properties of a rigid diaphragm
class diaphragm:
    def __init__(self, name:str, Z:float, CX:float, CY:float):
        self.name = name  # name of the diaphragm
        self.Z = Z  # elevation of the diaphragm
        self.CX = CX  # center of mass in the x direction
        self.CY = CY  # center of mass in the y direction
    def __str__(self):
        return 'Diaphragm: {}, Z: {}, CX: {}, CY: {}'.format(self.name, self.Z, self.CX, self.CY)
    def __repr__(self):
        return self.__str__()

# The self_weight_load class is used to store the self weight load direction
# direction is a tuple of 3 floats (Dx, Dy, Dz), that will be multiplied by the gravity and mass
# to obtain the self weight load in the x, y, and z directions
class self_weight_load:
    def __init__(self, direction:Tuple[float, float, float]):
        self.direction = direction
    def __str__(self):
        return 'Self Weight Load: {}'.format(self.direction)
    def __repr__(self):
        return self.__str__()

# nodal load
# where value is a tuple of 6 floats (Fx,Fy,Fz,Mx,My,Mz)
class nodal_load:
    def __init__(self, value:Tuple[float,float,float,float,float,float]):
        self.value = value
    def __str__(self):
        return 'Nodal Load: {}'.format(self.value)
    def __repr__(self):
        return self.__str__()
    def __hash__(self):
        return hash(self.value)
    def __eq__(self, other):
        if not isinstance(other, nodal_load):
            return NotImplemented
        return self.value == other.value
    def copy(self, scale_factor:float, add_factor:float):
        new_value = tuple(v * scale_factor + add_factor for v in self.value)
        return nodal_load(new_value)

# beam load
# where value is a tuple of 3 floats (Fx, Fy, Fz)
# local is a boolean that indicates if the load is in the local coordinate system of the element
class beam_load:
    def __init__(self, value:Tuple[float, float, float], local:bool):
        self.value = value  # value is a tuple of 3 floats (Fx, Fy, Fz)
        self.local = local  # if True, the load is in the local coordinate system of the element
    def __str__(self):
        return 'Beam Load: {}, Local: {}'.format(self.value, self.local)
    def __repr__(self):
        return self.__str__()
    def __hash__(self):
        return hash((self.value, self.local))
    def __eq__(self, other):
        if not isinstance(other, beam_load):
            return NotImplemented
        return self.value == other.value and self.local == other.local
    def copy(self, scale_factor:float, add_factor:float):
        new_value = tuple(v * scale_factor + add_factor for v in self.value)
        return beam_load(new_value, self.local)

# pressure load
# where value is a tuple of 3 floats (Fx, Fy, Fz)
# local is a boolean that indicates if the load is in the local coordinate system of the element
class pressure_load:
    def __init__(self, value:Tuple[float, float, float], local:bool):
        self.value = value  # value is a tuple of 3 floats (Fx, Fy, Fz)
        self.local = local  # if True, the load is in the local coordinate system of
    def __str__(self):
        return 'Pressure Load: {}, Local: {}'.format(self.value, self.local)
    def __repr__(self):
        return self.__str__()
    def __hash__(self):
        return hash((self.value, self.local))
    def __eq__(self, other):
        if not isinstance(other, pressure_load):
            return NotImplemented
        return self.value == other.value and self.local == other.local
    def copy(self, scale_factor:float, add_factor:float):
        new_value = tuple(v * scale_factor + add_factor for v in self.value)
        return pressure_load(new_value, self.local)

# floor load (actually equivalent nodal load)
# where value is a tuple of 6 floats (Fx,Fy,Fz)
class floor_nodal_load:
    def __init__(self, value:Tuple[float,float,float]):
        self.value = value
    def __str__(self):
        return 'Floor Nodal Load: {}'.format(self.value)
    def __repr__(self):
        return self.__str__()
    def __hash__(self):
        return hash(self.value)
    def __eq__(self, other):
        if not isinstance(other, floor_nodal_load):
            return NotImplemented
        return self.value == other.value
    def copy(self, scale_factor:float, add_factor:float):
        new_value = tuple(v * scale_factor + add_factor for v in self.value)
        return floor_nodal_load(new_value)

# The load case class is used to store the properties of a load case
class load_case:
    def __init__(self, name:str):
        self.name = name
        self.self_weight:self_weight_load = None  # self weight
        self.nodal_loads:DefaultDict[nodal_load, List[int]] = defaultdict(list)  # dictionary of nodal loads (key: nodal load, value: list of node IDs)
        self.beam_loads:DefaultDict[beam_load, List[int]] = defaultdict(list)  # dictionary of beam loads (key: beam load, value: list of element IDs)
        self.pressure_loads:DefaultDict[pressure_load, List[int]] = defaultdict(list)  # dictionary of pressure loads (key: pressure load, value: list of element IDs)
        self.floor_loads:DefaultDict[floor_nodal_load, List[int]] = defaultdict(list)  # dictionary of equivalent nodal loads (key: nodal load, value: list of node IDs)
    def __str__(self):
        return 'Load Case: {}, Self Weight: {}, Nodal Loads: {}, Beam Loads: {}, Pressure Loads: {}'.format(
            self.name, self.self_weight, self.nodal_loads, self.beam_loads, self.pressure_loads)
    def __repr__(self):
        return self.__str__()
    def check(self):
        # define a helper function to make a counter from a list
        def _make_counter(source:List[int]) -> DefaultDict[int, int]:
            counter = defaultdict(int)
            for item in source:
                counter[item] += 1
            return counter
        # process nodal loads
        for source_dict, label in zip((self.nodal_loads, self.beam_loads, self.pressure_loads, self.floor_loads), 
                                      ('Nodal Loads', 'Beam Loads', 'Pressure Loads', 'Floor Loads')):
            merged_dict:DefaultDict[nodal_load, List[int]] = defaultdict(list)
            for key, targets in source_dict.items():
                # check if the nodes are unique
                counter = _make_counter(targets)
                if len(counter) == len(targets): continue
                # compute a tolerance for the add_factor
                tol = 1e-6 * max(abs(v) for v in key.value)
                for i_target, count in counter.items():
                    if count > 1:
                        # remove this target from the source targets
                        targets.remove(i_target)
                        # compute a new key for this target scaling it by the count
                        # note that if another key may have the same value as the scaled one!
                        # so we iteratively add an add_factor to the key until it is unique
                        add_factor = 0.0
                        merged_key = None
                        for iter in range(100):
                            merged_key = key.copy(count, add_factor)
                            if merged_key not in merged_dict and merged_key not in source_dict:
                                # add the target to the merged nodal loads
                                merged_dict[merged_key].append(i_target)
                                break
                            add_factor += tol
                        if merged_key is None:
                            raise ValueError('Could not merge for key {} and target {} after 100 iterations'.format(key, i_target))
            # add the merged dict to the source dict
            for key, targets in merged_dict.items():
                if key in source_dict:
                    # this should never happen, but if it does, raise an error
                    raise ValueError('Key {} already exists in the source dictionary'.format(key))
                source_dict[key] = targets
            if len(merged_dict) > 0:
                print('   Handled {} duplicated targets in {}'.format(len(merged_dict), label))

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
        # The section scale factors list is used to store the scale factors for the sections
        self.section_scale_factors : List[section_scale_factors] = []
        # The thickness scale factors list is used to store the scale factors for the thicknesses
        self.thickness_scale_factors : List[thickness_scale_factors] = []
        # constraints (key = vertex id, value = list of restraint ids 1 or 0 for 6 DOFs)
        self.constraints : Dict[int, Tuple[int,int,int,int,int,int]] = {}
        # diaphragm dictionary is used to store the rigid diaphragm members, where the key is the name
        self.diaphragms : Dict[str, diaphragm] = {}
        # vertices not to be used in the diaphragm
        self.diaphragm_released_vertices : List[int] = []
        # load cases
        self.load_cases : Dict[str, load_case] = {}
        # nodal masses in X,Y, and Z directions and the associated nodes
        self.masses : DefaultDict[Tuple[float, float, float], List[int]] = defaultdict(list)
        # computed tolerance
        self.bbox = FxBndBox()
        self.tolerance = 1.0e-6
        # penalty value for the model
        self.penalty_hinges = 1.0e12

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
        # restraints (node id is in the key)
        _restraints = {i:j for i,j in self.constraints.items()}
        self.constraints = {}
        for k,v in _restraints.items():
            self.constraints[old_to_new[k]] = v
        # diaphragms-released vertices
        self.diaphragm_released_vertices = [old_to_new[n] for n in self.diaphragm_released_vertices]
        # laods in load cases
        for lc in self.load_cases.values():
            # nodal loads, node id is in the value
            for nl, nodes in lc.nodal_loads.items():
                new_nodes = [old_to_new[n] for n in nodes]
                lc.nodal_loads[nl] = new_nodes
            # floor nodal loads, node id is in the value
            for fl, nodes in lc.floor_loads.items():
                new_nodes = [old_to_new[n] for n in nodes]
                lc.floor_loads[fl] = new_nodes
        # nodal masses, vertex ids are in the value
        _masses = {k:v for k,v in self.masses.items()}
        self.masses = defaultdict(list)
        for k,v in _masses.items():
            self.masses[k] = [old_to_new[i] for i in v]        
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