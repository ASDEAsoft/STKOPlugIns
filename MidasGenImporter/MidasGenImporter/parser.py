from MidasGenImporter.document import *
from MidasGenImporter.stko_interface import stko_interface
from MidasGenImporter.units_utils import unit_system
from MidasGenImporter.section_utils import get_section_offset, test_section
from collections import defaultdict
import os
import math
from PyMpc import *

def _split_line(contents:str, skip_empty:bool=True) -> List[str]:
    """
    Splits a comma-separated string into a list of stripped substrings.

    Args:
        contents (str): The input string to split by commas.
        skip_empty (bool, optional): If True, empty substrings are omitted from the result. Defaults to True.

    Returns:
        List[str]: A list of non-empty (if skip_empty is True) or all (if skip_empty is False) stripped substrings.
    """
    return [pline for pline in (line.strip() for line in contents.split(',')) if (pline or not skip_empty)]

def _mgt_string_to_id_or_range(s: str) -> List[int]:
    """
    Converts a string like '1to3 5' to a list of integers [1, 2, 3, 5].
    It can handle ranges like '1to3' and single numbers like '5'.
    It should also handle ranges with strides like '197to722by175'
    """
    result = []
    for part in s.split():
        part = part.strip()
        if not part:
            continue
        if 'to' in part:
            start, end = part.split('to')
            start = int(start)
            if 'by' in end:
                end, stride = map(int, end.split('by'))
                # generate the range with stride
                result.extend(range(start, end + 1, stride))
            else:
                # generate the range without stride
                end = int(end)
                result.extend(range(start, end + 1))
        else:
            result.append(int(part))
    return result

class _globals:
    _load_commands = ('*SELFWEIGHT', '*CONLOAD', '*BEAMLOAD', '*PRESSURE')
    _xyz_to_index = {
        'X': 0,
        'Y': 1,
        'Z': 2,
    }

class _floor_polygon:
    class _edge:
        """
        A class to represent an edge of a polygon.
        It is used to store the start and end vertices of the edge,
        and compute the length of the edge.
        """
        def __init__(self, v1:int, v2:int, length:float, ratio:float):
            self.v1 = v1
            self.v2 = v2
            self.length = length
            self.ratio = ratio
    """
    A class to represent a polygon for floor loads.
    It is used to store the vertices of the polygon, assumed to be in the XY plane and
    ordered in a counter-clockwise manner.
    The vertices are stored as a list of tuples (ID, x, y).
    This class computes:
    - the total area of the polygon using the shoelace formula.
    - the total perimeter of the polygon.
    - the list of edges as a list of tuples of indices of the vertices.
    - for each edge, it computes the length and it's ratio to the total perimeter (this will be used to distribute the load along the edges).
    """
    def __init__(self, vertices:List[Tuple[int, float, float]]):
        """
        Initializes the polygon with the given vertices.
        The vertices are assumed to be in the XY plane and ordered in a counter-clockwise manner.
        """
        self.vertices:List[Tuple[int, float, float]] = vertices
        self.edges:List[_floor_polygon._edge] = []
        self.area = 0.0
        self.perimeter = 0.0
        self._compute_polygon_properties()
    def _compute_polygon_properties(self):
        """
        Computes the area, perimeter and edges of the polygon.
        The area is computed using the shoelace formula.
        The perimeter is computed as the sum of the lengths of the edges.
        The edges are stored as a list of tuples (v1, v2) where v1 and v2 are the indices of the vertices.
        For each edge, it computes the length and its ratio to the total perimeter.
        """
        n = len(self.vertices)
        if n < 3:
            raise Exception('A polygon must have at least 3 vertices')
        # compute area using the shoelace formula
        self.area = 0.0
        for i in range(n):
            id1, x1, y1 = self.vertices[i]
            id2, x2, y2 = self.vertices[(i + 1) % n]
            self.area += x1 * y2 - x2 * y1
        self.area = abs(self.area) / 2.0
        # compute total perimeter and edge structure
        self.perimeter = 0.0
        for i in range(n):
            v1 = i
            v2 = (i + 1) % n
            id1, x1, y1 = self.vertices[v1]
            id2, x2, y2 = self.vertices[v2]
            length = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            self.edges.append(_floor_polygon._edge(v1, v2, length, 0.0))
            self.perimeter += length
        # compute the ratio of each edge length to the total perimeter
        for edge in self.edges:
            edge.ratio = edge.length / self.perimeter if self.perimeter > 0 else 0.0
    def __str__(self):
        # let's print the edges' ratios
        edges_str = ', '.join(f'({self.vertices[edge.v1][0]}, {self.vertices[edge.v2][0]}): {edge.ratio:.4f}' for edge in self.edges)
        return f'Polygon with {len(self.vertices)} vertices, Area: {self.area:.3g}, Perimeter: {self.perimeter:.3g}, Edges: [{edges_str}]'
    def __repr__(self):
        return self.__str__()

# The parser class is used to parse the ETABS text file and build the document
class parser:
    
    def __init__(self, fname, interface:stko_interface=None):
        # save the filename
        self.fname = fname
        # get the filename without the path and extension
        fname_only = os.path.splitext(os.path.basename(fname))[0]
        # initialize document
        self.doc = document(name = fname_only)
        # save the interface
        self.interface = interface
        # initialize command dictionary
        # this will hold, for each command, a list of lines that belong to that command.
        self.commands = defaultdict(list)
        # the current load case, used to group commands
        # that belong to the same load case
        # this is used for *USE-STLD, <name of the load case>
        current_load_case = None
        # build command dictionary
        with open(fname, 'r') as ifile:
            # first handle the multi-line commands (\ at the end of the line)
            # read all, replace \\n with \n, then split by \n)
            contents = ifile.read().replace('\\\n', '').replace('\\\r\n', '')
            lines = [pline for pline in (line.strip() for line in contents.splitlines()) if (pline and not pline.startswith(';'))]
            value = None
            for line in lines:
                if line.startswith('*'):
                    # remove comments
                    the_command = line.split(';')[0].strip()
                    # handle the special case of *USE-STLD, <name of the load case>
                    if the_command.startswith('*USE-STLD'):
                        # get the load case name
                        words = _split_line(the_command)
                        if len(words) != 2:
                            raise Exception('Invalid *USE-STLD line: {}, expecting 2 values'.format(the_command))
                        current_load_case = words[1]
                        continue
                    # if we have a current load case, and the command is one of the load commands,
                    # we append the load case name to the command
                    if current_load_case is not None and the_command in _globals._load_commands:
                        the_command += f', {current_load_case}'
                    # get (or create) the command list
                    value = self.commands[the_command]
                else:
                    if value is not None:
                        value.append(line)
        # parse
        self._parse()
    
    def _parse(self):
        self._parse_units()
        self._parse_structure_type()
        self._parse_nodes()
        self._parse_elements()
        self._parse_elastic_links()
        self._parse_groups()
        self._parse_elastic_materials()
        self._parse_frame_sections()
        self._parse_thickness()
        self._parse_section_scale_factors()
        self._parse_thickness_scale_factors()
        self._parse_constraints()
        self._parse_diaphragms()
        self._parse_load_cases()
        self._parse_self_weight()
        self._parse_nodal_loads()
        self._parse_beam_loads()
        self._parse_pressure_loads()
        self._parse_floor_loads()
        self._check_load_cases()

        self._parse_load_patterns()
        self._parse_joint_loads()
        self._parse_joint_masses()
        self._parse_time_history_functions()
        self._parse_load_case_static()
        self._parse_load_case_time_history()
    
    # this function parses the units and adds them to the document
    # the units are stored in the document as a dictionary
    def _parse_units(self):
        for item in self.commands['*UNIT']:
            words = _split_line(item)
            if len(words) != 4:
                raise Exception('Invalid units line: {}, expecting 4 values'.format(item))
            self.doc.units = (words[0], words[1], words[3]) # skip the third value (heat)
            if self.interface is not None:
                self.interface.send_message(f'Units: {self.doc.units}', mtype=stko_interface.message_type.INFO)

    # this function parses the structure information and adds it to the document
    def _parse_structure_type(self):
        for item in self.commands['*STRUCTYPE']:
            words = _split_line(item)
            if len(words) != 10:
                raise Exception('Invalid structure type line: {}, expecting 10 values'.format(item))
            self.doc.struct_type = structure_type(words[4] == 'YES', float(words[5]))
            if self.interface is not None:
                self.interface.send_message(f'Structure type: {self.doc.struct_type}', mtype=stko_interface.message_type.INFO)

    # this function parses the nodes and adds them to the document
    # the nodes are stored in the document as vertices
    def _parse_nodes(self):
        for item in self.commands['*NODE']:
            words = _split_line(item)
            if len(words) != 4:
                raise Exception('Invalid node line: {}, expecting 4 values'.format(item))
            id = int(words[0])
            x = float(words[1])
            y = float(words[2])
            z = float(words[3])
            self.doc.vertices[id] = Math.vec3(x,y,z)
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.vertices)} nodes', mtype=stko_interface.message_type.INFO)
     
    # this function parses the elements and adds them to the document
    def _parse_elements(self):
        '''
        ; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB,                     ; Frame  Element
        ; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL, EXVAL2, bLMT ; Comp/Tens Truss
        ; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID , LCAXIS    ; Planar Element
        ; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8     ; Solid  Element
        what's iWID in plate and wall? (WALL ID, not used)
        '''
        all_types = set()
        for item in self.commands['*ELEMENT']:
            words = _split_line(item)
            # index and type + material and property for all elements
            if len(words) < 4:
                raise Exception('Invalid element line: {}, expecting at least 4 values'.format(item))
            id = int(words[0])
            type = words[1]
            mat = int(words[2])
            pro = int(words[3])
            # depending on the type, we have different number of nodes
            if type == 'BEAM':
                if len(words) != 8:
                    raise Exception('Invalid element line: {}, expecting 8 values for BEAM'.format(item))
                nodes = [int(words[i]) for i in range(4, 6)]
                angle = float(words[6])
                self.doc.frames[id] = frame(mat, pro, nodes, angle)
            elif type == 'PLATE' or type == 'WALL':
                if len(words) < 10:
                    raise Exception('Invalid element line: {}, expecting at least 10 values for PLATE'.format(item))
                nodes = [int(words[i]) for i in range(4, 8)]
                nodes = [i for i in nodes if i in self.doc.vertices]  # filter out nodes not in vertices (trianles have an extra node at 0)
                angle = float(words[-1]) # not sure is an angle... TODO: ask nicola
                self.doc.areas[id] = area(mat, pro, nodes, angle)
            else:
                raise Exception('Unknown element type: {}'.format(type))
        print('Parsed {} frames'.format(len(self.doc.frames)))
        print('Parsed {} areas'.format(len(self.doc.areas)))
    
    # this function parses the elastic links and adds them to the document
    def _parse_elastic_links(self):
        '''
        *ELASTICLINK    ; Elastic Link
        GEN:
        ; iNO, iNODE1, iNODE2, LINK, ANGLE, R_SDx, R_SDy, R_SDz, R_SRx, R_SRy, R_SRz, SDx, SDy, SDz, SRx, SRy, SRz, bSHEAR, DRy, DRz, GROUP
        RIGID:
        ; iNO, iNODE1, iNODE2, LINK, ANGLE, bSHEAR, DRy, DRz, GROUP
        '''
        # get the largest frame ID to assign the elastic links
        id = max(self.doc.frames.keys(), default=0) + 1  # start from the next ID
        for item in self.commands['*ELASTICLINK']:
            words = _split_line(item, skip_empty=False)  # we may have empty strings (for empty group)
            if len(words) < 5:
                raise Exception('Invalid elastic link line: {}, expecting at least 5 values'.format(item))
            link_id = int(words[0])
            node1 = int(words[1])
            node2 = int(words[2])
            link_type = words[3]
            angle = float(words[4])
            if link_type == 'GEN':
                # next 6 are rigid flags
                rigid_flags = tuple((words[i] == 'YES') for i in range(5, 11))  # R_SDx, R_SDy, R_SDz, R_SRx, R_SRy, R_SRz
                # next 6 are spring constants
                spring_constants = tuple(float(words[i]) for i in range(11, 17))  # SDx, SDy, SDz, SRx, SRy, SRz
                # next 3 are shear constants
                b_shear = words[17] == 'YES'  # bSHEAR
                dry = float(words[18])  # DRy
                drz = float(words[19])  # DRz
                # create the link_data object
                link = link_data(link_id, rigid_flags, spring_constants, b_shear, dry, drz)
            elif link_type == 'RIGID':
                rigid_flags = (True, True, True, True, True, True)  # all rigid
                spring_constants = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)  # no springs
                b_shear = words[5] == 'YES'  # bSHEAR
                dry = float(words[6])  # DRy
                drz = float(words[7])  # DRz
                # create the link_data object
                link = link_data(link_id, rigid_flags, spring_constants, b_shear, dry, drz)
            else:
                raise Exception('Unknown elastic link type: {}'.format(link_type))
            # create the elastic link object as a frame
            elastic_link = frame(0, 0, [node1, node2], angle, link=link)
            # add it to the document
            self.doc.frames[id] = elastic_link
            # increment the ID for the next elastic link
            id += 1
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.frames)} elastic links', mtype=stko_interface.message_type.INFO)

    # this function parses the groups and adds them to the document
    def _parse_groups(self):
        '''
        *GROUP    ; Group
        ; NAME, NODE_LIST, ELEM_LIST, PLANE_TYPE
        '''
        for item in self.commands['*GROUP']:
            words = _split_line(item, skip_empty=False) # we may have empty strings (for empty node or elem lists)
            if len(words) < 3:
                raise Exception('Invalid group line: {}, expecting at least 3 values'.format(item))
            name = words[0]
            node_list = _mgt_string_to_id_or_range(words[1])
            elem_list = _mgt_string_to_id_or_range(words[2])
            # create the group
            group_obj = group(name, node_list, elem_list)
            # add it to the document
            self.doc.groups[name] = group_obj
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.groups)} groups', mtype=stko_interface.message_type.INFO)    

    # this function parses the elastic material and adds it to the document as elastic_material
    def _parse_elastic_materials(self):
        '''
        *MATERIAL    ; Material
        ; iMAT, TYPE, MNAME, SPHEAT, HEATCO, PLAST, TUNIT, bMASS, DAMPRATIO, [DATA1]           ; STEEL, CONC, ALUMINUM, USER
        ; iMAT, TYPE, MNAME, SPHEAT, HEATCO, PLAST, TUNIT, bMASS, DAMPRATIO, [DATA2], [DATA2]  ; SRC
        ; [DATA1] : 1, STANDARD, CODE/PRODUCT, DB, USEELAST, ELAST
        ; [DATA1] : 2, ELAST, POISN, THERMAL, DEN, MASS
        ; [DATA1] : 3, Ex, Ey, Ez, Tx, Ty, Tz, Sxy, Sxz, Syz, Pxy, Pxz, Pyz, DEN, MASS         ; Orthotropic
        ; [DATA2] : 1, STANDARD, CODE/PRODUCT, DB, USEELAST, ELAST or 2, ELAST, POISN, THERMAL, DEN, MASS
        '''
        '''
        Now we support only the first line (STEEL, CONC, or USER):
        ; iMAT, TYPE, MNAME, SPHEAT, HEATCO, PLAST, TUNIT, bMASS, DAMPRATIO, [DATA1]           ; STEEL, CONC, ALUMINUM, USER
        and either one of:
        ; [DATA1] : 1, STANDARD, CODE/PRODUCT, DB, USEELAST, ELAST
        ; [DATA1] : 2, ELAST, POISN, THERMAL, DEN, MASS
        based on the first item in the DATA1 list.
        '''
        for item in self.commands['*MATERIAL']:
            words = _split_line(item, skip_empty=False)  # we may have empty strings (for empty DATA1)
            if len(words) < 10:
                raise Exception('Invalid material line: {}, expecting at least 10 values'.format(item))
            mat_id = int(words[0])
            mat_type = words[1]
            mat_name = words[2]
            if mat_type not in ('STEEL', 'CONC', 'USER'):
                raise Exception('Unknown material type: {}'.format(mat_type))
            # now we just need E, poisson, and sw from DATA1
            data1 = words[9:]
            data1_type = int(data1[0])
            if data1_type == 1:
                # from standard... we can take E from the last item
                E = float(data1[-1])
                if mat_type == 'STEEL':
                    poiss =  0.3  # default value for steel
                    sw = 78.5
                elif mat_type == 'CONC':
                    poiss = 0.2
                    sw = 25.0  # default value for concrete
                else:
                    raise Exception('Invalid material type: {} for DATA1 type = 1'.format(mat_type))
                sw *= unit_system.F(1.0, 'KN', self.doc.units[0])/unit_system.L(1.0, 'M', self.doc.units[1])**3  # specific weight in kN/m^3
            elif data1_type == 2:
                # from user data
                if len(data1) < 6:
                    raise Exception('Invalid material line: {}, expecting at least 6 values for DATA1 type = 2'.format(item))
                E = float(data1[1])
                poiss = float(data1[2])
                sw = float(data1[4])
            else:
                raise Exception('Unknown DATA1 type: {}'.format(data1_type))
            rho = sw / self.doc.struct_type.grav
            # create the elastic material
            self.doc.elastic_materials[mat_id] = elastic_material(mat_name, E, poiss, rho)
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.elastic_materials)} elastic materials', mtype=stko_interface.message_type.INFO)

    # this function parses the frame sections and adds them to the document
    def _parse_frame_sections(self):
        '''
        *SECTION    ; Section
        ; iSEC, TYPE, SNAME, [OFFSET], bSD, bWE, SHAPE, [DATA1]                   ; 1st line - DB/USER
        ; [DATA1] : DB-NAME or 2, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10
        ; [OFFSET] : OFFSET, iCENT, iREF, iHORZ, HUSER, iVERT, VUSER
        '''
        for item in self.commands['*SECTION']:
            words = _split_line(item, skip_empty=False)
            if len(words) < 3:
                raise Exception('Invalid section line: {}, expecting at least 3 values'.format(item))
            sec_id = int(words[0])
            sec_type = words[1]
            sec_name = words[2]
            if sec_type == 'DBUSER':
                if len(words) < 24:
                    raise Exception('Invalid section line: {}, expecting at least 23 values for DBUSER'.format(item))
                offset_data = words[3:10]
                shape = words[12]
                if shape == 'SB':
                    shape = section.shape_type.SB
                elif shape == 'T':
                    shape = section.shape_type.T
                elif shape == 'L':
                    shape = section.shape_type.L
                else:
                    raise Exception('Unknown shape type: {}'.format(shape))
                data1 = words[13:23]
                if data1[0] != '2':
                    raise Exception('Invalid DATA1 type: {}, expecting 2 for DBUSER'.format(data1[0]))
                shape_info = [float(i) for i in data1[1:]]
                offset_y, offset_z = get_section_offset(shape, shape_info, offset_data)
                self.doc.sections[sec_id] = section(sec_name, shape, shape_info, offset_y, offset_z)
            else:
                raise Exception('Unknown section type: {}'.format(sec_type))
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.sections)} frame sections', mtype=stko_interface.message_type.INFO)

    # this function parses the thickness of the area sections
    def _parse_thickness(self):
        '''
        mgt command:
        *THICKNESS    ; Thickness
        ; iTHK, TYPE, NAME, bSAME, THIK-IN, THIK-OUT, bOFFSET, OFFTYPE, VALUE ; TYPE=VALUE
        '''
        for item in self.commands['*THICKNESS']:
            words = _split_line(item, skip_empty=False)
            if len(words) < 3:
                raise Exception('Invalid thickness line: {}, expecting at least 3 values'.format(item))
            thk_id = int(words[0])
            thk_type = words[1]
            thk_name = words[2]
            if thk_type != 'VALUE':
                raise Exception('Unknown thickness type: {}'.format(thk_type))
            bsame = words[3] == 'YES'
            in_thickness = float(words[4])
            out_thickness = in_thickness if bsame else float(words[5])
            boffset = words[6] == 'YES'
            offset_type = int(words[7]) #TODO: ask nicola what this is
            offset = float(words[8]) if boffset else 0.0
            # create the thickness object
            thick = thickness(thk_name, in_thickness, out_thickness, offset, offset_type)
            # add it to the document
            self.doc.thicknesses[thk_id] = thick
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.thicknesses)} thicknesses', mtype=stko_interface.message_type.INFO)

    # this function parses section scale factors
    def _parse_section_scale_factors(self):
        '''
        *ELEM-SRF    ; Element Stiffness Scale Factor
        ; ELEM_LIST,  AREA_SF, ASY_SF, ASZ_SF, IXX_SF, IYY_SF, IZZ_SF, WGT_SF, GROUP, iPart
        '''
        for item in self.commands['*ELEM-SRF']:
            words = _split_line(item, skip_empty=False)
            if len(words) < 9:
                raise Exception('Invalid section scale factors line: {}, expecting at least 9 values'.format(item))
            elem_list = _mgt_string_to_id_or_range(words[0])
            area_sf = float(words[1])
            asy_sf = float(words[2])
            asz_sf = float(words[3])
            ixx_sf = float(words[4])
            iyy_sf = float(words[5])
            izz_sf = float(words[6])
            # not supported ones
            wgt_sf = float(words[7])  # weight scale factor, not supported
            if wgt_sf != 1.0 and self.interface is not None:
                self.interface.send_message(f'[ELEM-SRF]: Weight scale factor {wgt_sf} is not supported, ignoring it', mtype=stko_interface.message_type.ERROR)
            # create a section scale factor object
            section_sf = section_scale_factors(area_sf, asy_sf, asz_sf, ixx_sf, iyy_sf, izz_sf, elements=elem_list)
            self.doc.section_scale_factors.append(section_sf)
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.section_scale_factors)} section scale factors', mtype=stko_interface.message_type.INFO)

    # this function parses thickness scale factors
    def _parse_thickness_scale_factors(self):
        '''
        *WALL-SSRF    ; Wall Stiffness Scale Factor
        ; ELEM_LIST, SHEAR_FACTOR, BENDING_FACTOR, GROUP, AXIAL, OUTTORSION, OUTSHEAR, OUTBENDING
        '''
        # note: constraint SHEAR_FACTOR = BENDING_FACTOR
        # note: OUTTORSION, OUTSHEAR, OUTBENDING are not supported, must be 1.0
        # TODO: ask nicola if we can use the same scale factors for all walls
        for item in self.commands['*WALL-SSRF']:
            words = _split_line(item, skip_empty=False)
            if len(words) < 8:
                raise Exception('Invalid thickness scale factors line: {}, expecting at least 8 values'.format(item))
            elem_list = _mgt_string_to_id_or_range(words[0])
            shear_factor = float(words[1])
            bending_factor = float(words[2])
            if shear_factor != bending_factor and self.interface is not None:
                self.interface.send_message(f'[WALL-SSRF]: Shear factor {shear_factor} and bending factor {bending_factor} should be equal, ignoring the difference', 
                                            mtype=stko_interface.message_type.WARNING)
            axial = float(words[4]) 
            outtorsion = float(words[5]) 
            outshear = float(words[6]) 
            outbending = float(words[7]) 
            if axial != 1.0 or outtorsion != 1.0 or outshear != 1.0 or outbending != 1.0:
                if self.interface is not None:
                    self.interface.send_message(f'[WALL-SSRF]: Axial {axial}, out-of-plane torsion {outtorsion}, out-of-plane shear {outshear}, and out-of-plane bending {outbending} scale factors are not supported, ignoring them',
                                                mtype=stko_interface.message_type.ERROR)
            # now consider only in-plane scale and out-of-plane scale factors
            in_plane = axial
            out_of_plane = bending_factor
            # create a thickness scale factors object
            self.doc.thickness_scale_factors.append(thickness_scale_factors(in_plane, out_of_plane, elements=elem_list))
        '''
        *PLATE-SSRF    ; Plate Stiffness Scale Factor
        ; ELEM_LIST, Axial(Fxx), Axial(Fyy), Shear(Fxy), Bending(Mxx), Bending(Myy), Torsion(Mxy), Shear(Vxx), Shear(Vyy), GROUP
        '''
        # similar to the wall scale factors, but with different names, however we build the same thickness scale factors object
        for item in self.commands['*PLATE-SSRF']:
            words = _split_line(item, skip_empty=False)
            if len(words) < 10:
                raise Exception('Invalid thickness scale factors line: {}, expecting at least 10 values'.format(item))
            elem_list = _mgt_string_to_id_or_range(words[0])
            axial_fxx = float(words[1])
            axial_fyy = float(words[2])
            shear_fxy = float(words[3])
            bending_mxx = float(words[4])
            bending_myy = float(words[5])
            torsion_mxy = float(words[6])
            shear_vxx = float(words[7])
            shear_vyy = float(words[8])
            # fxx fyy and fxy should be equal to each other
            if (axial_fxx != axial_fyy or axial_fxx != shear_fxy) and self.interface is not None:
                self.interface.send_message(f'[PLATE-SSRF]: Axial Fxx {axial_fxx}, Fyy {axial_fyy} and Fxy {shear_fxy} should be equal to each other',
                                            mtype=stko_interface.message_type.WARNING)
            # mxx, myy, mxy, vxx and vyy should all be equal to each other
            if (bending_mxx != bending_myy or bending_mxx != torsion_mxy or bending_mxx != shear_vxx or bending_mxx != shear_vyy) and self.interface is not None:
                self.interface.send_message(f'[PLATE-SSRF]: Bending Mxx {bending_mxx}, Myy {bending_myy}, Torsion Mxy {torsion_mxy}, Shear Vxx {shear_vxx} and Vyy {shear_vyy} should be equal to each other',
                                            mtype=stko_interface.message_type.WARNING)
            # now consider only in-plane scale and out-of-plane scale factors
            in_plane = axial_fxx  # or axial_fyy or shear_fxy
            out_of_plane = bending_mxx  # or bending_myy or torsion_mxy
            # create a thickness scale factors object
            self.doc.thickness_scale_factors.append(thickness_scale_factors(in_plane, out_of_plane, elements=elem_list))
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.thickness_scale_factors)} thickness scale factors', mtype=stko_interface.message_type.INFO)

    # this function parses the constraints and adds them to the document
    def _parse_constraints(self):
        '''
        *CONSTRAINT    ; Supports
        ; NODE_LIST, CONST(Dx,Dy,Dz,Rx,Ry,Rz), GROUP
        Note: CONST(Dx,Dy,Dz,Rx,Ry,Rz) is a string like 111000
        '''
        for item in self.commands['*CONSTRAINT']:
            words = _split_line(item, skip_empty=False)
            if len(words) < 3:
                raise Exception('Invalid constraint line: {}, expecting at least 3 values'.format(item))
            node_list = _mgt_string_to_id_or_range(words[0])
            const_str = words[1]
            if len(const_str) != 6:
                raise Exception('Invalid constraint string: {}, expecting 6 characters'.format(const_str))
            const = tuple(int(c) for c in const_str)  # convert to tuple of integers (0 or 1)
            # now, for each node in the node list, add the constraint
            for node_id in node_list:
                self.doc.constraints[node_id] = const
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.constraints)} constraints', mtype=stko_interface.message_type.INFO)

    # this function parses the rigid diaphragms and adds them to the document
    def _parse_diaphragms(self):
        '''
        *STORY    ; Story
        ;  NAME=NAME, LEVEL, bFLDIAP, WINDWX, WINDWY, WINDCX, WINDCY, ECCX, ECCY, IECCX, IECCY, TAFX, TAFY, WTAFX, WTAFY, AreaPosNum ; line 1
        we just need, name, level, WTAFX, WTAFY, only if bFLDIAP is YES
        '''
        for item in self.commands['*STORY']:
            words = _split_line(item, skip_empty=False)
            if len(words) < 15:
                raise Exception('Invalid story line: {}, expecting at least 15 values'.format(item))
            name = words[0].split('=')[1]
            level = float(words[1])
            bfldiap = words[2] == 'YES'
            if not bfldiap:
                # if the diaphragm is not rigid, we skip it
                continue
            wtafx = float(words[13])
            wtafy = float(words[14])
            # create a diaphragm object
            diaphragm_obj = diaphragm(name, level, wtafx, wtafy)
            # add it to the document
            self.doc.diaphragms[name] = diaphragm_obj
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.diaphragms)} diaphragms', mtype=stko_interface.message_type.INFO)
        '''
        *RLS-DIAP    ; Diaphragm Disconnect
        ; NODE_LIST
        '''
        # this command is used to define the vertices that are not part of the diaphragm
        for item in self.commands['*RLS-DIAP']:
            # there is only one line with the node list
            nodes = _mgt_string_to_id_or_range(item)
            self.doc.diaphragm_released_vertices.extend(nodes)
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.diaphragm_released_vertices)} diaphragm released vertices', mtype=stko_interface.message_type.INFO)
    
    # this function parses the load cases and adds them to the document
    def _parse_load_cases(self):
        '''
        *STLDCASE    ; Static Load Cases
        ; LCNAME, LCTYPE, DESC
        '''
        for item in self.commands['*STLDCASE']:
            words = _split_line(item, skip_empty=False)
            if len(words) < 3:
                raise Exception('Invalid load case line: {}, expecting at least 3 values'.format(item))
            lc_name = words[0]
            lc_type = words[1]
            lc_desc = words[2] if len(words) > 2 else ''
            # create a load case object
            load_case_obj = load_case(lc_name)
            # add it to the document
            self.doc.load_cases[lc_name] = load_case_obj
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.load_cases)} load cases', mtype=stko_interface.message_type.INFO)

    # this function parses the self weight load and adds it to the document
    def _parse_self_weight(self):
        '''
        *SELFWEIGHT, LCASE_NAME    ; Self Weight
        ; X, Y, Z, GROUP
        '''
        for item, values in self.commands.items():
            if not '*SELFWEIGHT' in item:
                continue
            # get the load case from the command name
            command_words = _split_line(item)
            if len(command_words) < 2:
                raise Exception('Invalid self weight command: {}, expecting at least 2 values'.format(item))
            lc_name = command_words[1]
            lc = self.doc.load_cases.get(lc_name, None)
            if lc is None:
                raise Exception('Load case {} not found for self weight'.format(lc_name))
            # make sure only one line for the self weight is present
            if len(values) > 1:
                raise Exception('Multiple self weight lines found for load case {}: {}'.format(lc_name, values))
            # parse the values
            words = _split_line(values[0])
            if len(words) < 3:
                raise Exception('Invalid self weight line: {}, expecting at least 3 values'.format(values[0]))
            # get the self weight direction
            dx = float(words[0])
            dy = float(words[1])
            dz = float(words[2])
            # create the self weight load object
            lc.self_weight = self_weight_load((dx, dy, dz))
        if self.interface is not None:
            for _, lc in  self.doc.load_cases.items():
                if lc.self_weight is not None:
                    self.interface.send_message(f'Parsed self weight for load case {lc.name}: {lc.self_weight.direction}', mtype=stko_interface.message_type.INFO)

    # this function parses the nodal loads and adds them to the document
    def _parse_nodal_loads(self):
        '''
        *CONLOAD, LCASE_NAME    ; Nodal Loads
        ; NODE_LIST, FX, FY, FZ, MX, MY, MZ, GROUP,STRTYPENAME
        '''
        for item, values in self.commands.items():
            if not '*CONLOAD' in item:
                continue
            # get the load case from the command name
            command_words = _split_line(item)
            if len(command_words) < 2:
                raise Exception('Invalid nodal load command: {}, expecting at least 2 values'.format(item))
            lc_name = command_words[1]
            lc = self.doc.load_cases.get(lc_name, None)
            if lc is None:
                raise Exception('Load case {} not found for nodal loads'.format(lc_name))
            # parse the values
            for value in values:
                words = _split_line(value, skip_empty=False)  # we may have empty strings (for empty group)
                if len(words) < 7:
                    raise Exception('Invalid nodal load line: {}, expecting at least 7 values'.format(value))
                node_list = _mgt_string_to_id_or_range(words[0])
                fx = float(words[1])
                fy = float(words[2])
                fz = float(words[3])
                mx = float(words[4])
                my = float(words[5])
                mz = float(words[6])
                # create the nodal load object
                nodal_load_obj = nodal_load((fx, fy, fz, mx, my, mz))
                # add it to the load case
                lc.nodal_loads[nodal_load_obj].extend(node_list)  # add the node IDs to the list of nodes for this load
        if self.interface is not None:
            for lc_name, lc in self.doc.load_cases.items():
                if len(lc.nodal_loads) > 0:
                    self.interface.send_message(f'Parsed {len(lc.nodal_loads)} nodal loads for load case {lc_name}', mtype=stko_interface.message_type.INFO)

    # this function parses the beam loads and adds them to the document
    def _parse_beam_loads(self):
        '''
        *BEAMLOAD, LCASE_NAME    ; Element Beam Loads
        ; ELEM_LIST, CMD, TYPE, DIR, bPROJ, [ECCEN], [VALUE], GROUP
        where:
        ECCEN = bECCEN, ECCDIR, I-END, J-END, bJ-END
        VALUE = D1, P1, D2, P2, D3, P3, D4, P4
        '''
        for item, values in self.commands.items():
            if not '*BEAMLOAD' in item:
                continue
            # get the load case from the command name
            command_words = _split_line(item)
            if len(command_words) < 2:
                raise Exception('Invalid beam load command: {}, expecting at least 2 values'.format(item))
            lc_name = command_words[1]
            lc = self.doc.load_cases.get(lc_name, None)
            if lc is None:
                raise Exception('Load case {} not found for beam loads'.format(lc_name))
            # parse the values
            for value in values:
                words = _split_line(value, skip_empty=False)
                if len(words) < 18:
                    raise Exception('Invalid beam load line: {}, expecting at least 18 values'.format(value))
                elem_list = _mgt_string_to_id_or_range(words[0])
                cmd = words[1]
                if cmd != 'BEAM':
                    raise Exception('Invalid beam load command: {}, expecting BEAM'.format(cmd))
                load_type = words[2]
                if load_type != 'UNILOAD':
                    raise Exception('Invalid beam load type: {}, expecting UNILOAD'.format(load_type))
                direction = words[3]
                if len(direction) != 2:
                    raise Exception('Invalid beam load direction: {}, expecting 2 characters'.format(direction))
                is_local = direction[0] == 'L'
                component = _globals._xyz_to_index.get(direction[1], None)
                if component is None:
                    raise Exception('Invalid beam load direction component: {}, expecting X, Y or Z'.format(direction[1]))
                # make sure no projection is used
                proj = words[4]
                if proj != 'NO':
                    raise Exception('Beam loads with projection are not supported, found: {}'.format(words[4]))
                # take the first value of eccen, we don't support eccen loads yet, so make sure it's = NO
                eccen = words[5]
                if eccen != 'NO':
                    raise Exception('Beam loads with eccen are not supported, found: {}'.format(eccen))
                # take the second value of VALUE, that's the constant value in case of UNILOAD
                load_component = float(words[11])
                load_value = (load_component if component == 0 else 0.0,
                              load_component if component == 1 else 0.0,
                              load_component if component == 2 else 0.0)
                # create the beam load object
                beam_load_obj = beam_load(load_value, is_local)
                # add the load to the load case
                lc.beam_loads[beam_load_obj].extend(elem_list)  # add the element IDs to the list of elements for this load
        if self.interface is not None:
            for lc_name, lc in self.doc.load_cases.items():
                if len(lc.beam_loads) > 0:
                    self.interface.send_message(f'Parsed {len(lc.beam_loads)} beam loads for load case {lc_name}', mtype=stko_interface.message_type.INFO)

    # this function parses the pressure loads and adds them to the document
    def _parse_pressure_loads(self):
        '''
        *PRESSURE, LCASE_NAME    ; Pressure Loads
        ; ELEM_LIST, CMD, ETYP, LTYP, DIR, VX, VY, VZ, bPROJ, PU, P1, P2, P3, P4, GROUP,PSLTKEY  ; ETYP=PLATE, LTYP=FACE
        Notes:
        - CMD must be PRES
        - ETYP must be PLATE
        - LTYP must be FACE
        - VX, VY, VZ must be 0
        - bPROJ must be NO
        - PU is the pressure value
        - P1, P2, P3, P4 must be 0.0
        '''
        for item, values in self.commands.items():
            if not '*PRESSURE' in item:
                continue
            # get the load case from the command name
            command_words = _split_line(item)
            if len(command_words) < 2:
                raise Exception('Invalid pressure load command: {}, expecting at least 2 words'.format(item))
            lc_name = command_words[1]
            lc = self.doc.load_cases.get(lc_name, None)
            if lc is None:
                raise Exception('Load case {} not found for pressure loads'.format(lc_name))
            # parse the values
            for value in values:
                words = _split_line(value, skip_empty=False)
                if len(words) < 14:
                    raise Exception('Invalid pressure load line: {}, expecting at least 14 values'.format(value))
                elem_list = _mgt_string_to_id_or_range(words[0])
                cmd = words[1]
                if cmd != 'PRES':
                    raise Exception('Invalid pressure load command: {}, expecting PRES'.format(cmd))
                etyp = words[2]
                if etyp != 'PLATE':
                    raise Exception('Invalid pressure load element type: {}, expecting PLATE'.format(etyp))
                ltyp = words[3]
                if ltyp != 'FACE':
                    raise Exception('Invalid pressure load load type: {}, expecting FACE'.format(ltyp))
                direction = words[4]
                if len(direction) != 2:
                    raise Exception('Invalid pressure load direction: {}, expecting 2 characters'.format(direction))
                is_local = direction[0] == 'L'
                component = _globals._xyz_to_index.get(direction[1], None)
                if component is None:
                    raise Exception('Invalid pressure load direction component: {}, expecting X, Y or Z'.format(direction[1]))
                vx = float(words[5])
                vy = float(words[6])
                vz = float(words[7])
                if vx != 0.0 or vy != 0.0 or vz != 0.0:
                    raise Exception('Invalid VX, VY, VZ components: {}, expecting (0.0, 0.0, 0.0)'.format((vx, vy, vz)))
                proj = words[8]
                if proj != 'NO':
                    raise Exception('Pressure loads with projection are not supported, found: {}'.format(proj))
                pu = float(words[9])
                p1 = float(words[10])
                p2 = float(words[11])
                p3 = float(words[12])
                p4 = float(words[13])
                if p1 != 0.0 or p2 != 0.0 or p3 != 0.0 or p4 != 0.0:
                    raise Exception('Invalid P1, P2, P3, P4 values: {}, expecting (0.0, 0.0, 0.0, 0.0)'.format((p1, p2, p3, p4)))
                if pu == 0.0:
                    continue  # no pressure load
                # make the load vector
                load_vector = (
                    pu if is_local and component == 0 else 0.0,  # X component
                    pu if is_local and component == 1 else 0.0,  # Y component
                    pu if is_local and component == 2 else 0.0,  # Z component
                )
                # create the pressure load object
                pressure_load_obj = pressure_load(load_vector, is_local)
                # add the pressure load to the document
                lc.pressure_loads[pressure_load_obj].extend(elem_list)
        if self.interface is not None:
            for lc_name, lc in self.doc.load_cases.items():
                num_pressure_loads = len(lc.pressure_loads)
                if num_pressure_loads > 0:
                    self.interface.send_message(f'Load case {lc_name} has {num_pressure_loads} pressure loads', mtype=stko_interface.message_type.INFO)

    # this function parses the floor loads and adds them to the document as equivalent loads
    def _parse_floor_loads(self):
        # Step 1: parse the floor load type.
        # build a dictionary of floor load types
        # with: key = floor load type name, value = List[load_case_name:str, load_value:float]
        '''
        *FLOADTYPE    ; Define Floor Load Type
        ; NAME, DESC                                             ; 1st line
        ; LCNAME1, FLOAD1, bSBU1, ..., LCNAME_N, FLOAD_N, bSBU_N  ; 2nd line
        '''
        floor_load_types : Dict[str, List[Tuple[str, float]]] = {}
        lines = self.commands['*FLOADTYPE']
        if len(lines) % 2 != 0:
            raise Exception('Invalid floor load type command, expecting an even number of lines')
        for i in range(0, len(lines), 2):
            line_1 = _split_line(lines[i], skip_empty=False)
            if len(line_1) < 2:
                raise Exception('Invalid floor load type command, expecting at least 2 words in the first line')
            key = line_1[0]
            line_2 = _split_line(lines[i + 1], skip_empty=False)
            # multiple of 3
            if len(line_2) % 3 != 0:
                raise Exception('Invalid floor load type command, expecting a multiple of 3 words in the second line')
            # parse the second line
            value = []
            for j in range(0, len(line_2) - 1, 3):
                lc_name = line_2[j]
                fload = float(line_2[j + 1])
                value.append((lc_name, fload))
            floor_load_types[key] = value
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(floor_load_types)} floor load types', mtype=stko_interface.message_type.INFO)
        # Step 2: parse the floor loads and add them to the document
        '''
        *FLOORLOAD    ; Floor Loads
        ; LTNAME, iDIST, ANGLE, iSBEAM, SBANG, SBUW, DIR, bPROJ, DESC, bEX, bAL, GROUP, NODE1, ..., NODEn
        ; Note: we only suppport iDIST 2=Two Way
        '''
        for item in self.commands['*FLOORLOAD']:
            words = _split_line(item, skip_empty=False)
            if len(words) < 15:
                raise Exception('Invalid floor load command, expecting at least 15 words')
            lt_name = words[0]  # floor load type name
            idist = int(words[1])  # distance type, we only support 2=Two Way
            if idist != 2:
                raise Exception(f'Invalid distance type {idist}, only 2=Two Way is supported')
            # ANGLE, iSBEAM, SBANG, SBUW must be zero
            angle = float(words[2])
            isbeam = int(words[3])
            sbang = float(words[4])
            sbuw = float(words[5])
            if angle != 0.0 or isbeam != 0 or sbang != 0.0 or sbuw != 0.0:
                raise Exception(f'Invalid floor load command, expecting angle, isbeam, sbang and sbuw to be zero, got {angle}, {isbeam}, {sbang}, {sbuw}')
            direction = words[6] 
            if len(direction) != 2:
                raise Exception(f'Invalid floor load command, expecting direction to be a 2-character string, got {direction}')
            is_local = direction[0] == 'L'
            if is_local:
                raise Exception('Local direction is not supported for floor loads')
            component = _globals._xyz_to_index.get(direction[1], None)
            if component is None:
                raise Exception('Invalid floor load direction component: {}, expecting X, Y or Z'.format(direction[1]))
            bproj = words[7]
            if bproj != 'NO':
                raise Exception('Invalid floor load command, expecting bproj to be NO, got {}'.format(bproj))
            # skip DESC, bEX, bAL, GROUP
            # get all nodes (comma separated)
            node_ids = [int(i) for i in words[12:] if i.isdigit()]
            # collect list of tuples(id,x,y) for the polygon vertices
            # and check they are all at the same Z coordinate within a tolerance
            points: List[Tuple[int, float, float]] = []
            z_coord = None
            z_error = 0.0
            bbox = FxBndBox()
            for node_id in node_ids:
                vertex = self.doc.vertices.get(node_id, None)
                if vertex is None:
                    raise Exception(f'Node {node_id} not found in the document')
                bbox.add(vertex)
                # check if the node is at the same Z coordinate
                if z_coord is None:
                    z_coord = vertex.z
                z_error = max(z_error, abs(vertex.z - z_coord))
                points.append((node_id, vertex.x, vertex.y))
            bbox_diag = bbox.maxPoint - bbox.minPoint
            bbox_diag_avg = (bbox_diag.x + bbox_diag.y + bbox_diag.z) / 3.0
            tolerance = 1.0e-12*bbox_diag_avg
            if z_error > tolerance:
                raise Exception(f'Nodes {node_ids} are not at the same Z coordinate, max error = {z_error}, tolerance = {tolerance}')
            # now we can create the polygon
            polygon = _floor_polygon(points)
            # now, we can create the equivalent beam loads
            # for each edge of the polygon
            # and for each item in the associated floor load type
            lt_obj = floor_load_types.get(lt_name, None)
            if lt_obj is None:
                raise Exception(f'Floor load type {lt_name} not found in the document')
            for lc_name, area_load in lt_obj:
                # get the load case from the document
                lc = self.doc.load_cases.get(lc_name, None)
                if lc is None:
                    raise Exception(f'Load case {lc_name} not found in the document')
                # from area load (F/A) to concentrated load on the whole floor (F)
                total_load = area_load * polygon.area
                # now we have to 1) compute the distributed load (F/L) for each edge
                # and 2) lump it to the 2 end-nodes of the edge
                # so we can create equivalent nodal loads
                eq_nodal_loads = defaultdict(float)
                for edge in polygon.edges:
                    # get the portion of the total load to apply on the 2 nodes this edge (F)
                    edge_nodal_load = total_load * edge.ratio / 2.0
                    # accumulate the load on the 2 end-nodes of the edge
                    N1 = polygon.vertices[edge.v1][0]
                    N2 = polygon.vertices[edge.v2][0]
                    eq_nodal_loads[N1] += edge_nodal_load
                    eq_nodal_loads[N2] += edge_nodal_load
                # now we can create the equivalent nodal loads
                for node_id, edge_nodal_load in eq_nodal_loads.items():
                    # create a nodal load object
                    load_value = (edge_nodal_load if component == 0 else 0.0,
                                  edge_nodal_load if component == 1 else 0.0,
                                  edge_nodal_load if component == 2 else 0.0)
                    nodal_load_obj = floor_nodal_load(load_value)
                    # add it to the load case
                    lc.floor_loads[nodal_load_obj].append(node_id)
        if self.interface is not None:
            for lc_name, lc in self.doc.load_cases.items():
                # send a message with the number of equivalent nodal loads
                num_eq_nodal_loads = sum(len(nodes) for _, nodes in lc.floor_loads.items())
                self.interface.send_message(f'Load case {lc_name} has {num_eq_nodal_loads} equivalent nodal loads from floor loads', mtype=stko_interface.message_type.INFO)
            
    # this function checks if the load cases are valid
    def _check_load_cases(self):
        for lc_name, lc in self.doc.load_cases.items():
            if self.interface is not None:
                self.interface.send_message(f'Checking load case {lc_name}', mtype=stko_interface.message_type.INFO)
            lc.check()






    # this function parses the load patterns and adds them to the document
    def _parse_load_patterns(self):
        for item in self.commands['* LOAD_PATTERNS']:
            words = item.split(',')
            name = words[0]
            is_auto = words[1] == 'Yes'
            type = words[2]
            self_wt_mult = float(words[3])
            self.doc.load_patterns[name] = load_pattern(name, is_auto, type, self_wt_mult)
    
    # this function parses the joint loads and adds them to the document
    def _parse_joint_loads(self):
        for item in self.commands['* JOINT_LOADS']:
            words = item.split(',')
            id = int(words[0])
            load_pattern = words[1]
            if load_pattern not in self.doc.load_patterns:
                raise Exception('Load pattern {} not found'.format(load_pattern))
            value = tuple([float(words[i]) for i in range(2, 8)])
            self.doc.joint_loads[id] = nodal_load(load_pattern, value)
    
    # this function parses the joint masses and adds them to the document
    def _parse_joint_masses(self):
        for item in self.commands['* JOINT_ADDITIONAL_MASS']:
            words = item.split(',')
            id = int(words[0])
            mass_xy = float(words[1])
            mass_z = float(words[2])
            mmi_x = float(words[3])
            mmi_y = float(words[4])
            mmi_z = float(words[5])
            self.doc.joint_masses[id] = joint_mass(mass_xy, mass_z, mmi_x, mmi_y, mmi_z)
    
    # this function parses the time history functions and adds them to the document
    def _parse_time_history_functions(self):
        for item in self.commands['* TH_FUNCTIONS']:
            words = item.split(',')
            name = words[0]
            dt = float(words[1])
            file_name = words[2]
            # the input file name might be absolute on the user's machine
            # if we can find the file using the absolute path, we will use it
            # otherwise we will use the relative path
            # to the current self.fname
            if os.path.exists(file_name):
                # it exists on this machine, either relative or absolute
                file_name = os.path.abspath(file_name)
                print('Time history function file {} found'.format(file_name))
            else:
                # it doesn't exist, let's find the file name
                file_name_only = os.path.basename(file_name)
                print('Time history function file {} not found, searching for it in the same folder as the input file...'.format(file_name_only))
                # let's search for the file in the same folder as the input file or in the sub-directories
                # let's get the directory of the input file
                input_file_dir = os.path.dirname(self.fname)
                # let's search for the file in the same folder as the input file or in the sub-directories
                for root, dirs, files in os.walk(input_file_dir):
                    if file_name_only in files:
                        file_name = os.path.join(root, file_name_only)
                        print('   Time history function file {} found'.format(file_name))
                        break
                else:
                    # we didn't find the file
                    raise Exception('Time history function file {} not found'.format(file_name_only))
            # read the file
            with open(file_name, 'r') as ifile:
                values = [float(pline) for pline in (line.strip() for line in ifile.readlines()) if pline]
            # save it
            self.doc.th_functions[name] = th_function(name, dt, values)
    
    # parses static load cases
    def _parse_load_case_static(self):
        '''
        command format:
        * LOAD_CASES_LINEAR_STATIC
        # Name,Group,MassSource,StiffType,LoadType,LoadName,LoadSF,DesignType
        '''
        for item in self.commands['* LOAD_CASES_LINEAR_STATIC']:
            words = item.split(',')
            name = words[0]
            group = words[1]
            mass_source = words[2]
            stiff_type = words[3]
            load_type = words[4]
            load_name = words[5]
            load_sf = float(words[6])
            design_type = words[7] if len(words) > 7 else None
            # add the load case to the document or get existing one
            lc = self.doc.load_cases_static.get(name, None)
            if lc is None:
                # create a new load case
                lc = load_case_static(name)
                self.doc.load_cases_static[name] = lc
            # add the load pattern to the load case
            lc.load_patterns.append((load_name, load_sf))
            # check kinematics
            if self.doc.kinematics: # if already defined...
                if stiff_type != self.doc.kinematics:
                    self.interface.send_message(f'[LOAD_CASES_LINEAR_STATIC {name}]: kinematics type {stiff_type} does not match the existing kinematics type {self.doc.kinematics}', 
                                                mtype=stko_interface.message_type.WARNING)
            else:
                # set the kinematics type
                self.doc.kinematics = stiff_type
    
    # parse time history load cases
    def _parse_load_case_time_history(self):
        '''
        command format:
        * LOAD_CASES_TH_NL_DIRECT_INTEGRATION
        # Name,MassSource,InitialCond,NonlinCase,LoadType,LoadName,Function,TransAccSF,TimeFactor,ArrivalTime,Angle,GeoNonlin,
        # NumSteps,StepSize,ProBy,MassCoeff,StiffCoeff,ProTimeVal1,ProDamping1,ProTimeVal2,ProDamping2,Mode4Ratio,ModalCase,
        # IntType,Gamma,Beta,Alpha,SolScheme,MaxSize,MinSize,EventTol,
        # MaxEvents,MaxCSIters,MaxNRIters,IterTol,LineSearch,MaxSearches,SearchTol,SearchFact,DesignType,GUID,Notes

        Now we just need:
        - Name
        - LoadType (Accel or TODO: ask kristijan if we need to support other types)
        - LoadName (U1, U2, U3, R1, R2, R3, direction, TODO: ask kristijan if the last 3 are correct)
        - Function (name of the time history function)
        - TransAccSF (scale factor for the acceleration time history function)
        - GeoNonlin (geometric nonlinearity)
        - NumSteps (number of steps)
        - StepSize (step size)
        - ProBy (Damping: proportional by, Direct, Period Ratio, or Frequency Ratio)
        - MassCoeff (mass coefficient, only if ProBy is Direct)
        - StiffCoeff (stiffness coefficient, only if ProBy is Direct)
        - ProTimeVal1 (proportional time value 1, only if ProBy is Period Ratio or Frequency Ratio)
        - ProDamping1 (proportional damping 1, only if ProBy is Period Ratio or Frequency Ratio)
        - ProTimeVal2 (proportional time value 2, only if ProBy is Period Ratio or Frequency Ratio)
        - ProDamping2 (proportional damping 2, only if ProBy is Period Ratio or Frequency Ratio)
        - Mode4Ratio (reference mode index for ratio, only if ProBy is Period Ratio or Frequency Ratio)
        '''
        # iterate over the load cases
        for item in self.commands['* LOAD_CASES_TH_NL_DIRECT_INTEGRATION']:
            words = item.split(',')
            name = words[0]
            # get existing or generare a new load case
            lc = self.doc.load_cases_dynamic.get(name, None)
            if lc is None:
                # parse other parameters
                load_type = words[4]
                kinematics = words[11]
                num_steps = int(words[12])
                step_size = float(words[13])
                pro_by = words[14]
                if pro_by == 'Direct':
                    mass_coeff = float(words[15])
                    stiff_coeff = float(words[16])
                    pro_time_val1 = 0.0
                    pro_damping1 = 0.0
                    pro_time_val2 = 0.0
                    pro_damping2 = 0.0
                    mode4_ratio = 0
                else:
                    mass_coeff = 0.0
                    stiff_coeff = 0.0
                    pro_time_val1 = float(words[17])
                    pro_damping1 = float(words[18])
                    pro_time_val2 = float(words[19])
                    pro_damping2 = float(words[20])
                    mode4_ratio = int(words[21])
                # create the load case
                lc = load_case_dynamic(name, load_type, num_steps, step_size,
                                    pro_by, mass_coeff, stiff_coeff, pro_time_val1, pro_damping1,
                                    pro_time_val2, pro_damping2, mode4_ratio)
                # add the load case to the document
                self.doc.load_cases_dynamic[name] = lc
                # check kinematics
                if self.doc.kinematics: # if already defined...
                    if kinematics != self.doc.kinematics:
                        self.interface.send_message(f'[LOAD_CASES_TH_NL_DIRECT_INTEGRATION {name}]: kinematics type {kinematics} does not match the existing kinematics type {self.doc.kinematics}', 
                                                    mtype=stko_interface.message_type.WARNING)
            # parse direction, function and scale factor
            load_name = words[5]
            function = words[6]
            scale_factor = float(words[7])
            # add the load pattern to the load case
            lc.functions.append((function, load_name, scale_factor))
