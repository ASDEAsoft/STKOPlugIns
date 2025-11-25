from MidasGenImporter.document import *
from MidasGenImporter.stko_interface import stko_interface
from MidasGenImporter.units_utils import unit_system
from MidasGenImporter.section_utils import get_section_offset_and_area
from MidasGenImporter.localaxes_utils import MIDASLocalAxesConvention
from MidasGenImporter.area_load_utils import floor_polygon, floor_polygon_generator
from MidasGenImporter.beam_load_utils import compute_equivalent_constant_load
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
        # if a frame overlaps with a fully rigid elastic link, we need to remove the frame.
        # with this map, we can get the link id from the removed frame id (useful for loads)
        self.removed_frames_to_elastic_link : Dict[int, int] = dict()
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
        self._parse_masses()

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

        # if we have fully rigid elastic links overlapping with frames, we need to remove the frames
        frames_to_remove = []

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
            # check overlapping with existing frames
            for old_frame_id, old_frame in self.doc.frames.items():
                if old_frame.link is None: # if old frame is not a link...
                    if (old_frame.nodes[0] == node1 and old_frame.nodes[1] == node2) or (old_frame.nodes[0] == node2 and old_frame.nodes[1] == node1):
                        # overlapping with an existing frame, make sure it is fully rigid
                        if not all(elastic_link.link.rigid_flags):
                            raise Exception(f'Elastic link {id} is not fully rigid ({elastic_link.link.rigid_flags}), but it overlaps with frame {old_frame_id} which is elastic')
                        # add the frame to the list to remove
                        frames_to_remove.append(old_frame_id)
                        self.removed_frames_to_elastic_link[old_frame_id] = id
            # increment the ID for the next elastic link
            id += 1
        if self.interface is not None:
            self.interface.send_message(f'Parsed {len(self.doc.frames)} elastic links', mtype=stko_interface.message_type.INFO)
        # if we have frames to remove, we need to remove them from the document
        if len(frames_to_remove) > 0:
            if self.interface is not None:
                self.interface.send_message(f'Removing {len(frames_to_remove)} frames that overlap with fully rigid elastic links', mtype=stko_interface.message_type.WARNING)
            # remove the frames from the document
            for frame_id in frames_to_remove:
                del self.doc.frames[frame_id]
            # remove the frame ids from the section_scale_factors
            for section_sf in self.doc.section_scale_factors:
                section_sf.elements = [elem for elem in section_sf.elements if elem not in frames_to_remove]

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
                offset_y, offset_z, area = get_section_offset_and_area(shape, shape_info, offset_data)
                self.doc.sections[sec_id] = section(sec_name, shape, shape_info, offset_y, offset_z, area)
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
        map_mod_to_ele : DefaultDict[Tuple[float,float,float,float,float,float], List[int]] = defaultdict(list)
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
            # map it so that we don't have multiple section scale factors with the same values
            map_mod_to_ele[(area_sf, asy_sf, asz_sf, ixx_sf, iyy_sf, izz_sf)].extend(elem_list)
        # create section scale factors objects from unique values
        for (area_sf, asy_sf, asz_sf, ixx_sf, iyy_sf, izz_sf), elem_list in map_mod_to_ele.items():
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
        map_mod_to_ele : DefaultDict[Tuple[float,float], List[int]] = defaultdict(list)
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
            # map it so that we don't have multiple section scale factors with the same values
            map_mod_to_ele[(in_plane, out_of_plane)].extend(elem_list)
        # create thickness scale factors objects from unique values
        for (in_plane, out_of_plane), elem_list in map_mod_to_ele.items():
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
            # map it so that we don't have multiple section scale factors with the same values
            map_mod_to_ele[(in_plane, out_of_plane)].extend(elem_list)
        # create thickness scale factors objects from unique values
        for (in_plane, out_of_plane), elem_list in map_mod_to_ele.items():
            self.doc.thickness_scale_factors.append(thickness_scale_factors(in_plane, out_of_plane, elements=elem_list))
        # print
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
        def _get_frame_or_link_id(elem_id: int) -> int:
            link_id = self.removed_frames_to_elastic_link.get(elem_id, None)
            if link_id is not None:
                return link_id
            return elem_id

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
                # NOTE: this plain conversion is fine IF we have constant uniform loads only!!
                mgt_ele_list = _mgt_string_to_id_or_range(words[0])
                elem_list = [] # for beam load
                link_list = [] # for equivalent nodal load (they don't support beam loads on elastic links)
                for ele_id in mgt_ele_list:
                    link_id = self.removed_frames_to_elastic_link.get(ele_id, None)
                    if link_id is not None:
                        link_list.append(link_id)
                    else:
                        elem_list.append(ele_id)
                cmd = words[1]
                if cmd != 'BEAM' and cmd != 'FLOOR':
                    raise Exception('Invalid beam load command: {}, expecting BEAM or FLOOR'.format(cmd))
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
                load_component = compute_equivalent_constant_load([float(v) for v in words[10:18]])
                # make load vector
                load_value = (load_component if component == 0 else 0.0,
                              load_component if component == 1 else 0.0,
                              load_component if component == 2 else 0.0)
                # create the beam load object
                beam_load_obj = beam_load(load_value, is_local)
                # add the load to the load case
                lc.beam_loads[beam_load_obj].extend(elem_list)  # add the element IDs to the list of elements for this load
                # now, if we have link elements, we need to convert the beam load to equivalent nodal loads
                if len(link_list) > 0:
                    if is_local:
                        raise Exception('Beam loads on elastic links cannot be in local direction, found local beam load on links')
                    for link_id in link_list:
                        link_elem = self.doc.frames.get(link_id, None)
                        if link_elem is None:
                            raise Exception(f'Elastic link ID: {link_id} not found in document')
                        n1, n2 = link_elem.nodes
                        node1 = self.doc.vertices.get(n1, None)
                        node2 = self.doc.vertices.get(n2, None)
                        if node1 is None or node2 is None:
                            raise Exception(f'Elastic link ID: {link_id} has invalid node IDs: {n1}, {n2}')
                        L = (node2 - node1).norm()
                        # compute equivalent nodal loads for the link element
                        eqLoad = load_component * L / 2.0
                        fx = eqLoad if component == 0 else 0.0
                        fy = eqLoad if component == 1 else 0.0
                        fz = eqLoad if component == 2 else 0.0
                        # create the nodal load object
                        nodal_load_obj = nodal_load((fx, fy, fz, 0.0, 0.0, 0.0))
                        # add it to the load case
                        lc.nodal_loads[nodal_load_obj].extend([n1, n2])  # add the node IDs to the list of nodes for this load
                
        if self.interface is not None:
            for lc_name, lc in self.doc.load_cases.items():
                if len(lc.beam_loads) > 0:
                    self.interface.send_message(f'Parsed {len(lc.beam_loads)} beam loads for load case {lc_name}', mtype=stko_interface.message_type.INFO)

        # beam load summary
        print('Beam Loads Summary:')
        Q_SUM_ALL = 0.0
        N_SUM_ALL = 0.0
        for lc_name, lc in self.doc.load_cases.items():
            print(f' Load Case: {lc_name}')
            Q_SUM = 0.0 # sum load per unit length
            N_SUM = 0.0 # sum of total load from beam loads
            for beam_load_obj, elem_list in lc.beam_loads.items():
                value = beam_load_obj.value[2]
                for elem_id in elem_list:
                    elem = self.doc.frames.get(elem_id, None)
                    if elem is None:
                        raise Exception(f'  Element ID: {elem_id} not in document frames')
                    n1,n2 = elem.nodes
                    node1 = self.doc.vertices.get(n1, None)
                    node2 = self.doc.vertices.get(n2, None)
                    if node1 is None or node2 is None:
                        raise Exception(f'  Element ID: {elem_id} has invalid node IDs: {n1}, {n2}')
                    L = (node2 - node1).norm()
                    Q_SUM += value
                    N_SUM += value * L
            print(f'    Total uniform load per unit length (Z dir): {Q_SUM} kN/m')
            print(f'    Total uniform load from beam loads (Z dir): {N_SUM} kN')
            Q_SUM_ALL += Q_SUM
            N_SUM_ALL += N_SUM
        print(f' Overall total uniform load per unit length (Z dir): {Q_SUM_ALL} kN/m')
        print(f' Overall total uniform load from beam loads (Z dir): {N_SUM_ALL} kN')
                    
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
                if is_local and component != 2:
                    raise Exception('Invalid pressure load direction: {}, expecting LZ for local Z direction'.format(direction))
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
                    pu if component == 0 else 0.0,  # X component
                    pu if component == 1 else 0.0,  # Y component
                    pu if component == 2 else 0.0,  # Z component
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
            points: Dict[int, Math.vec3] = {}
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
                points[node_id] = vertex
            bbox_diag = bbox.maxPoint - bbox.minPoint
            bbox_diag_avg = (bbox_diag.x + bbox_diag.y + bbox_diag.z) / 3.0
            tolerance = 1.0e-12*bbox_diag_avg
            if z_error > tolerance:
                raise Exception(f'Nodes {node_ids} are not at the same Z coordinate, max error = {z_error}, tolerance = {tolerance}')
            # now we can create the polygon
            polygon_array = floor_polygon_generator(points, self.doc)
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
                # now we have to 1) compute the distributed load (F/L) for each edge
                # and 2) lump it to the 2 end-nodes of the edge
                # so we can create equivalent nodal loads
                eq_nodal_loads = defaultdict(float)
                for polygon in polygon_array:
                    total_load = area_load * polygon.area
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
                self.interface.send_message(f'Load case {lc_name} has {num_eq_nodal_loads} equivalent nodal loads from floor loads', 
                                            mtype=stko_interface.message_type.INFO)

    # this function checks if the load cases are valid
    def _check_load_cases(self):
        for lc_name, lc in self.doc.load_cases.items():
            if self.interface is not None:
                self.interface.send_message(f'Checking load case {lc_name}', mtype=stko_interface.message_type.INFO)
            lc.check()
    
    # parses all the mass sources and makes them nodal masses
    def _parse_masses(self):
        # nodal masses, then we group them
        node_mass : DefaultDict[int, Math.vec3] = defaultdict(Math.vec3)
        # first build the masses from densities.
        # parse frame elements
        ele_length_map : Dict[int, float] = {}
        for ele_id, ele in self.doc.frames.items():
            # skip link elements
            # UPDATE: we now consider link elements for mass calculation (in case beam loads are applied on them)
            #if ele.link is not None:
            #    continue
            # get element length
            n1, n2 = ele.nodes
            node1 = self.doc.vertices.get(n1, None)
            node2 = self.doc.vertices.get(n2, None)
            if node1 is None or node2 is None:
                raise Exception(f'Frame element {ele_id} has invalid nodes: {n1}, {n2}')
            length = (node2 - node1).norm()
            ele_length_map[ele_id] = length
            # now we can skip link elements (we only need to keep track of their length)
            if ele.link is not None:
                continue
            # get material density
            mat:elastic_material = self.doc.elastic_materials.get(ele.mat, None)
            if mat is None:
                raise Exception(f'Frame element {ele_id} has invalid material: {ele.mat}')
            density = mat.rho
            # get the section area
            sec:section = self.doc.sections.get(ele.sec, None)
            if sec is None:
                raise Exception(f'Frame element {ele_id} has invalid section: {ele.sec}')
            area = sec.area
            # lumped mass
            mass = abs(density * length * area / 2.0)
            mass = Math.vec3(mass, mass, mass)
            # add the mass to the 2 nodes of the element
            node_mass[n1] += mass
            node_mass[n2] += mass
        # parse area elements
        ele_area_map : Dict[int, float] = {}
        for ele_id, ele in self.doc.areas.items():
            # get nodes
            nodes = []
            for node_id in ele.nodes:
                node = self.doc.vertices.get(node_id, None)
                if node is None:
                    raise Exception(f'Area element {ele_id} has invalid node: {node_id}')
                nodes.append(node)
            # compute element area (they can have 3 or 4 nodes)
            p01 = nodes[1] - nodes[0]  # vector from node 0 to node 1
            p02 = nodes[2] - nodes[0]  # vector from node 0 to node 2
            area = p01.cross(p02).norm() / 2.0
            if len(nodes) == 4:
                p03 = nodes[3] - nodes[0]
                area += p02.cross(p03).norm() / 2.0
            ele_area_map[ele_id] = area
            # get material density
            mat:elastic_material = self.doc.elastic_materials.get(ele.mat, None)
            if mat is None:
                raise Exception(f'Area element {ele_id} has invalid material: {ele.mat}')
            density = mat.rho
            # get the section thickness
            sec:thickness = self.doc.thicknesses.get(ele.sec, None)
            if sec is None:
                raise Exception(f'Area element {ele_id} has invalid section: {ele.sec}')
            thick = sec.in_thick
            # lumped mass
            mass = abs(density * area * thick / len(nodes))  # divide by number of nodes to distribute the mass evenly
            mass = Math.vec3(mass, mass, mass)
            # add the mass to the nodes of the element
            for node_id in ele.nodes:
                node_mass[node_id] += mass
        print(f'Parsed {len(node_mass)} nodal masses from densities')


        STRUCT_MASS = 0.0
        L2M_NODE = 0.0
        L2M_BEAM = 0.0
        L2M_BEAM_LOAD = 0.0
        L2M_BEAM_LOAD_DIST = 0.0
        L2M_FLOOR = 0.0
        L2M_PRES = 0.0
        for _, mass in node_mass.items():
            STRUCT_MASS += mass.z  # sum Z component as total mass

        # copy this dict so keep only densities
        node_mass_dens:Dict[int,float] = {node_id : mass.x for node_id, mass in node_mass.items()}
        # now parse the load to masses
        '''
        *LOADTOMASS    ; Load to Mass
        ; DIR, bNODAL, bBEAM, bFLOOR, bPRES, GRAV
        ; LCNAME1, FACTOR1, LCNAME2, FACTOR2, ...   ; from line 1
        '''
        # scale factors
        direction = Math.vec3(0.0, 0.0, 0.0)
        b_nodal = False
        b_beam = False
        b_floor = False
        b_pres = False
        grav = 0.0
        for line in self.commands['*LOADTOMASS']:
            words = _split_line(line, skip_empty=False)
            lc_name = words[0]
            if lc_name not in self.doc.load_cases:
                # this is a new line DIR, bNODAL, bBEAM, bFLOOR, bPRES, GRAV
                if len(words) != 6:
                    raise Exception('Invalid LOADTOMASS line: {}, expecting 6 values'.format(line))
                xyz = words[0]
                direction.x = 1.0 if 'X' in xyz else 0.0
                direction.y = 1.0 if 'Y' in xyz else 0.0
                direction.z = 1.0 if 'Z' in xyz else 0.0
                b_nodal = words[1] == 'YES'
                b_beam = words[2] == 'YES'
                b_floor = words[3] == 'YES'
                b_pres = words[4] == 'YES'
                grav = float(words[5])
            else:
                # this is a new line CNAME1, FACTOR1, LCNAME2, FACTOR2, ...
                if len(words) % 2 != 0:
                    raise Exception('Invalid LOADTOMASS line: {}, expecting an even number of values'.format(line))
                # parse the load cases and factors
                for i in range(0, len(words), 2):
                    lc_name = words[i]
                    factor = float(words[i + 1])
                    lc = self.doc.load_cases.get(lc_name, None)
                    if lc is None:
                        raise Exception(f'Load case {lc_name} not found for LOADTOMASS')
                    # now we can convert the load in Z direction to mass:
                    # mass_vector = load_z:float * factor / grav * direction (*length or area for beam loads or pressure loads)
                    if b_nodal:
                        for load, nodes in lc.nodal_loads.items():
                            Fz = load.value[2]
                            if Fz == 0.0: continue
                            mass = abs(Fz * factor / grav) * direction
                            for node_id in nodes:
                                node_mass[node_id] += mass
                                L2M_NODE += mass.z
                    if b_floor:
                        for load, nodes in lc.floor_loads.items():
                            Fz = load.value[2]
                            if Fz == 0.0: continue
                            mass = abs(Fz * factor / grav) * direction
                            for node_id in nodes:
                                node_mass[node_id] += mass
                                L2M_FLOOR += mass.z
                    if b_beam:
                        for load, elems in lc.beam_loads.items():
                            if load.local:
                                raise Exception('Local beam loads are not supported for LOADTOMASS')
                            Fz = load.value[2]
                            if Fz == 0.0: continue
                            # get the length of the element
                            for elem_id in elems:
                                ele = self.doc.frames.get(elem_id, None)
                                if ele is None:
                                    raise Exception(f'Element {elem_id} not found for LOADTOMASS')
                                length = ele_length_map.get(elem_id, None)
                                if length is None:
                                    raise Exception(f'Length of element {elem_id} not found for LOADTOMASS')
                                mass = abs(Fz * factor / grav) * direction * length / len(ele.nodes)  # divide by number of nodes to distribute the mass evenly
                                # get the nodes of the element
                                for ele_node in ele.nodes:
                                    node_mass[ele_node] += mass
                                    L2M_BEAM += mass.z
                                L2M_BEAM_LOAD_DIST += Fz
                                L2M_BEAM_LOAD += Fz * length
                    if b_pres:
                        for load, elems in lc.pressure_loads.items():
                            # get the area of the element
                            for elem_id in elems:
                                # compute the component in the global Z direction
                                if load.local:
                                    ele = self.doc.areas.get(elem_id, None)
                                    T = MIDASLocalAxesConvention.get_local_axes([self.doc.vertices[i] for i in ele.nodes])  # this is just to check if the load is local
                                    global_value = T * Math.vec3(*load.value)
                                    Fz = global_value.z  # global Z direction
                                else:
                                    Fz = load.value[2]  # global Z direction
                                if Fz == 0.0: continue
                                area = ele_area_map.get(elem_id, None)
                                if area is None:
                                    raise Exception(f'Area of element {elem_id} not found for LOADTOMASS')
                                mass = abs(Fz * factor / grav) * direction * area / len(ele.nodes)  # divide by number of nodes to distribute the mass evenly
                                # get the nodes of the element
                                ele = self.doc.areas.get(elem_id, None)
                                if ele is None:
                                    raise Exception(f'Element {elem_id} not found for LOADTOMASS')
                                for ele_node in ele.nodes:
                                    node_mass[ele_node] += mass
                                    L2M_PRES += mass.z

        # merge nodal masses with same value
        for node_id, mass in node_mass_dens.items():
            self.doc.masses_dens[mass].append(node_id)  # add the node ID to the list of nodes for this mass
        for node_id, mass in node_mass.items():
            key = (mass.x, mass.y, mass.z)  # use a tuple as key to avoid duplicates
            self.doc.masses[key].append(node_id)  # add the node ID to the list of nodes for this mass
        print(f'reduced to {len(self.doc.masses)} grouped nodal masses from densities and load-to-mass')

        print('Mass summary:')
        print(f'  Structural Mass from densities: {STRUCT_MASS} tons')
        print(f'  Load to Mass from Nodal Loads: {L2M_NODE} tons')
        print(f'  Load to Mass from Beam Loads: {L2M_BEAM} tons (from total load {L2M_BEAM_LOAD} kN) (from distributed load {L2M_BEAM_LOAD_DIST} kN/m)')
        print(f'  Load to Mass from Floor Loads: {L2M_FLOOR} tons')
        print(f'  Load to Mass from Pressure Loads: {L2M_PRES} tons')
        total_mass = STRUCT_MASS + L2M_NODE + L2M_BEAM + L2M_FLOOR + L2M_PRES
        print(f'  Total Mass: {total_mass} tons')
        





