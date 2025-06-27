from MidasGenImporter.document import *
from MidasGenImporter.stko_interface import stko_interface
from MidasGenImporter.units_utils import unit_system
from MidasGenImporter.section_utils import get_section_offset, test_section
from collections import defaultdict
import os
import math
from PyMpc import *

def _split_line(contents:str, skip_empty:bool=True) -> List[str]:
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
