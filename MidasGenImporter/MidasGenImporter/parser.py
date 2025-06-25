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
        # initialize command map
        self.commands = defaultdict(list)
        # build command map
        # TODO: handle mult-line commands
        with open(fname, 'r') as ifile:
            # first handle the multi-line commands (\ at the end of the line)
            # read all, replace \\n with \n, then split by \n)
            contents = ifile.read().replace('\\\n', '').replace('\\\r\n', '')
            lines = [pline for pline in (line.strip() for line in contents.splitlines()) if (pline and not pline.startswith(';'))]
            value = None
            for line in lines:
                if line.startswith('*'):
                    the_command = line.split(';')[0].strip()  # remove comments
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
        self._parse_groups()
        self._parse_elastic_materials()
        self._parse_frame_sections()
        
        self._parse_area_sections()
        self._parse_area_sections_assignment()
        
        self._parse_frame_sections_assignment()
        self._parse_frame_nonlinear_hinges()
        self._parse_frame_nonlinear_hinges_assignment()
        self._parse_diaphragms()
        self._parse_restraints()
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

    # called by the _parse_area_sections function if the MVLEM Notes is used
    def _parse_area_sections_mvlem(self, asec:area_section, notes:str):
        '''
        typical command:
        MVLEM 4 -thick 36*4 -width 9*4 -rho 0.0145*4 -matConcrete Conc_for_Walls*4 -matSteel Steel_for_Walls*4 -ShearK 1684540.78005847
        '''
        import re
        # first find all * chars, and if found, remove all spaces before and after it using regex
        notes = re.sub(r'\s*\*\s*', '*', notes)
        # now in a pattern like A*N, for each * found, get the string before (A) and after (N) it until a whitespace or end of string.
        # then replace the A*N with A A A ... A (N times)
        notes = re.sub(r'([^\s*]+)\*(\d+)', lambda m: ' '.join([m.group(1)] * int(m.group(2))), notes)
        # now we can split the notes by whitespaces (skipping empty strings)
        words = [i.strip() for i in notes.split() if i.strip()]
        nfibers = int(words[1])
        # each keyword starts with a dash.
        # let's create a dictionary with key = keyword (with dash) and value = list of strings found after that keyword
        keyword_values = {}
        keyword_values['type'] = words[0]  # the first word is the type, e.g. MVLEM
        keyword_values['nfibers'] = nfibers  # the number of fibers is the second word
        current_keyword = None
        for word in words[2:]:
            if word.startswith('-'):
                current_keyword = word
            elif current_keyword is not None:
                if current_keyword in ('-thick', '-width', '-rho', '-ShearK'):
                    word = float(word)
                values = keyword_values.get(current_keyword, None)
                if values is None:
                    values = []
                    keyword_values[current_keyword] = values
                values.append(word)
        # add it to the area section
        asec.conversion_info = keyword_values

    # this function parses the area section and adds it to the document as area_material     
    def _parse_area_sections(self):
        for is_wall in (False, True):
            cmd_name = 'AREA_SECTION_PROPERTIES_{}'.format('WALL' if is_wall else 'SLAB')
            for item in self.commands[f'* {cmd_name}']:
                words = item.split(',')
                name = words[0]
                type = words[1]
                material = words[2]
                try:
                    thickness = float(words[3])
                except:
                    thickness = 0.0
                f11 = float(words[4])
                f22 = float(words[5])
                f12 = float(words[6])
                m11 = float(words[7])
                m22 = float(words[8])
                m12 = float(words[9])
                v13 = float(words[10])
                v23 = float(words[11])
                mmod = float(words[12])
                wmod = float(words[13])
                # f11 f22 and f12 should be equal to each other
                # m11 m22 and m12 should be equal to each other
                # v13 and v23 should be one
                # mmod and wmod should be zero
                if (f11 != f22 or f11 != f12) and (self.interface is not None):
                    self.interface.send_message(f'[AREA_SECTION_PROPERTIES {name}]: f11, f22 and f12 should be equal to each other',
                                                 mtype=stko_interface.message_type.WARNING)
                if (m11 != m22 or m11 != m12) and (self.interface is not None):
                    self.interface.send_message(f'[AREA_SECTION_PROPERTIES {name}]: m11, m22 and m12 should be equal to each other', 
                                                mtype=stko_interface.message_type.WARNING)
                if (v13 != f11 or v23 != f11) and (self.interface is not None):
                    self.interface.send_message(f'[AREA_SECTION_PROPERTIES {name}]: v13 and v23 should be equal to f11', 
                                                mtype=stko_interface.message_type.WARNING)
                if (mmod != 1.0 or wmod != 1.0) and (self.interface is not None):
                    self.interface.send_message(f'[AREA_SECTION_PROPERTIES {name}]: mmod and wmod should be one', 
                                                mtype=stko_interface.message_type.WARNING)
                # creae the area material
                asec = area_section(name, type, material, thickness, f11, m11, is_wall=is_wall)
                # add the area material to the document
                self.doc.area_sections[name] = asec
                # parse notes if available
                if len(words) > 14:
                    notes = words[14].strip()
                    if notes.startswith('MVLEM'):
                        # parse MVLEM notes
                        self._parse_area_sections_mvlem(asec, notes)
                    else:
                        # just store the notes
                        self.interface.send_message(f'[{cmd_name} {name}]: unrecognized notes: {notes}',
                                                    mtype=stko_interface.message_type.WARNING)

    # this function parses the area section assignments and adds them to the document
    def _parse_area_sections_assignment(self):
        for item in self.commands['* AREA_SECTION_PROPERTIES_ASSIGNMENT']:
            words = item.split(',')
            area_id = int(words[0])
            material_name = words[1]
            # words[2] is not used here! wall or slab determined by the _parse_area_material command
            self.doc.area_sections_assignment[material_name].append(area_id)
            # fill the inverse map
            self.doc.area_sections_assignment_inverse[area_id] = material_name

    # this function parses the frame section assignments and adds them to the document
    def _parse_frame_sections_assignment(self):
        for item in self.commands['* FRAME_SECTION_PROPERTIES_ASSIGNMENT']:
            words = item.split(',')
            frame_id = int(words[0])
            section_name = words[1]
            self.doc.frame_sections_assignment[section_name].append(frame_id)
            # fill the inverse map
            self.doc.frame_sections_assignment_inverse[frame_id] = section_name

    # this function parses the frame nonlinear hinges and adds them to the document
    def _parse_frame_nonlinear_hinges(self):
        for item in self.commands['* FRAME_NONLINEAR_HINGE_PROPERTIES_V2']:
            words = [i for i in item.split(',') if i] # may contain empty strings
            name = words[0]
            # next values are a list of displacements and a list of forces. they must have the same length
            # and their length must be >= 2 and <= 7
            num_values = len(words) - 1
            # must be even
            if num_values % 2 != 0:
                raise Exception('Invalid (odd) number of values for frame nonlinear hinge {}: {}'.format(name, num_values))
            num_values //= 2
            if num_values < 2 or num_values > 7:
                raise Exception('Invalid (must be 2 <= N <= 7) number of values for frame nonlinear hinge {}: {}'.format(name, num_values))
            # read the forces (num_values)
            F = [float(words[i]) for i in range(1, num_values + 1)]
            # read the displacements (num_values)
            D = [float(words[i]) for i in range(num_values + 1, num_values * 2 + 1)]
            # some processing
            if F[0] <= 0.0:
                # skip the first value, the material in STKO assumed the first point not at origin
                F = F[1:]
                D = D[1:]
                self.interface.send_message(f'[FRAME_NONLINEAR_HINGE_PROPERTIES_V2 {name}]: the first value is not at origin, skipping it', mtype=stko_interface.message_type.WARNING)
                num_values -= 1
                if num_values < 2 or num_values > 7:
                    raise Exception('Invalid (must be 2 <= N <= 7) number of values for frame nonlinear hinge {}: {}'.format(name, num_values))
            # add
            self.doc.frame_nonlinear_hinges[name] = frame_nonlinear_hinge(name, D, F)

    # this function parses the frame nonlinear hinges assignments and adds them to the document
    def _parse_frame_nonlinear_hinges_assignment(self):
        for item in self.commands['* FRAME_NONLINEAR_HINGE_ASSIGNMENT']:
            words = item.split(',')
            frame_id = int(words[0])
            hinge_name = words[1]
            # words[2] is not used here! TODO: use relative distance once multi hinge are supported
            # add the hinge to the document
            self.doc.frame_nonlinear_hinges_assignment[hinge_name].append(frame_id)
            # fill the inverse map
            self.doc.frame_nonlinear_hinges_assignment_inverse[frame_id] = hinge_name

    # this function parses the rigid diaphragms and adds them to the document
    def _parse_diaphragms(self):
        # TODO: add the rigid diaphragm to the document only if it's rigid in ETABS (we need a flag)
        for item in self.commands['* JOINT_RIGID_DIAPHRAGMS']:
            words = item.split(',"')
            name = words[0]
            # the last one is the retained id
            items = [int(i.strip()) for i in words[1].replace('"', '').replace('[', '').replace(']', '').replace("'", '').split(',')]
            self.doc.diaphragms[name] = items

    # this function parses the restraints and adds them to the document
    def _parse_restraints(self):
        for item in self.commands['* JOINT_RESTRAINTS']:
            words = item.split(',')
            id = int(words[0])
            restraints = tuple([int(words[i]=='Yes') for i in range(1, 7)])
            self.doc.restraints[id] = restraints
    
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
            self.doc.joint_loads[id] = joint_load(load_pattern, value)
    
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
