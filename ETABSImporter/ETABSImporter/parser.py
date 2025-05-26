from ETABSImporter.document import *
from ETABSImporter.stko_interface import stko_interface
from collections import defaultdict
import os
import math
from PyMpc import *

# The parser class is used to parse the ETABS text file and build the document
class parser:
    
    def __init__(self, fname, interface:stko_interface=None):
        # save the filename
        self.fname = fname
        # initialize document
        self.doc = document()
        # save the interface
        self.interface = interface
        # initialize command map
        self.commands = defaultdict(list)
        # build command map
        with open(fname, 'r') as ifile:
            lines = [pline for pline in (line.strip() for line in ifile.readlines()) if (pline and not pline.startswith('#'))]
            value = None
            for line in lines:
                if line.startswith('*'):
                    value = self.commands[line]
                else:
                    if value is not None:
                        value.append(line)
        # parse
        self._parse()
    
    def _parse(self):
        self._parse_units()
        self._parse_nodes()
        self._parse_frames()
        self._parse_areas()
        self._parse_elastic_materials()
        self._parse_nonlinear_materials()
        self._parse_area_sections()
        self._parse_area_sections_assignment()
        self._parse_frame_sections()
        self._parse_frame_sections_assignment()
        self._parse_frame_nonlinear_hinges()
        self._parse_frame_nonlinear_hinges_assignment()
        self._parse_diaphragms()
        self._parse_restraints()
        self._parse_load_patterns()
        self._parse_joint_loads()
        self._parse_joint_masses()
        self._parse_time_history_functions()
    
    def _parse_units(self):
        # this function parses the units and adds them to the document
        # the units are stored in the document as a dictionary
        for item in self.commands['* UNITS']:
            words = tuple(i.strip() for i in item.split(','))
            if len(words) != 3:
                raise Exception('Invalid units line: {}, expecting 3 values'.format(item))
            self.doc.units = words

    # this function parses the nodes and adds them to the document
    # the nodes are stored in the document as vertices
    def _parse_nodes(self):
        for item in self.commands['* JOINT_COORDINATES']:
            words = item.split(',')
            id = int(words[0])
            x = float(words[1])
            y = float(words[2])
            z = float(words[3])
            self.doc.vertices[id] = Math.vec3(x,y,z)
    
    # this function parses the frames and adds them to the document
    def _parse_frames(self):
        for item in self.commands['* FRAMES_CONNECTIVITY']:
            words = item.split(',')
            id = int(words[0])
            nodes = [int(words[i]) for i in range(1, 3)]
            angle = float(words[-1])
            self.doc.frames[id] = frame(nodes, angle)
    
    # this function parses the areas and adds them to the document
    def _parse_areas(self):
        for item in self.commands['* AREAS_CONNECTIVITY']:
            words = item.split(',')
            id = int(words[0])
            nodes = [int(words[i]) for i in range(1, 5)]
            angle = float(words[-1])
            self.doc.areas[id] = area(nodes, angle)
    
    # this function parses the rigid diaphragms and adds them to the document
    def _parse_diaphragms(self):
        # TODO: add the rigid diaphragm to the document only if it's rigid in ETABS (we need a flag)
        for item in self.commands['* JOINT_RIGID_DIAPHRAGMS']:
            words = item.split(',"')
            name = words[0]
            # the last one is the retained id
            items = [int(i.strip()) for i in words[1].replace('"', '').replace('[', '').replace(']', '').replace("'", '').split(',')]
            self.doc.diaphragms[name] = items
    
    # this function parses the elastic material and adds it to the document as elastic_material
    def _parse_elastic_materials(self):
        for item in self.commands['* MATERIALS_ELASTIC_PROPERTIES']:
            words = item.split(',')
            name = words[0]
            density_type = words[1]
            unit_weight = float(words[2])
            unit_mass = float(words[3])
            e1 = float(words[4])
            try:
                g12 = float(words[5])
                u12 = float(words[6])
            except:
                g12 = 0.0
                u12 = 0.0
            a1 = float(words[7])
            # add the material to the document
            self.doc.elastic_materials[name] = elastic_material(name, density_type, unit_weight, unit_mass, e1, g12, u12, a1)

    # this function parses the nonlinear materials and adds it to the document as nonlinear_material
    def _parse_nonlinear_materials(self):
        # store a temporary dictionary of nonlinear materials with point,strain,stress
        # then we will order the points
        temp : Dict[Tuple[str,str], List[Tuple[float,float]]] = {}
        for item in self.commands['* MATERIALS_NONLINEAR_PROPERTIES']:
            words = item.split(',')
            name = words[0]
            mat_type = words[1]
            key = (name, mat_type)
            strain = float(words[3])
            stress = float(words[4])
            item = temp.get(key, None)
            if item is None:
                item = []
                temp[key] = item
            item.append((strain, stress))
        # now we have to:
        # 2. split them into positive and negative (and make them all positive)
        # 3. sort them
        # 3. create the nonlinear material
        # 4. add it to the document
        for (name, mat_type), points in temp.items():
            # split the points into positive and negative
            pos_points = []
            neg_points = []
            for strain, stress in points:
                if strain == 0.0:
                    # skip the zero point
                    continue
                if strain > 0.0:
                    pos_points.append((strain, stress))
                else:
                    neg_points.append((-strain, -stress)) # make them positive
            # sort the points
            pos_points.sort(key=lambda x: x[0])
            neg_points.sort(key=lambda x: x[0])
            # create the nonlinear material
            nmat = nonlinear_material(name, mat_type, pos_points, neg_points)
            # add it to the document
            self.doc.nonlinear_materials[name] = nmat

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

    # this function parses the frame sections and adds them to the document
    def _parse_frame_sections(self):
        # first parse shape information if available
        # key is the name of the frame section, value is a tuple of shape type and shape info data
        all_shape_info : Dict[str, Tuple[frame_section.shape_type, List[float]]] = {}
        # TODO: we assume t2,t3 but ETABS has a bug on it
        # rectangular
        for item in self.commands['* FRAME_SECTION_DIMENSIONS_RECTANGULAR']:
            words = item.split(',')
            name = words[0]
            values = [float(i) for i in words[1:]]
            all_shape_info[name] = (frame_section.shape_type.rectangle, values)
        # TODO: add other shapes
        # now parse the frame sections
        for item in self.commands['* FRAME_SECTION_PROPERTIES']:
            # this is the line, use it for parsing:
            # Name,Material,Shape,Area,J,I33,I22,I23,As2,As3,S33Pos,S33Neg,S22Pos,S22Neg,Z33,Z22,R33,R22,CGOffset3,CGOffset2,PNAOffset3,PNAOffset2,AMod,A2Mod,A3Mod,JMod,I3Mod,I2Mod,MMod,WMod
            # Note: some of the values are not used in the document, skipped here...
            words = item.split(',')
            name = words[0]
            material = words[1]
            _shape = words[2]
            A = float(words[3])
            J = float(words[4])
            I33 = float(words[5])
            I22 = float(words[6])
            As2 = float(words[8])
            As3 = float(words[9])
            CGOffset3 = float(words[18])
            CGOffset2 = float(words[19])
            AMod = float(words[22])
            A2Mod = float(words[23])
            A3Mod = float(words[24])
            JMod = float(words[25])
            I3Mod = float(words[26])
            I2Mod = float(words[27])
            # check shape
            # TODO: wait for FRAME_SECTION_DIMENSIONS_RECTANGULAR
            shape = frame_section.shape_type.generic
            shape_info = None
            shape_data = all_shape_info.get(name, None)
            if shape_data is not None:
                shape, shape_info = shape_data
            # compute shear correction factors
            S2 = As2 / A if A > 0.0 else 0.0
            S3 = As3 / A if A > 0.0 else 0.0
            # create the frame section
            fsec = frame_section(name, shape, material, A, I22, I33, J, S2, S3, CGOffset2, CGOffset3, AMod, A2Mod, A3Mod, I2Mod, I3Mod, JMod, shape_info=shape_info)
            # add the frame section to the document
            self.doc.frame_sections[name] = fsec

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