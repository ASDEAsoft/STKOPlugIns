from ETABSImporter.document import *
from collections import defaultdict
import os
from PyMpc import *

# The parser class is used to parse the ETABS text file and build the document
class parser:
    
    def __init__(self, fname):
        # save the filename
        self.fname = fname
        # initialize document
        self.doc = document()
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
        self._parse_nodes()
        self._parse_frames()
        self._parse_areas()
        self._parse_diaphragm()
        self._parse_restraints()
        self._parse_load_patterns()
        self._parse_joint_loads()
        self._parse_joint_masses()
        self._parse_time_history_functions()
    
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
    def _parse_diaphragm(self):
        # todo: add the rigid diaphragm to the document only if it's rigid in ETABS (we need a flag)
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
            print(id, restraints)
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
            file_name = words[1]
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
            self.doc.th_functions[name] = th_function(name, values)