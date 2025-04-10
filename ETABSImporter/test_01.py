import sys
sys.path.append('C:/Develop/Utility/STKO-install/data/win64')

from PyMpc import *

import importlib
import ETABSImporter.parser
import ETABSImporter.builder
import ETABSImporter.stko_interface
import ETABSImporter.document
importlib.reload(ETABSImporter.document)
importlib.reload(ETABSImporter.parser)
importlib.reload(ETABSImporter.builder)
importlib.reload(ETABSImporter.stko_interface)

from ETABSImporter.parser import parser
from ETABSImporter.builder import builder
from ETABSImporter.stko_interface import stko_interface
#p = parser('Examples/Ex00/Model_Info.txt')
p = parser('Examples/Ex01/STKO_IMPORT_40_Story_Rc_Core_Wall_Bldng_e4.3.txt')

print(p.doc)

stko = stko_interface(App.caeDocument())
b = builder(p.doc, stko)



'''
Assumptions
areas are generic flat polygons
- they can be split in multiple lines
- holes use another command *AREA_OPENINGS_CONNECTIVITY
'''