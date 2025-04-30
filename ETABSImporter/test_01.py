import sys
sys.path.append('C:/Develop/Utility/STKO-install/data/win64')

from PyMpc import *

import importlib
import ETABSImporter.parser
import ETABSImporter.builder
import ETABSImporter.stko_interface_plain
import ETABSImporter.stko_interface_gui
import ETABSImporter.document
importlib.reload(ETABSImporter.document)
importlib.reload(ETABSImporter.parser)
importlib.reload(ETABSImporter.builder)
importlib.reload(ETABSImporter.stko_interface_plain)
importlib.reload(ETABSImporter.stko_interface_gui)
from ETABSImporter.stko_interface_plain import stko_interface_plain
from ETABSImporter.stko_interface_gui import stko_interface_gui

# create the stko interface
# stko = stko_interface_plain(
#     doc = App.caeDocument(), 
#     etabs_filename = 'C:/Develop/STKOPlugIns/ETABSImporter/Examples/Ex01/STKO_IMPORT_40_Story_Rc_Core_Wall_Bldng_e4.3.txt')
stko = stko_interface_gui()
stko.perform()