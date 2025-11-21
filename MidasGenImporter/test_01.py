import sys
sys.path.append('C:/Develop/Utility/STKO-install/data/win64')
sys.path.append(r'C:\Python310\Lib\site-packages')

from PyMpc import *

import importlib
import MidasGenImporter.area_load_utils
importlib.reload(MidasGenImporter.area_load_utils)
import MidasGenImporter.parser
import MidasGenImporter.builder
import MidasGenImporter.stko_interface
import MidasGenImporter.stko_interface_plain
import MidasGenImporter.stko_interface_gui
import MidasGenImporter.document
import MidasGenImporter.section_utils
import MidasGenImporter.localaxes_utils
importlib.reload(MidasGenImporter.document)
importlib.reload(MidasGenImporter.section_utils)
importlib.reload(MidasGenImporter.localaxes_utils)
importlib.reload(MidasGenImporter.parser)
importlib.reload(MidasGenImporter.builder)
importlib.reload(MidasGenImporter.stko_interface)
importlib.reload(MidasGenImporter.stko_interface_plain)
importlib.reload(MidasGenImporter.stko_interface_gui)
from MidasGenImporter.stko_interface_plain import stko_interface_plain
from MidasGenImporter.stko_interface_gui import stko_interface_gui

App.setColorWorkTreeByUsageFlag(False)

# create the stko interface
# stko = stko_interface_plain(
#     doc = App.caeDocument(), 
#     etabs_filename = 'C:/Develop/STKOPlugIns/ETABSImporter/Examples/Ex01/STKO_IMPORT_40_Story_Rc_Core_Wall_Bldng_e4.3.txt')
stko = stko_interface_gui()
stko.perform()