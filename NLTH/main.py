
import importlib
import NLTH.odb_utils.mpco_document_utils
import NLTH.gui.odb_widget
import NLTH.odb_utils.mpco_evaluator
import NLTH.gui.isolator_widget
import NLTH.gui.main_window

importlib.reload(NLTH.odb_utils.mpco_document_utils)
importlib.reload(NLTH.gui.odb_widget)
importlib.reload(NLTH.odb_utils.mpco_evaluator)
importlib.reload(NLTH.gui.isolator_widget)
importlib.reload(NLTH.gui.main_window)

from PySide2.QtWidgets import QApplication
from NLTH.gui.main_window import NLTHMainWindow

app = QApplication.instance()
if app is None:
    raise RuntimeError("No running QApplication instance found.")
mwin = app.activeWindow()
dialog = NLTHMainWindow(mwin)
dialog.exec_()