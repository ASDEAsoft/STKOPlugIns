from PyMpc import *
from typing import List
from MidasGenImporter.stko_interface import stko_interface
from MidasGenImporter.parser import parser
from MidasGenImporter.builder import builder
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
	QApplication,
	QDialog,
	QFileDialog,
    QInputDialog,
	QVBoxLayout,
	QLabel,
	QProgressBar,
	QPlainTextEdit)

# A custom waiting dialog
class _WDialog(QDialog):
    def __init__(self, parent=None):
        QDialog.__init__(self, parent)
        self.canClose = False
        # build dialog
        self.setWindowTitle('Importing Tcl...')
        lo = QVBoxLayout(self)
        self.label = QLabel('Process status:')
        lo.addWidget(self.label)
        self.pbar = QProgressBar(self)
        self.pbar.setMaximum(100)
        self.pbar.setTextVisible(True)
        self.pbar.setValue(0)
        lo.addWidget(self.pbar)
        self.log = QPlainTextEdit(self)
        self.log.setReadOnly(True)
        self.log.setOverwriteMode(False)
        self.log.setTextInteractionFlags(Qt.NoTextInteraction)
        self.log.setUndoRedoEnabled(False)
        self.log.setBackgroundVisible(False)
        lo.addWidget(self.log)
        self.setLayout(lo)
    # handle the non-closable stuff
    def closeNow(self):
        self.canClose = True
        self.accept()
        self.canClose = False
    def accept(self):
        if self.canClose:
            QDialog.accept(self)
    def reject(self):
        if self.canClose:
            QDialog.reject(self)
    # print and percentage
    def sendPercentage(self, perc : float):
        self.pbar.setValue(max(0, min(100, int(perc*100.0))))
    def sendLog(self, info : str):
        self.log.appendPlainText(info)

# ask the user for the input file
def _selectInputFile():
    filename = QFileDialog.getOpenFileName(
        QApplication.activeWindow(), "Open Midas Gen File", ".", 
        "MGT files (*.mgt);;All files (*.*)")
    if filename is not None:
        filename = filename[0]
    return filename

# the GUI interface to the STKO document
class stko_interface_gui(stko_interface):
    def __init__(
			self,
            doc : MpcCaeDocument = None,
            etabs_filename : str = None):
        super().__init__(doc=doc, etabs_filename=etabs_filename)
        self.dialog = None
        self.warnings : List[str] = []
        self.errors : List[str] = []

    def perform(self):
        try:
            # start the STKO interface
            self.start()
            # check document
            if self.doc is None:
                App.runCommand('NewPreDocument')
                self.doc = App.caeDocument()
            # parse the input file
            if self.etabs_filename is None:
                self.etabs_filename = _selectInputFile()
            if self.etabs_filename is None:
                raise Exception('No input file selected.')
            # build the dialog
            App.clearTerminal()
            self.dialog = _WDialog(QApplication.activeWindow())
            self.dialog.show()
            # parse the input file
            self.send_percentage(0.0)
            self.send_message('Parsing input file...')
            the_parser = parser(self.etabs_filename, interface=self)
            the_parser.doc.process()
            # build the document
            self.send_percentage(0.5)
            self.send_message('Building document...')
            the_builder = builder(the_parser.doc, self)
            # run a regenerate command
            self.send_percentage(0.9)
            self.send_message('Regenerating document...')
            self.regenerate()
            self.send_percentage(1.0)
            self.send_message('Done!')
        except Exception as e:
            raise e
        finally:
            # stop the STKO interface
            self.stop()
            if self.dialog is not None:
                # close the dialog
                self.dialog.closeNow()
                self.dialog.deleteLater()
    
    # send a percentage to the dialog
    def send_percentage(self, percentage : float):
        if self.dialog is not None:
            self.dialog.sendPercentage(percentage)
        # process events
        App.processEvents()

    # send a message to the terminal window
    def send_message(self, msg : str, mtype : stko_interface.message_type = stko_interface.message_type.INFO):
        print(msg)
        if self.dialog is not None:
            self.dialog.sendLog(msg)
        if mtype == stko_interface.message_type.WARNING:
            self.warnings.append(msg)
        elif mtype == stko_interface.message_type.ERROR:
            self.errors.append(msg)
        # process events
        App.processEvents()
    
    # let the user select an input from a list
    def select_from_list(self, title : str, items : List[str], default_index : int = 0) -> int:
        if len(items) == 0:
            return -1
        if default_index < 0 or default_index >= len(items):
            default_index = 0
        # this default implementation is for GUI interfaces
        index, ok = QInputDialog.getItem(
            QApplication.activeWindow(), title, 
            "Select an item:", items, default_index, False)
        if ok and index:
            return items.index(index)
        return -1
