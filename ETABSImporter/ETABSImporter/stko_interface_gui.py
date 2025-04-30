from PyMpc import *
from ETABSImporter.stko_interface import stko_interface
from ETABSImporter.parser import parser
from ETABSImporter.builder import builder
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
	QApplication,
	QDialog,
	QFileDialog,
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
	# print
	def sendLog(self, perc, info):
		self.pbar.setValue(max(0, min(100, int(perc*100.0))))
		self.log.appendPlainText(info)

# ask the user for the input file
def _selectInputFile():
    filename = QFileDialog.getOpenFileName(
        QApplication.activeWindow(), "Open ETABS File", ".", 
        "Text files (*.txt);;All files (*.*)")
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
            self.dialog.sendLog(0.0, 'Parsing input file...')
            the_parser = parser(self.etabs_filename)
            the_parser.doc.process()
            # build the document
            self.dialog.sendLog(0.5, 'Building document...')
            the_builder = builder(the_parser.doc, self)
            # run a regenerate command
            self.dialog.sendLog(0.9, 'Regenerating document...')
            self.regenerate()
            self.dialog.sendLog(1.0, 'Done!')
        except Exception as e:
            raise e
        finally:
            # stop the STKO interface
            self.stop()
            if self.dialog is not None:
                # close the dialog
                self.dialog.closeNow()
                self.dialog.deleteLater()
    
    # send a message to the terminal window
    def send_message(self, msg : str):
        print(msg)
        if self.dialog is not None:
            self.dialog.sendLog(0.0, msg)
        App.processEvents()
