from PyMpc import *
from ETABSImporter.stko_interface import stko_interface
from ETABSImporter.parser import parser
from ETABSImporter.builder import builder

# the plain interface to the STKO document
class stko_interface_plain(stko_interface):
    def __init__(
			self,
            doc : MpcCaeDocument = None,
            etabs_filename : str = None):
        super().__init__(doc=doc, etabs_filename=etabs_filename)

    def perform(self):
        try:
            # start the STKO interface
            self.start()
            # check document
            if self.doc is None:
                raise Exception('No document provided.')
            # parse the input file
            if self.etabs_filename is None:
                raise Exception('No input file provided.')
            the_parser = parser(self.etabs_filename)
            the_parser.doc.process()
            # build the document
            the_builder = builder(the_parser.doc, self)
            # run a regenerate command
            self.regenerate()
        except Exception as e:
            raise e
        finally:
            # stop the STKO interface
            self.stop()