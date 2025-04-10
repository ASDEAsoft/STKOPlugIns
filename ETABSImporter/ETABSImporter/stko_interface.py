from PyMpc import *

class stko_interface:
    def __init__(
            self,
            doc : MpcCaeDocument):
        # the STKO document
        self.doc = doc
        # todo: add some settings
        
    # start the STKO interface
    def start(self):
        # turn off multithreading on stko
        App.setRegenerateOnWorkingThreadFlag(False)

    # stop the STKO interface
    def stop(self):
        # turn on multithreading on stko
        App.setRegenerateOnWorkingThreadFlag(True)

    # get the geometry id of the last geometry in the document + 1
    def new_geometry_id(self) -> int:
        return self.doc.geometries.getlastkey(0) + 1

    # adds a new geometry to the STKO document
    def add_geometry(self, geom : MpcGeometry):
        self.doc.addGeometry(geom)
        self.doc.commitChanges()
        self.doc.dirty = True
        App.processEvents()
    
    # regenerate the document
    def regenerate(self):
        App.runCommand('Regenerate')
        App.processEvents()