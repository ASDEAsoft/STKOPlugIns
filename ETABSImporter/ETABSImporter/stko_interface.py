from PyMpc import *

# The interface to the STKO document
class stko_interface:
    def __init__(
            self,
            doc : MpcCaeDocument = None,
            etabs_filename : str = None):
        self.doc = doc
        self.etabs_filename = etabs_filename
        
    # the main method
    def perform(self):
        raise Exception('perform() method not implemented.')
 
    # start the STKO interface
    def start(self):
        # turn off multithreading on stko
        App.setRegenerateOnWorkingThreadFlag(False)

    # stop the STKO interface
    def stop(self):
        # turn on multithreading on stko
        App.setRegenerateOnWorkingThreadFlag(True)

    # get the local axes id of the last local axes in the document + 1
    def new_local_axes_id(self) -> int:
        return self.doc.localAxes.getlastkey(0) + 1

    # get the geometry id of the last geometry in the document + 1
    def new_geometry_id(self) -> int:
        return self.doc.geometries.getlastkey(0) + 1

    # get the interaction id of the last interaction in the document + 1
    def new_interaction_id(self) -> int:
        return self.doc.interactions.getlastkey(0) + 1

    # adds a new local axes to the STKO document
    def add_local_axes(self, locax : MpcLocalAxes):
        # adds a new local axes to the STKO document
        self.doc.addLocalAxes(locax)
        self.doc.commitChanges()
        self.doc.dirty = True
        App.processEvents()

    # adds a new geometry to the STKO document
    def add_geometry(self, geom : MpcGeometry):
        self.doc.addGeometry(geom)
        self.doc.commitChanges()
        self.doc.dirty = True
        App.processEvents()
    
    # adds a new interaction to the STKO document
    def add_interaction(self, inter : MpcInteraction):
        self.doc.addInteraction(inter)
        self.doc.commitChanges()
        self.doc.dirty = True
        App.processEvents()

    # regenerate the document
    def regenerate(self):
        App.runCommand('Regenerate', '2')
        App.processEvents()
    
    # send a message to the terminal window
    def send_message(self, msg : str):
        print(msg)