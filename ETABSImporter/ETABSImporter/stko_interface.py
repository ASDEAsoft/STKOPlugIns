from PyMpc import *

# The interface to the STKO document
class stko_interface:

    # the message types
    class message_type:
        INFO = 0
        WARNING = 1
        ERROR = 2

    # the constructor
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
    
    # get the definition id of the last definition in the document + 1
    def new_definition_id(self) -> int:
        return self.doc.definitions.getlastkey(0) + 1

    # get the interaction id of the last interaction in the document + 1
    def new_interaction_id(self) -> int:
        return self.doc.interactions.getlastkey(0) + 1

    # get the id of the last physical property in the document + 1
    def new_physical_property_id(self) -> int:
        return self.doc.physicalProperties.getlastkey(0) + 1

    # the the id of the last element property in the document + 1
    def new_element_property_id(self) -> int:
        return self.doc.elementProperties.getlastkey(0) + 1

    # get the condition id of the last condition in the document + 1
    def new_condition_id(self) -> int:
        return self.doc.conditions.getlastkey(0) + 1

    # get the analysis step id of the last analysis step in the document + 1
    def new_analysis_step_id(self) -> int:
        return self.doc.analysisSteps.getlastkey(0) + 1

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

    # adds a new definition to the STKO document
    def add_definition(self, defn : MpcDefinition):
        self.doc.addDefinition(defn)
        self.doc.commitChanges()
        self.doc.dirty = True
        App.processEvents()

    # adds a new physical property to the STKO document
    def add_physical_property(self, prop : MpcProperty):
        # adds a new physical property to the STKO document
        self.doc.addPhysicalProperty(prop)
        self.doc.commitChanges()
        self.doc.dirty = True
        App.processEvents()
    
    # adds a new element property to the STKO document
    def add_element_property(self, prop : MpcElementProperty):
        # adds a new element property to the STKO document
        self.doc.addElementProperty(prop)
        self.doc.commitChanges()
        self.doc.dirty = True
        App.processEvents()

    # adds a new condition to the STKO document
    def add_condition(self, cond : MpcCondition):
        self.doc.addCondition(cond)
        self.doc.commitChanges()
        self.doc.dirty = True
        App.processEvents()

    # adds a new analysis step to the STKO document
    def add_analysis_step(self, step : MpcAnalysisStep):
        self.doc.addAnalysisStep(step)
        self.doc.commitChanges()
        self.doc.dirty = True
        App.processEvents()

    # assign custom colors
    def assign_custom_colors(self):
        # coloring
        self.doc.randomizeMaterialColors()

    # regenerate the document
    def regenerate(self):
        App.runCommand('Regenerate', '2')
        App.processEvents()
    
    # send a percentage to the terminal window
    def send_percentage(self, percentage : float):
        ...

    # send a message to the terminal window
    def send_message(self, msg : str, mtype : message_type = message_type.INFO):
        print(msg)