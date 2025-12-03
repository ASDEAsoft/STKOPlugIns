from PyMpc import *
from NLTH.odb_utils.mpco_document_utils import extract_mesh, extract_isolator_nodes
import numpy as np
from typing import List, Dict, Tuple

class _progress_data:
    def __init__(self):
        self.progress = 0.0
        self.current_db_start = 0.0
        self.current_db_end = 0.0
        self.current_time_increment = 0.0
    def increment(self):
        db_span = self.current_db_end - self.current_db_start
        delta = db_span * self.current_time_increment
        self.progress += delta

class NLTHEvaluatorCallbackBase:
    def __init__(self):
        self._db_id = 0
        self._index = 0
        self._count = 0
        self._update_requested = False
        self._update_period = 0.1
        self._last_update_percentage = 0.0

    def on_progress(self, progress : float):
        """Called to report progress.
        
        Parameters
        ----------
        progress : float
            Progress value between 0.0 and 1.0
        """
        print(f"Progress: {progress*100:.2f}%")
        App.processEvents()

    def set_update_requested(self):
        """Called to request an update of the range.
        """
        self._update_requested = True

    def update_range_requested(self) -> bool:
        """Called when the evaluator requests an update of the range.
        """
        return self._update_requested
    
    def set_update_range(self, db_id : int, index : int, count : int):
        """Called to set the range of data to be updated.
        
        Parameters
        ----------
        db_id : int
            Database ID
        index : int
            Starting index
        count : int
            Number of items to update
        """
        self._db_id = db_id
        self._index = index
        self._count = count
        self._update_requested = False

    def get_update_range(self) -> Tuple[int, int, int]:
        """Get the current update range.
        
        Returns
        -------
        Tuple[int, int, int]
            (db_id, index, count)
        """
        return (self._db_id, self._index, self._count)

class NLTHEvaluator:

    def __init__(self,
            db_ids : List[int],
            isolator_sset_id : int,
            callback : NLTHEvaluatorCallbackBase = None):
        # store parameters
        self.db_ids = [i for i in db_ids]
        print(f"NLTHEvaluator initialized with DB IDs: {self.db_ids}")
        self.isolator_sset_id = isolator_sset_id
        # store data for progress
        self.progress = _progress_data()
        # callback, make sure the type inherits from base
        self.callback = callback
        if self.callback is None:
            self.callback = NLTHEvaluatorCallbackBase()
        if not isinstance(self.callback, NLTHEvaluatorCallbackBase):
            raise TypeError("Callback must be an instance of NLTHEvaluatorCallbackBase")
        self.callback.on_progress(0.0)
        self.callback.set_update_range(0, 0, 0)
        # output data
        self.db_steps : Dict[int, List[int]] = {} # {db_id: [time_steps]}
        self.db_times : Dict[int, np.ndarray] = {} # {db_id: np.array([time_values])}
        self.isolator_displacements : Dict[int, Dict[int, np.ndarray]] = {}  # {db_id: {node_id: (ux, uy, uz)}}
        self.isolator_reactions : Dict[int, Dict[int, np.ndarray]] = {}      # {db_id: {node_id: (rx, ry, rz)}}
    
    def _get_time_info(self, db : MpcOdb):
        stages = db.getStageIDs()
        if len(stages) == 0:
            raise RuntimeError("No stages found in the MPCO database.")
        last_stage_id = stages[-1]
        steps = db.getStepIDs(last_stage_id)
        times = db.getStepTimes(last_stage_id)
        if len(steps) == 0:
            raise RuntimeError(f"No steps found in stage '{last_stage_id}'.")
        return last_stage_id, steps, times

    def _get_reduced_meshes(self, db : MpcOdb) -> MpcMesh:
        # get the full mesh
        mesh = extract_mesh(db)
        # the isolator meshes
        bottom_nodes, top_nodes = extract_isolator_nodes(db, mesh, self.isolator_sset_id)
        iso_mesh_bot = MpcMesh()
        for node_id in bottom_nodes:
            iso_mesh_bot.addNode(mesh.getNode(node_id))
        iso_mesh_top = MpcMesh()
        for node_id in top_nodes:
            iso_mesh_top.addNode(mesh.getNode(node_id))
        # done
        return iso_mesh_bot, iso_mesh_top

    def _evaluate_db(self, db : MpcOdb):
        # get time info
        stage_id, steps, times = self._get_time_info(db)
        Nsteps = len(steps)
        self.progress.current_time_increment = 1.0 / Nsteps # same for all time steps
        # save time info in output for this DB
        OUT_STEPS = [i for i in steps]
        OUT_TIMES = np.array(times)
        self.db_steps[db.id] = OUT_STEPS
        self.db_times[db.id] = OUT_TIMES
        # get reduced meshes
        iso_mesh_bot, iso_mesh_top = self._get_reduced_meshes(db)
        # get results
        U = db.getNodalResult('Displacement')
        R = db.getNodalResult('Reaction Force')
        # allocate output data structures for each node
        OUT_U_DICT : Dict[int, np.ndarray] = {node_id: np.zeros((Nsteps, 3)) for node_id in iso_mesh_top.nodes.keys()}
        OUT_R_DICT : Dict[int, np.ndarray] = {node_id: np.zeros((Nsteps, 3)) for node_id in iso_mesh_bot.nodes.keys()}
        self.isolator_displacements[db.id] = OUT_U_DICT
        self.isolator_reactions[db.id] = OUT_R_DICT
        for current_step, step_id in enumerate(steps):
            # evaluate displacements at top isolator nodes
            opt = MpcOdbVirtualResultEvaluationOptions()
            opt.stage = stage_id
            opt.step = step_id
            # displacements on top isolator nodes
            opt.mesh = iso_mesh_top
            U_field = U.evaluate(opt)
            # reactions on bottom isolator nodes
            opt.mesh = iso_mesh_bot
            R_field = R.evaluate(opt)
            # process results as needed...
            for node_id in iso_mesh_top.nodes.keys():
                row = MpcOdbResultField.node(node_id)
                value = U_field[row]
                NODE_U = OUT_U_DICT[node_id]
                NODE_U[current_step, 0] = value[0]
                NODE_U[current_step, 1] = value[1]
                NODE_U[current_step, 2] = value[2]
            for node_id in iso_mesh_bot.nodes.keys():
                row = MpcOdbResultField.node(node_id)
                value = R_field[row]
                NODE_R = OUT_R_DICT[node_id]
                NODE_R[current_step, 0] = value[0]
                NODE_R[current_step, 1] = value[1]
                NODE_R[current_step, 2] = value[2]
            # update progress
            self.progress.increment()
            self.callback.on_progress(self.progress.progress)
            # check if update range requested
            if self.callback.update_range_requested() or (current_step == Nsteps - 1):
                last_db_id, last_index, last_count = self.callback.get_update_range()
                if last_db_id != db.id:
                    last_index = 0
                    last_count = 0
                new_index = min(last_index + last_count, Nsteps - 1)
                new_count = current_step - new_index + 1
                self.callback.set_update_range(db.id, new_index, new_count)
    
    def evaluate(self):
        doc = App.postDocument()
        if doc is None:
            raise RuntimeError("No post-processing document is currently open.")
        NDB = len(self.db_ids)
        if NDB == 0:
            raise RuntimeError("No database IDs provided for evaluation.")
        db_span = 1.0 / NDB # same for all DBs
        for i, db_id in enumerate(self.db_ids):
            db = doc.databases.get(db_id, None)
            self.progress.current_db_start = i * db_span
            self.progress.current_db_end = (i+1) * db_span
            if db is None:
                print(f"Database ID {db_id} not found.")
                continue
            self._evaluate_db(db)
        # just to make sure we reach 100%
        self.progress.progress = 1.0