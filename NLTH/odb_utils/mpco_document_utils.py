'''
This module provides utility functions for handling and processing
documents in the MPCO format.
'''

from PyMpc import *
from typing import List, Dict

def get_all_databases() -> Dict[int, str]:
    '''
    Retrieves all available MPCO databases in the current STKO session.

    Returns
    -------
    Dict[int, str]
        A dictionary mapping database IDs to their names.
    '''
    doc = App.postDocument()
    if doc is None:
        raise RuntimeError("No post-processing document is currently open.")
    return {db_id: db.fileName for db_id, db in doc.databases.items()}

def get_all_selection_sets(db : MpcOdb = None) -> Dict[int, str]:
    '''
    Retrieves all selection sets from the given MPCO database.

    Parameters
    ----------
    db : MpcOdb
        The MPCO database object.
        If none, this function will find the active database in the current STKO session.

    Returns
    -------
    Dict[int, str]
        A dictionary mapping selection set IDs to their names.
    '''
    if db is None:
        for _, _db in App.postDocument().databases.items():
            if _db:
                db = _db
                break
    if db is None:
        raise RuntimeError("No MPCO database found in the current post-processing document.")
    return {sset_id: sset.name for sset_id, sset in db.selectionSets.items()}

def extract_mesh(db : MpcOdb) -> MpcMesh:
    '''
    Extracts and returns the mesh from the given MPCO database.
    NOTE: This function returns the mesh at the last step of 
    the last stage.

    Parameters
    ----------
    db : MpcOdb
        The MPCO database object.

    Returns
    -------
    MpcMesh
        The extracted mesh object.
    '''
    if db is None:
        raise ValueError("The provided MPCO database is None.")
    stages = db.getStageIDs()
    if len(stages) == 0:
        raise ValueError("No stages found in the MPCO database.")
    last_stage_id = stages[-1]
    steps = db.getStepIDs(last_stage_id)
    if len(steps) == 0:
        raise ValueError(f"No steps found in stage '{last_stage_id}'.")
    last_step_id = steps[-1]
    U = db.getNodalResult('Displacement')
    opt = MpcOdbVirtualResultEvaluationOptions()
    opt.stage = last_stage_id
    opt.step = last_step_id
    field = U.evaluate(opt)
    mesh = field.mesh
    return mesh

def extract_isolator_nodes(db : MpcOdb, mesh : MpcMesh, sset_id : int) -> List[int]:
    '''
    Extracts and returns a list of nodes associated with isolators
    from the specified set in the MPCO database.

    Parameters
    ----------
    db : MpcOdb
        The MPCO database object.
    mesh : MpcMesh
        The mesh object containing nodes and elements.
    sset_id : int
        The ID of the set containing isolator elements.

    Returns
    -------
    Tuple[List[int], List[int]]
        2 lists of nodes IDs associated with isolators.
        The first list contains the bottom nodes (for reactions),
        and the second list contains the top nodes (for displacements).
    '''
    bottom_nodes : List[int] = []
    top_nodes : List[int] = []
    try:
        sset = db.selectionSets[sset_id]
    except KeyError:
        raise ValueError(f"Selection Set '{sset_id}' not found in the MPCO database.")
    for element_id in sset.elements:
        element = mesh.getElement(element_id)
        if element is None:
            raise ValueError(f"Element with ID '{element_id}' not found in the mesh.")
        if len(element.nodes) != 2:
            raise ValueError(f"Element with ID '{element_id}' does not have exactly 2 nodes.")
        N1, N2 = element.nodes
        if N1.position.z < N2.position.z:
            bottom_nodes.append(N1.id)
            top_nodes.append(N2.id)
        else:
            bottom_nodes.append(N2.id)
            top_nodes.append(N1.id)
    return bottom_nodes, top_nodes