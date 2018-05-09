"""
FireWorks implementation for the computation of electronic structure for
molecules in the Metatlas database.

-----------    ----------    -----------
| Create  |    | Run    |    | Process |
| Orca    | => | Orca   | => | Output  |
| Input   |    | Calc   |    | File    |
-----------    ----------    -----------
"""
import pybel
import sys
import numpy as np
from multiprocessing import Pool
from glob import glob
from fireworks import Firework, FWorker
from fireworks.core.rocket_launcher import rapidfire
from fireworks.core.launchpad import LockedWorkflowError
from metatlas import ProtonateMolecule, create_launchpad


def multirapidfire(nthreads):
    pool = Pool(processes=nthreads)
    pool.map(rapid, range(nthreads))

def rapid(dummy):
    lpad = create_launchpad(LOCAL_DB_CONFIG)
    rapidfire(lpad, FWorker(), nlaunches=25000)

def add_fws(reset=False):
    if reset:
        lpad.reset('', require_password=False)
    lpad = create_launchpad(LOCAL_DB_CONFIG)

    for fname in glob(QM9_DATA_HOME+'/*'):
        pm3 = Firework(ProtonateMolecule(xyzparent=fname),
                       name=fname.split('.')[0])

        lpad.add_wf(pm3)
    return

def add_neutral_fws(reset=False):
    if reset:
        lpad.reset('', require_password=False)
    lpad = create_launchpad(LOCAL_DB_CONFIG)

    for fname in glob(QM9_DATA_HOME+'/*'):
        smiles = row['original_smiles']
        formula = row['formula']

        if type(smiles) is float:
            continue

        uff = OBUFFOptimize(smiles_string=smiles)
        pm3 = OrcaOptimize()
        fw = Firework([uff, pm3], name=formula)

        lpad.add_wf(fw)
    return

def multi_update(nthreads):
    batch_ids = map(int, np.linspace(0, len(completed_ids), nthreads))

    batches = [completed_ids[batch_ids[i]:batch_ids[i+1]]
               for i, _ in enumerate(batch_ids[:-1])]

    pool = Pool(processes=nthreads)
    pool.map(rapid_update, batches)

def rapid_update(batch):
    lpad = create_launchpad(LOCAL_DB_CONFIG)
    for fw_id in batch:
        try:
            lpad.mark_fizzled(fw_id)
            lpad.rerun_fw(fw_id, recover_mode='prev_dir')

            fwdict = lpad.get_fw_dict_by_id(fw_id)
            parent_xyz = str(fwdict['name']+'.xyz')
            mol = pybel.readfile('xyz', parent_xyz).next()

            lpad.update_spec([fw_id], {'_tasks.0._fw_name': 'ProtonateMolecule',
                                    '_tasks.0.xyzparent': parent_xyz})
        except LockedWorkflowError:
            continue

if __name__ == "__main__":
    LOCAL_DB_CONFIG = '/home/bkrull/.fireworks/qm9_local.ini'
    QM9_DATA_HOME = '/home/bkrull/Documents/data/qm9'
    PROJECT_HOME = 'scr/'

    # lpad = create_launchpad(LOCAL_DB_CONFIG)
    # completed_ids = lpad.get_fw_ids({'state': 'FIZZLED'})

    nthreads = int(sys.argv[1])
    # multi_update(nthreads)
    multirapidfire(nthreads)
