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
from glob import glob
from fireworks import Firework
from metatlas import ProtonateMolecule, create_launchpad


if __name__ == "__main__":
    METATLAS_DB_CONFIG = '/home/bkrull/.fireworks/local_db.ini'
    QM9_DATA_HOME = '/home/bkrull/Documents/data/qm9'
    PROJECT_HOME = 'scr/'

    lpad = create_launchpad(METATLAS_DB_CONFIG)
    lpad.reset('2018-03-14')

    i=0
    for fname in glob(QM9_DATA_HOME+'/*'):
        pm3 = Firework(ProtonateMolecule(xyzparent=fname),
                       name=fname.split('.')[0])

        lpad.add_wf(pm3)
        i+=1

        if i>20:break
