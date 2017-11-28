"""
FireWorks implementation for the computation of electronic structure for
molecules in the Metatlas database.

-----------    ----------    -----------
| Create  |    | Run    |    | Process |
| Orca    | => | Orca   | => | Output  |
| Input   |    | Calc   |    | File    |
-----------    ----------    -----------
"""
from configparser import SafeConfigParser
from fireworks import Firework, LaunchPad, Workflow, FWorker
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter
from fireworks.user_objects.metatlas import ComputeEnergyTask, AddCalculationtoDBTask
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from metatlas import read_molecules_from_csv, create_pybel_molecule, create_orca_input_string


def create_launchpad(db_config_file):
    """use to create a FW launchpad using mongodb creds from file"""
    config = SafeConfigParser()
    config.read(db_config_file)
    db = config['db']

    lpad = LaunchPad(
        host=db['host'],
        port=int(db['port']),
        name=db['name'],
        username=db['username'],
        password=db['password'])

    return lpad


def create_fworker(name):
    fworker_config = '/home/bkrull/.fireworks/' + name.lower() + '.yaml'
    fworker = FWorker().from_file(fworker_config)

    return fworker


def create_queue_adapater(q_type):
    slurm_adapter = CommonAdapter(
        q_type=q_type,
        template_file='/home/bkrull/.fireworks/slurm.yaml',
        reserve=True)

    return slurm_adapter


if __name__ == "__main__":
    METATLAS_DB_CONFIG = '/home/bkrull/.fireworks/metatlas.ini'
    CSV_FILE = 'metatlas_inchi_inchikey.csv'
    PROJECT_HOME = 'scr/'

    metatlas_lpad = create_launchpad(METATLAS_DB_CONFIG)
    mols = read_molecules_from_csv(CSV_FILE)

    mol = 'C44H76NO8P'
    fw = Firework(AddCalculationtoDBTask(input_file=mol + '.out'), name=mol)
    wf = Workflow([fw])
    metatlas_lpad.append_wf(wf, [179949])

#    for _, mol_string in mols.iteritems():
#        molecule = create_pybel_molecule(mol_string)
#        orca_string, calc_details = create_orca_input_string(molecule)
#
#        fw1 = Firework(ComputeEnergyTask(input_string=orca_string,
#                                         calc_details=calc_details),
#                       name=molecule.formula)
#        fw2 = Firework(AddCalculationtoDBTask(input_file=molecule.formula+'.out'),
#                       name=molecule.formula)
#
#        wf = Workflow([fw1, fw2], {fw1: [fw2]})
#
#        metatlas_lpad.add_wf(wf)
#        print molecule.formula
#        raise
#perform_work(mols['UCMIRNVEIXFBKS-UHFFFAOYSA-N'])
