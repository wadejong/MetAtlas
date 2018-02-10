"""
FireWorks implementation for the computation of electronic structure for
molecules in the Metatlas database.

-----------    ----------    -----------
| Create  |    | Run    |    | Process |
| Orca    | => | Orca   | => | Output  |
| Input   |    | Calc   |    | File    |
-----------    ----------    -----------
"""
import pandas as pd
from configparser import SafeConfigParser
from fireworks import Firework, LaunchPad, Workflow, FWorker
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter
from fireworks.user_objects.metatlas import ComputeEnergyTask, AddCalculationtoDBTask
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from metatlas import create_orca_input_string, make_df_with_molecules_from_csv
from tqdm import tqdm
from subprocess import check_output


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
    CSV_FILE = './chebi_and_metacyc_molecules-short.csv'
    PROJECT_HOME = 'scr/'

    metatlas_lpad = create_launchpad(METATLAS_DB_CONFIG)
    molecules = make_df_with_molecules_from_csv(CSV_FILE, reset=True)

    for index, row in tqdm(molecules.iterrows()):
        orca_string, calc_details = create_orca_input_string(row['molecule'])

        fw = Firework(ComputeEnergyTask(input_string=orca_string,
                                        calc_details=calc_details)
        metatlas_lpad.add_wf(fw)
