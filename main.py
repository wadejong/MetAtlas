"""
FireWorks implementation for the computation of electronic structure for
molecules in the Metatlas database.

-----------    ----------    -----------
| Create  |    | Run    |    | Process |
| Orca    | => | Orca   | => | Output  |
| Input   |    | Calc   |    | File    |
-----------    ----------    -----------
"""
import re
import csv
import subprocess
import pybel
from mendeleev import element
from configparser import SafeConfigParser
from fireworks import Firework, LaunchPad, Workflow, FiretaskBase, FWAction, FWorker
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter
from fireworks.user_objects.metatlas import ComputeEnergyTask
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.utilities.fw_utilities import explicit_serialize


def read_molecules_from_csv(fname):
    """ given a csv file, return dict of inchikey to inchistring """
    mols = {}
    with open(fname) as csvFile:
        csv_reader = csv.reader(csvFile)
        for row in csv_reader:
            _, inchi_string, inchi_key = row[0], row[1], row[2]
            mols[inchi_key] = inchi_string
    return mols


def create_pybel_molecule(inchi_string, lprint=False):
    """create an openbabel molecule from an inchistring"""
    if lprint:
        print 'inchi string in create mol is ', inchi_string

    try:
        molecule = pybel.readstring('inchi', inchi_string, opt={})
    except TypeError:
        print 'Unable to convert inchi string to pybel.Molecule'
        quit()
    else:
        molecule.title = molecule.formula
        molecule.addh()
        molecule.make3D()

        return molecule


def create_orca_input_string(molecule):
    charge = molecule.charge
    nelectrons = get_n_electrons(molecule)
    mult = (1 if (nelectrons + charge) % 2 == 0 else 2)
    calc_details = {
        "molecular_formula": molecule.formula,
        "molecularSpinMultiplicity": mult,
        "charge": charge,
        "numberOfElectrons": nelectrons,
        "waveFunctionTheory": "PM3"
    }

    orca_string = ''
    orca_string += '%MaxCore 6000\n'
    orca_string += '!SlowConv\n'
    orca_string += '!NOSOSCF\n'
    orca_string += '!PM3 Opt \n%coords \n  CTyp xyz\n'
    orca_string += ' Charge ' + str(charge) + '\n'
    orca_string += ' Mult ' + str(mult) + '\n coords\n'

    for atom in molecule.atoms:
        tmp = ''
        tmp += ' ' + str(element(atom.atomicnum).symbol)
        tmp += ' ' + str(atom.coords[0])
        tmp += ' ' + str(atom.coords[1])
        tmp += ' ' + str(atom.coords[2]) + ' \n'
        orca_string += tmp

    orca_string += ' end\nend\n'
    orca_string += '%geom\n MaxIter 200\n end\n'
    orca_string += '%scf\n MaxIter 1500\n end\n'

    return orca_string, calc_details


def get_n_electrons(molecule):
    elec_count = [atom.atomicnum for atom in molecule.atoms]
    return sum(elec_count)


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
    slurm_adapter = CommonAdapter(q_type=q_type,
                                  template_file='/home/bkrull/.fireworks/slurm.yaml',
                                  reserve=True)

    return slurm_adapter


if __name__ == "__main__":
    METATLAS_DB_CONFIG = '/home/bkrull/.fireworks/metatlas.ini'
    CSV_FILE = 'metatlas_inchi_inchikey.csv'
    PROJECT_HOME = 'scr/'

    metatlas_lpad = create_launchpad(METATLAS_DB_CONFIG)
    edison = create_fworker(name='Edison')
    slurm_adapter = create_queue_adapater(q_type='SLURM')

    metatlas_lpad.reset('2017-04-20')

    mols = read_molecules_from_csv(CSV_FILE)

    count = 0
    for _, mol_string in mols.iteritems():
        molecule = create_pybel_molecule(mol_string)
        orca_string, calc_details = create_orca_input_string(molecule)

        fw = Firework(ComputeEnergyTask(input_string=orca_string,
                                        calc_details=calc_details),
                      name=molecule.formula)

        metatlas_lpad.add_wf(fw)

        if count > 100:
            break
        else:
            count += 1
        # id = metatlas_lpad.get_new_launch_id()
        # launch_rocket_to_queue(metatlas_lpad,
        #                        edison,
        #                        slurm_adapter,
        #                        launcher_dir=PROJECT_HOME,
        #                        create_launcher_dir=True,
        #                        fw_id=id)

#perform_work(mols['UCMIRNVEIXFBKS-UHFFFAOYSA-N'])
