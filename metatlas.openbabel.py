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
import openbabel as ob
from configparser import SafeConfigParser
from fireworks import Firework, LaunchPad, Workflow, FiretaskBase
from fireworks.core.rocket_launcher import launch_rocket
from fireworks.utilities.fw_utilities import explicit_serialize


@explicit_serialize
class CreateOrcaInputTask(FiretaskBase):
    required_params = ['molecule_string']
    optional_params = ['level_of_theory']

    def _create_orca_input_string(self, molecule):
        orca_string = ''
        charge = molecule.GetTotalCharge()
        mult = (1 if (self._get_n_electrons(molecule) + charge) % 2 == 0 else 2)
        orca_string += '%MaxCore 6000\n'
        orca_string += '!SlowConv\n'
        orca_string += '!NOSOSCF\n'
        orca_string += '!PM3 Opt \n%coords \n  CTyp xyz\n'
        orca_string += ' Charge ' + str(charge) + '\n'
        orca_string += ' Mult ' + str(mult) + '\n coords\n'

        for atom in ob.OBMolAtomIter(molecule):
            orca_string += ' ' + obet.GetSymbol(atom.GetAtomicNum())
            orca_string += ' ' + str(atom.GetX())
            orca_string += ' ' + str(atom.GetY())
            orca_string += ' ' + str(atom.GetZ()) + ' \n'

        orca_string += ' end\nend\n'
        orca_string += '%geom\n MaxIter 200\n end\n'
        orca_string += '%scf\n MaxIter 1500\n end\n'
        calc_setup = {
            "molecular_formula": molecule.GetFormula(),
            "molecularSpinMultiplicity": mult,
            "charge": charge,
            "numberOfElectrons": get_n_electrons(m),
            "waveFunctionTheory": "PM3"
        }
        return orca_string, calc_setup

    def _get_n_electrons(self, molecule):
        elec_count = [
            atom.GetAtomicNum() for atom in ob.OBMolAtomIter(molecule)
        ]
        return sum(elec_count)

    def run_task(self, fw_spec):

        molecule = create_openbabel_molecule(self['molecule_string'])
        orca_string, calc_setup = self._create_orca_input_string(molecule)
        run_calculation = Firework(ComputeEnergyTask(),
                                   {'orca_string': orca_string,
                                    'formula': calc_setup['molecular_formula']})
        return FWAction(
            stored_data={
                'orca_string': orca_string,
                'calculation_setup': calc_setup
            },
            update_spec={'orca_string': orca_string},
            additions=run_calculation)


@explicit_serialize
class ComputeEnergyTask(FiretaskBase):
    _fw_name = 'ComputeEnergyTask'

    def _write_string_to_orca_file(self):
        input_name = 'scr/' + self.formula + '.inp'
        with open(input_name, 'w') as f:
            f.write(self.orca_string)

        return

    def _create_slurm_file(self):
        pass

    def _calculate_energy(self):
        path_to_output = self.formula + '.out'
        with open(path_to_output, 'w') as f:
            subprocess.call(['../../orca_3_0_3/orca', fname], stdout=f)

        return path_to_output

    def run_task(self, fw_spec):
        self.orca_string = fw_spec['orca_string']
        self.formula = fw_spec['formula']
        self._write_string_to_orca_file()

        try:
            path_to_output = self._calculate_energy()
        except:  # some kind of fault error
            rerun_fw = Firework(ComputeEnergyTask(),
                                {'molecule': self['molecule']})
            return FWAction(
                stored_data={'path_to_output': path_to_output}, detour=rerun_fw)
        else:
            process_fw = Firework(AddCalculationtoDBTask(),
                                  {'path_to_output': path_to_output})
            return FWAction(
                stored_data={'path_to_output': path_to_output}, addition=process_fw)


@explicit_serialize
class AddCalculationtoDBTask(FiretaskBase):
    required_params = ['input_file']

    def _get_optimized_coords(self, input_file):
        atomic_coords = []
        with open(input_file, 'r') as output:
            for line in output:
                if 'CARTESIAN' in line and 'ANGSTROEM' in line:
                    start_atoms = True

                if start_atoms and 'CARTESIAN' not in line:
                    if '----------' in line:
                        pass
                    elif line == '\n':
                        start_atoms = False
                    else:
                        atom = {}
                        match = re.search(r'\s*(?P<atom>[A-Z][a-z]*)' +
                                          r'\s*(?P<x>\-*[0-9]+\.[0-9]+)' +
                                          r'\s*(?P<y>\-*[0-9]+\.[0-9]+)' +
                                          r'\s*(?P<z>\-*[0-9]+\.[0-9]+)', line)

                        atom['elementSymbol'] = match.group('atom')
                        coords = [match.group(i) for i in ['x', 'y', 'z']]
                        atom['cartesianCoordinates'] = \
                                {'value' : coords, 'units' : 'Angstrom'}
                        atomic_coords.append(atom)

            return atomic_coords

    def _get_energy(self, input_file):
        with open(input_file, 'r') as output:
            print input_file, 'file opened for parsing'
            for line in output:
                if 'Total Energy' in line:
                    egy = line.split()[3]
                    print "breaking out!"
                    break
            print input_file, 'file closed'

        print 'get_energy is returning properly'
        return egy

    def run_task(self, fw_spec):
        energy = self._get_energy(self['input_file'])
        coords = self._get_optimized_coords(self['input_file'])

        return FWAction(
            stored_data={'energy': {
                'value': energy,
                'units': 'Hartree'
            }})


def read_molecules_from_csv(fname):
    mols = {}
    with open(fname) as csvFile:
        csv_reader = csv.reader(csvFile)
        for row in csv_reader:
            _, inchi_string, inchi_key = row[0], row[1], row[2]
            mols[inchi_key] = inchi_string
    return mols


def create_openbabel_molecule(molecule_string):
    m = ob.OBMol()
    m.SetTitle(molecule_string)
    obc.ReadString(m, molecule_string)
    m.AddHydrogens()
    b.Build(m)

    return m


def protonate(m, atom, db):
    mm = ob.OBMol(m)
    atom.IncrementImplicitValence()
    mm.AddHydrogens(atom)
    print 'protonating atom', mm.GetFormula()
    total_charge = m.GetTotalCharge()
    m.SetTotalCharge(total_charge + 1)
    egy = getEnergy(m, db, mongoDB)
    if egy < protonationEnergy:
        protonated_energy = egy
        protonated_atom = atom.GetIdx()
    atom_to_delete = m.getAtom(m.NumAtoms())
    m.DeleteAtom(atom_to_delete)
    m.SetTotalCharge(total_charge)
    print('protonation', egy)


def deprotonate(m, atom, db):
    mm = ob.OBMol(m)
    mm.DeleteAtom(atom)
    print 'deprotonating atom', mm.GetFormula()
    mm.SetTotalCharge(m.GetTotalCharge() - 1)
    try:
        egy = get_energy(mm, db)
    except:
        print "failed to assign deprot energy"

    deprotonated_energy = 0
    if egy < deprotonated_energy:
        deprotonated_energy = egy
    for connected_atom in ob.OBAtomAtomIter(atom):
        deprotonated_atom = connectedAtom.GetIdx()
        print('deprotonation', egy)


def perform_work(mol_string):
    obc.WriteFile(m, 'xyz/' + m.GetFormula() + '.xyz')
    print m.GetFormula(), str(m.NumAtoms()) + ' atoms'

    db['_id'] = fname.split('/')[1].strip('.inp')

    c = neutrals.find({"_id": db['_id']}).count()
    if c > 0:
        print 'database entry ', db['_id'], ' found. moving on.'
    else:
        egy_neutral = get_energy(m, fname, db)

        print 'energy =', egy_neutral

        try:
            insert = neutrals.insert_one(db)
        except pymongo.errors.OperationError:
            print "failed to add energy to database"
        else:
            print "successfully added energy to database"


if __name__ == "__main__":
    config = SafeConfigParser()
    config.read('/home/bkrull/.fireworks/metatlas.ini')
    db = config['db']

    lpad = LaunchPad(
        host=db['host'],
        port=int(db['port']),
        name=db['name'],
        username=db['username'],
        password=db['password'])

    obc = ob.OBConversion()
    obc.SetInAndOutFormats('inchi', 'xyz')
    obet = ob.OBElementTable()
    b = ob.OBBuilder()

    CSV_FILE = 'metatlas_inchi_inchikey.csv'
    MOLS = read_molecules_from_csv(CSV_FILE)

    for _, mol_string in MOLS.iteritems():
        molecule = create_openbabel_molecule(mol_string)
        fw = Firework(CreateOrcaInputTask(molecule_string=mol_string),
                      name=molecule.GetFormula())

        lpad.add_wf(fw)
        id = lpad.get_new_launch_id()
        print lpad.get_fw_dict_by_id(id)
        launch_rocket(lpad, fw_id=id)
        raise

#perform_work(mols['UCMIRNVEIXFBKS-UHFFFAOYSA-N'])
