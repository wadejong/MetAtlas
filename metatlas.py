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
import pybel
from configparser import SafeConfigParser
from fireworks import Firework, LaunchPad, Workflow, FiretaskBase
from fireworks.core.rocket_launcher import launch_rocket
from fireworks.utilities.fw_utilities import explicit_serialize


@explicit_serialize
class CreateOrcaInputTask(FiretaskBase):
    required_params = ['molecule_string']
    optional_params = ['level_of_theory']

    def _create_orca_input_string(self, molecule):
        obet = ob.OBElementTable()

        charge = molecule.charge
        nelectrons = self._get_n_electrons(molecule)
        mult = (1 if (nelectrons + charge) % 2 == 0 else 2)
        calc_setup = {
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
            tmp += ' ' + obet.GetSymbol(atom.atomicnum)
            tmp += ' ' + str(atom.coords[0])
            tmp += ' ' + str(atom.coords[1])
            tmp += ' ' + str(atom.coords[2]) + ' \n'
            print tmp
            orca_string += tmp

        orca_string += ' end\nend\n'
        orca_string += '%geom\n MaxIter 200\n end\n'
        orca_string += '%scf\n MaxIter 1500\n end\n'

        return orca_string, calc_setup

    def _get_n_electrons(self, molecule):
        elec_count = [atom.atomicnum for atom in molecule.atoms]
        return sum(elec_count)

    def run_task(self, fw_spec):
        molecule = create_pybel_molecule(self['molecule_string'], lprint=True)
        print 'type of molecule in run task is ', type(molecule)
        orca_string, calc_setup = self._create_orca_input_string(molecule)
        run_calculation = Firework(ComputeEnergyTask(),
                                   {'orca_string': orca_string,
                                    'formula': molecule.formula})
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


# def protonate(m, atom, db):
#     mm = ob.OBMol(m)
#     atom.IncrementImplicitValence()
#     mm.AddHydrogens(atom)
#     print 'protonating atom', mm.GetFormula()
#     total_charge = m.GetTotalCharge()
#     m.SetTotalCharge(total_charge + 1)
#     egy = getEnergy(m, db, mongoDB)
#     if egy < protonationEnergy:
#         protonated_energy = egy
#         protonated_atom = atom.GetIdx()
#     atom_to_delete = m.getAtom(m.NumAtoms())
#     m.DeleteAtom(atom_to_delete)
#     m.SetTotalCharge(total_charge)
#     print('protonation', egy)


# def deprotonate(m, atom, db):
#     mm = ob.OBMol(m)
#     mm.DeleteAtom(atom)
#     print 'deprotonating atom', mm.GetFormula()
#     mm.SetTotalCharge(m.GetTotalCharge() - 1)
#     try:
#         egy = get_energy(mm, db)
#     except:
#         print "failed to assign deprot energy"

#     deprotonated_energy = 0
#     if egy < deprotonated_energy:
#         deprotonated_energy = egy
#     for connected_atom in ob.OBAtomAtomIter(atom):
#         deprotonated_atom = connectedAtom.GetIdx()
#         print('deprotonation', egy)


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


def create_launchpad():
    """use to create a FW launchpad using mongodb creds from file"""
    config = SafeConfigParser()
    config.read('/home/bkrull/.fireworks/metatlas.ini')
    db = config['db']

    lpad = LaunchPad(
        host=db['host'],
        port=int(db['port']),
        name=db['name'],
        username=db['username'],
        password=db['password'])

    return lpad


if __name__ == "__main__":
    csv_file = 'metatlas_inchi_inchikey.csv'
    mols = read_molecules_from_csv(csv_file)

    lpad = create_launchpad()

    lpad.reset('2017-04-12')
    for _, mol_string in mols.iteritems():
        molecule = create_pybel_molecule(mol_string)
        fw = Firework(CreateOrcaInputTask(molecule_string=mol_string),
                      name=molecule.formula)

        lpad.add_wf(fw)
        id = lpad.get_new_launch_id()
        print lpad.get_fw_dict_by_id(id)
        launch_rocket(lpad, fw_id=id)

#perform_work(mols['UCMIRNVEIXFBKS-UHFFFAOYSA-N'])