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
import re
import subprocess
import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from os import remove
from mendeleev import element
from fireworks import Firework, FiretaskBase, FWAction
import psi4

try:
    termtype = get_ipython().__class__.__name__
except:
    from tqdm import tqdm as tqdm
else:
    if 'ZMQ' in termtype:
        from tqdm import tqdm_notebook as tqdm
    else:
        from tqdm import tqdm as tqdm


class ComputeEnergyTask(FiretaskBase):
    _fw_name = 'ComputeEnergyTask'
    required_params = ['input_string', 'calc_details']

    def _write_string_to_orca_file(self, formula, input_string):
        input_name = formula + '.inp'
        with open(input_name, 'w') as f:
            f.write(input_string)

        return

    def _optimize_with_orca(self):
        formula = self['calc_details']['molecular_formula']
        fname = formula + '.inp'
        path_to_output = formula + '.out'
        with open(path_to_output, 'w') as f:
            p = subprocess.Popen(['srun', 'orca', formula+'.inp'], stdout=f)
            p.wait()

        return path_to_output

    def run_task(self, fw_spec):
        formula = self['calc_details']['molecular_formula']
        input_string = self['input_string']
        self._write_string_to_orca_file(formula, input_string)

        try:
            with open(formula+'.out', 'r') as f:
                content = f.read()
        except IOError:
            print 'No output file found yet. Running!'
            try:
                path_to_output = self._optimize_with_orca()
                with open(path_to_output, 'r') as f:
                    contents = f.read()

                if 'TERMINATED NORMALLY' not in contents:
                    raise
            except:  # some kind of fault error
                rerun_fw = Firework(ComputeEnergyTask(input_string=self['input_string'],
                                                      calc_details=self['calc_details']),
                                    name=formula)
                return FWAction(detours=rerun_fw)
        else:
            if 'OPTIMIZATION RUN DONE' not in content:
                print 'Not optimized. Running!'
                try:
                    path_to_output = self._calculate_energy()
                    with open(path_to_output, 'r') as f:
                        contents = f.read()

                    if 'TERMINATED NORMALLY' not in contents:
                        raise
                except:  # some kind of fault error
                    rerun_fw = Firework(ComputeEnergyTask(input_string=self['input_string'],
                                                          calc_details=self['calc_details']),
                                        name=formula)
                    return FWAction(detours=rerun_fw)

        Results = ParseResults(formula+'.trj')

        return FWAction(
            stored_data={
                'energy': {
                    'value': Results.energy,
                    'units': 'Hartree'
                },
                'optimized_coords': Results.opt_coords
            })


class ParseResults(object):
    def __init__(self, path_to_calc_file):
        with open(path_to_calc_file, 'r') as output:
            content = output.readlines()

        natoms = int(content[0].strip())
        # back up 2 additional items, energy, natoms in xyz format
        content = content[-natoms-2:]

        self.opt_coords = ''.join(content)
        self.energy = self._get_energy(content[1])

    def _get_energy(self, contents):
        match = re.search(r'\-[0-9]+\.[0-9]+', contents)
        return match.group(0)


class AddCalculationtoDBTask(FiretaskBase):
    required_params = ['path_to_calc_output']
    def run_task(self, fw_spec):
        path_to_calc_output = ' '
        return FWAction()


class ProtonateMolecule(ComputeEnergyTask):

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

    def run_task(self):
        pass


class DeprotonateMolecule(ComputeEnergyTask):
    def deprotonate(m, atom, db):
        mm = ob.OBMol(m) #     mm.DeleteAtom(atom) #     print 'deprotonating atom', mm.GetFormula()
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

    def run_task(self):
        pass


def make_df_with_molecules_from_csv(csv_file, reset=False):
    if reset:
        try:
            remove('molecules.pkl')
        except:
            pass

    try:
        df = pd.read_pickle('molecules.pkl')
    except IOError:
        print 'No pickle to recover'
    else:
        return df

    df = pd.read_csv(csv_file)
    molecules = []
    kekulized = []
    key_used = []
    for index, row in tqdm(df.iterrows(), total=len(df)):
        mol = None
        for key in ['sanitized_smiles', 'original_smiles']:#, 'sanitized_inchikey']:
            mol_string = row[key]
            try:
                mol = Chem.MolFromSmiles(mol_string)
            except TypeError:
                continue
            else:
                if mol:
                    key_used.append(key)
                    break

        if not mol:
            molecules.append('no')
            key_used.append(None)
            kekulized.append('no')
            continue

        try:
            Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
            kekulized.append('yes')
        except:
            kekulized.append('no')

        mol = Chem.AddHs(mol)
        molecules.append(mol)

    df['molecule'] = molecules
    df['kekulized'] = kekulized
    df['key_used'] = key_used

    df.to_pickle('molecules.pkl')

    return df

def read_molecules_from_csv(fname):
    """ given a csv file, return dict of inchikey to inchistring """
    mols = {}
    with open(fname) as csvFile:
        csv_reader = csv.reader(csvFile)
        for row in csv_reader:
            _, inchi_string, inchi_key = row[0], row[1], row[2]
            mols[inchi_key] = inchi_string
    return mols

def read_molecules_from_csv_new(fname):
    mols = {}
    with open(fname) as csvFile:
        csv_reader = csv.reader(csvFile)
        for row in csv_reader:
            smiles, formula = row[5], row[1]
            mols[formula] = smiles

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

def psi4_xyzfile_to_psi4mol(fname):
    qmol = psi4.qcdb.Molecule.init_with_xyz(fname)
    lmol = psi4.geometry(qmol.create_psi4_string_from_molecule())
    lmol.update_geometry()

    return lmol
