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
import pybel
import openbabel
import numpy as np
from configparser import SafeConfigParser
from rdkit import Chem
from rdkit.Chem import AllChem
from os import remove
from mendeleev import element
from subprocess import Popen, PIPE
from fireworks import Firework, LaunchPad, FiretaskBase, FWAction
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


class RDKitUFFOptimize(FiretaskBase):
    _fw_name = 'RDKitUFFOptimize'
    required_params = ['smiles_string']

    def run_task(self, fw_spec):
        smiles_string = self['smiles_string']
        try:
            mol = Chem.MolFromSmiles(smiles_string)
        except TypeError:
            raise

        try:
            Chem.rdmolops.Kekulize(mol, clearAromaticFlags=True)
        except:
            pass

        mol = Chem.AddHs(mol)

        try:
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
        except RuntimeError:
            mol = None

        return FWAction(
            stored_data={'uff_optimized_mol': mol},
            mod_spec=[{'_push': {'input_mol': mol}}]
        )


class OBUFFOptimize(FiretaskBase):
    _fw_name = 'OBUFFOptimize'
    required_params = ['smiles_string']

    @staticmethod
    def _create_pybel_molecule(smiles_string, lprint=False):
        """create an openbabel molecule from an inchistring"""
        smiles = str(smiles_string)
        try:
            molecule = pybel.readstring('smi', smiles, opt={})
        except TypeError:
            print type(smiles_string)
            print 'Unable to convert smiles string {} to pybel.Molecule'.format(smiles_string)
            raise
        else:
            molecule.title = molecule.formula
            molecule.addh()

            return molecule

    def run_task(self, fw_spec):
        smiles = self['smiles_string']

        obmol = self._create_pybel_molecule(smiles)
        try:
            obmol.make3D()
            obmol.localopt()
        except:
            pass

        xyzfile = '{}-uff.xyz'.format(obmol.formula)
        obmol.write('xyz', filename=xyzfile, overwrite=True)

        orca_string = self.create_orca_input_string(obmol)
        formula = obmol.formula

        return FWAction(
            stored_data={'mol': obmol,
                         'moltype': 'openbabel'},
            mod_spec=[{'_push': {'orca_string': orca_string,
                                 'formula': formula}}]
        )


class Psi4Optimize(FiretaskBase):
    _fw_name = 'Psi4Optimize'
    required_params = ['xyzfile']

    @staticmethod
    def _xyzfile_to_psi4mol(fname):
        qmol = psi4.qcdb.Molecule.init_with_xyz(fname)
        lmol = psi4.geometry(qmol.create_psi4_string_from_molecule())
        lmol.update_geometry()

        return lmol

    @staticmethod
    def _optimize_with_psi4(self, xyzfile):
        psi4mol = xyzfile_to_psi4mol(xyzfile)
        e, wfn = psi4.optimize('pbeh3c/def2-svp', molecule=psi4mol)

        output = wfn.gradient().print_out()

        return output

    def run_task(self, fw_spec):
        xyzfile = self['xyzfile']

        psi4mol = self._xyzfile_to_psi4mol(xyzfile)
        output = self._optimize_with_psi4(psi4mol)

        return FWAction(
            stored_data={'optcoords': output}
        )


class OrcaOptimize(FiretaskBase):
    _fw_name = 'OrcaOptimize'
    optional_params = ['orca_string', 'formula']
    # required_params = ['input_string', 'calc_details']

    def run_task(self, fw_spec):
        orca_string = fw_spec['orca_string'][0][0]
        formula = fw_spec['formula'][0]
        self._write_string_to_orca_file(formula, orca_string)

        try:
            output = optimize_with_orca(formula)
        except ValueError:  # some kind of fault error
            # DON"T KNOW WHAT GOES HERE YET"
            # rerun_fw = Firework(ComputeEnergyTask(input_string=self['input_string'],
            #                                         calc_details=self['calc_details']),
            #                     name=formula)
            # return FWAction(detours=rerun_fw)
            pass

            # Parse results needs a new format to better store more info
        try:
            Results = ParseResults(formula)
        except IOError:
            raise

        return FWAction(
            stored_data={
                'coords': Results.coords,
                'grads': Results.grads,
                'energies': Results.energies,
                'atom_list': Results.atom_list
            })


class ParseResults(object):
    def __init__(self, formula):
        try:
            with open(formula+'.opt') as f:
                contents = f.readlines()
        except IOError:
            raise IOError, 'no output file to parse'
        self.coords, self.grads, self.energies = self._parse_orca_opt_file(contents)

        try:
            with open(formula+'.out') as f:
                contents = f.readlines()
        except IOError:
            raise IOError, 'no output file to parse'

        natoms = len(self.coords[1,:,1])
        self.atom_list = self._get_orca_atom_list(natoms, contents)

    def _parse_orca_opt_file(self, contents):
        steps = {'coords': [], 'gradients': [], 'energies': []}
        coords = []
        grads = []
        energies = []
        get_dims = False
        start = ''
        dims = [1111]

        for i, line in enumerate(contents):
            if 'coordinates' in line:
                start = 'coords'
                get_dims = True
                continue
            elif 'energies' in line:
                start = 'energies'
                get_dims = True
                continue
            elif 'gradients' in line:
                start = 'gradients'
                get_dims = True
                continue

            if get_dims:
                dims = map(int, line.split())
                get_dims = False
                continue

            if len(steps['gradients']) == dims[0]:
                break

            if 'coords' in start:
                coords.extend(map(np.float_, line.split()))
                if len(coords) == dims[1]:
                    steps['coords'].append(coords)
                    coords = []
            elif 'energies' in start:
                energies.extend(map(np.float_, line.split()))
                if len(energies) == dims[0]:
                    steps['energies'].append(energies)
                    energies = []
            elif 'gradients' in start:
                grads.extend(map(np.float_, line.split()))
                if len(grads) == dims[1]:
                    steps['gradients'].append(grads)
                    grads = []

        grads = np.reshape(np.asarray(steps['gradients']), (dims[0], dims[1]/3, 3))
        coords = np.reshape(np.asarray(steps['coords']), (dims[0], dims[1]/3, 3))
        energies = np.reshape(np.asarray(steps['energies']), (dims[0]))

        return coords, grads, energies

    def _get_orca_atom_list(self, natoms, contents):
        start = False
        atom_list = []
        for line in contents:
            if 'CARTESIAN COORDINATES (A.U.)' in line:
                start = True
                continue

            if start:
                match = re.search(r'\ +([0-9]+)\ +([A-Z][a-z]*)\ +([0-9]+\.[0-9]+)\ +([0-9]*)\ +([0-9]+\.[0-9]+)', line)
                if match:
                    num = int(match.group(1))
                    lb = match.group(2)
                    mass = float(match.group(5))
                    atom_list.append((lb, mass))

                    if num+1 == natoms:
                        break
        return atom_list

    def _get_energy(self, contents):
        match = re.search(r'\-[0-9]+\.[0-9]+', contents)
        return match.group(0)


class ProtonateMolecule(OrcaOptimize):
    _fw_name = 'ProtonateMolecule'
    required_params = ['xyzparent']

    def _single_protonations(self, pymol):
        orca_strings = []

        for atom in pymol.atoms:
            if atom.atomicnum in [7,8,15,16]:
                mol = create_obmol(pymol)
                a = mol.NewAtom()
                a.SetAtomicNum(1)
                a.SetVector(atom.coords[0]+0.2,
                atom.coords[1]+0.2,
                atom.coords[2]+0.2)
                mol.AddBond(atom.idx, pymol.atoms[-1].idx+1, 1)
                mol.SetTotalCharge(pymol.charge+1)
                tmp_pymol = pybel.Molecule(mol)
                tmp_pymol.localopt()

                orca_string, _ = create_orca_input_string(tmp_pymol)
                orca_strings.append(orca_string)

        return orca_strings

    def run_task(self, fw_spec):

        print self['xyzparent']
        pymol = pybel.readfile('xyz', self['xyzparent']).next()

        orca_strings = self._single_protonations(pymol)

        molecule_list = []
        for inp in orca_string:
            try:
                write_string_to_orca_file(pymol.formula+'.inp', orca_string)
                output = optimize_with_orca(pymol.formula)
            except ValueError:  # some kind of fault error
                # DON"T KNOW WHAT GOES HERE YET"
                continue

            try:
                Results = ParseResults(pymol.formula)
            except IOError:
                continue

            molecule_list.append(Results)

        return FWAction(
            stored_data={
                'coords': [r.coords for r in molecule_list],
                'grads': [r.grads for r in molecule_list],
                'energies': [r.energies for r in molecule_list],
                'atom_list': [r.atom_list for r in molecule_list]
            })


def make_df_with_smiles_only_from_csv(csv_file, reset=False):
    """
    This function takes in a csv file like the ones Ben has been sending
    around.
    1 creates DataFrame from csv
    2 iterate rows in DataFrame, creating an rdkit.Chem.rdchem.Mol object
        If 'original_smiles' is available, this will be used to generate the
        rdkit.Chem.rdchem.Mol using Chem.MolFromSmiles.
    It will save the entire DataFrame as 'molecules.pkl' and can try to recover
    it as well.
    """
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

    df.to_pickle('molecules.pkl')

    return df


def make_df_with_molecules_from_csv(csv_file, reset=False):
    """
    This function takes in a csv file like the ones Ben has been sending
    around.
    1 creates DataFrame from csv
    2 iterate rows in DataFrame, creating an rdkit.Chem.rdchem.Mol object
        If 'original_smiles' is available, this will be used to generate the
        rdkit.Chem.rdchem.Mol using Chem.MolFromSmiles.
    It will save the entire DataFrame as 'molecules.pkl' and can try to recover
    it as well.
    """
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
        try:
            AllChem.EmbedMolecule(mol)
        except RuntimeError:
            mol = None

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


def get_n_electrons(molecule, rdkit=True):
    if rdkit:
        elec_count = [atom.GetAtomicNum() for atom in molecule.GetAtoms()]
    else:
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
    slurm_adapter = CommonAdapter(
        q_type=q_type,
        template_file='/home/bkrull/.fireworks/slurm.yaml',
        reserve=True)

    return slurm_adapter


def create_orca_input_string(molecule):
    charge = molecule.charge
    nelectrons = get_n_electrons(molecule, rdkit=False)
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
    orca_string += '!COPT\n'
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


def create_obmol(pymol):
    mol = openbabel.OBMol()

    for atom in pymol.atoms:
        a = mol.NewAtom()
        a.SetAtomicNum(atom.atomicnum)
        a.SetVector(atom.coords[0], atom.coords[1], atom.coords[2])

    return mol


def optimize_with_orca(formula):
    iter = 0
    output = ''

    try:
        with open(formula+'.out', 'r') as f:
            output = f.read()
    except IOError:
        while 'OPTIMIZATION RUN DONE' not in output and \
                'TERMINATED NORMALLY' not in output:

            iter += 1
            # process = Popen(['srun', 'orca', formula+'.inp'], stdout=f)

            process = Popen(['orca', formula+'.inp'], stdout=PIPE)
            output, err = process.communicate()
            if iter > 4:
                raise ValueError, "unable to reach opt convergence"

        with open(formula+'.out', 'w') as f:
            f.write(output)

    return output


def write_string_to_orca_file(formula, input_string):
    input_name = formula + '.inp'
    with open(input_name, 'w') as f:
        f.write(input_string)

    return

