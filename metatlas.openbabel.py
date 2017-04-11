import openbabel as ob
import re
import csv
import subprocess
from fireworks import Firework, LaunchPad, Workflow, FiretaskBase
from fireworks.core.rocket_launcher import launch_rocket
from fireworks.utilities.fw_utilities import explicit_serialize

@explicit_serialize
class CreateOrcaInputTask(FiretaskBase):
    required_params = ['molecule_string']
    optional_params = ['level_of_theory']

    def _create_openbabel_molecule(molecule_string):
        m = ob.OBMol()
        m.SetTitle(molecule_string)
        obc.ReadString(m, molecule_string)
        m.AddHydrogens()
        b.Build(m)

        return m

    def _create_orca_input_string(molecule_string, db):
        molecule = self._create_openbabel_molecule(molecule_string)

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
        db['calculationSetup'] = {
            "molecularSpinMultiplicity": mult,
            "charge": charge,
            "numberOfElectrons": get_n_electrons(m),
            "waveFunctionTheory": "PM3"
        }
        return orca_string

    def _get_n_electrons(molecule):
        elec_count = [atom.GetAtomicNum() 
                for atom in ob.OBMolAtomIter(molecule)]
        return sum(elec_count)

    def run_task(self, fw_spec):
        orca_string = self._create_orca_input_string(self['molecule_string'])
        return FWAction(stored_data={'orca_string': orca_string},
                update_spec={'orca_string': orca_string})

@explicit_serialize
class ComputeEnergyTask(FiretaskBase):
    _fw_name = 'ComputeEnergyTask'

    def _write_string_to_orca_file(orca_string):
        input_name = 'scr/' + m.GetFormula() + '.inp'
        with open(input_name, 'w') as f:
            f.write(orca_string)

        return

    def _create_slurm_file():
        pass

    def _calculate_energy(fname):
        outfile = fname.strip('inp') + 'out'
        out = open(outfile, 'w')
        if np > 1:
            subprocess.call(['mpirun.openmpi', '-np', str(np), \
             '../../orca_3_0_3/orca', fname], stdout=out)
        else:
            subprocess.call(['../../orca_3_0_3/orca', fname], stdout=out)
        return outfile

    def run_task(self, fw_spec):
        try:
            self._calculate_energy(self['molecule'])
        except: # some kind of fault error
            rerun_fw = Firework(ComputeEnergyTask(), 
                    {'molecule': self['molecule']})
            return FWAction(stored_data={'output': output_file}, detour=rerun_fw)
        else:
            return FWAction()

@explicit_serialize
class AddCalculationtoDBTask(FiretaskBase):
    _fw_name = 'AddCalctoDBTask'

    def get_energy(m, f, db):
        output_fname = calculate_energy(f)
        mol = {}
        atoms = []
        egy = 0.0
        converged = False
        start_atoms = False
        legy = False
        with open(output_fname, 'r') as output:
            if 'THE OPTIMIZATION HAS CONVERGED' in output.read():
                print "optimization converged"
                converged = True
            else:
                egy = get_energy(m, f, db)

        with open(output_fname, 'r') as output:
            print output_fname, 'file opened for parsing'
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
                        atoms.append(atom)

                if 'Total Energy' in line:
                    egy = line.split()[3]
                    legy = True
                    print "breaking out!"
                    break

            print output_fname, 'file closed'

        if legy:
            print 'thank god an energy was found -- adding to DB'
            #    db['parentInchi'] = m.GetTitle(mol_string)
            db['childInchi'] = obc.WriteString(m)
            obc.SetOutFormat('smi')
            db['childSmiles'] = obc.WriteString(m)
            obc.SetOutFormat('inchi')
            db['molecularFormula'] = m.GetFormula()
            mol['atoms'] = atoms
            db['molecule'] = mol
            db['totalEnergy'] = {"value": egy, "units": "Hartree"}

        else:
            print "no energy was found for ", m.GetFormula()
            egy = get_energy(m, f, db)

        print 'get_energy is returning properly'
        return egy


    def run_task(self):
        get_energy()
def read_molecules_from_csv(fname):
    mols = {}
    with open(fname) as csvFile:
        csvReader = csv.reader(csvFile)
        for row in csvReader:
            _, inchiString, inchiKey = row[0], row[1], row[2]
            mols[inchiKey] = inchiString
    return mols

def protonate(m, atom, db):
    mm = ob.OBMol(m)
    atom.IncrementImplicitValence()
    mm.AddHydrogens(atom)
    print 'protonating atom', mm.GetFormula()
    totalCharge = m.GetTotalCharge()
    m.SetTotalCharge(totalCharge + 1)
    egy = getEnergy(m, db, mongoDB)
    if egy < protonationEnergy:
        protonatedEnergy = egy
        protonatedAtom = atom.GetIdx()
    atomToDelete = m.getAtom(m.NumAtoms())
    m.DeleteAtom(atomToDelete)
    m.SetTotalCharge(totalCharge)
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
        deprotonatedAtom = connectedAtom.GetIdx()
        print('deprotonation', egy)


def perform_work(mol_string):
    obc.WriteFile(m, 'xyz/' + m.GetFormula() + '.xyz')
    print m.GetFormula(), str(m.NumAtoms()) + ' atoms'

    if 'Mo' in m.GetFormula():
        return

    fname = create_orca_input(m, db)
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

#    Parallel(n_jobs=8)(delayed(deprotonate)(m, atom) \
#            for atom in ob.OBMolAtomIter(m) if atom.IsHydrogen)

#    Parallel(n_jobs=8)(delayed(protonate)(m, atom) \
#            for atom in ob.OBMolAtomIter(m) if not atom.IsHydrogen)

if __name__ == "__main__":
    obc = ob.OBConversion()
    obc.SetInAndOutFormats('inchi', 'xyz')
    obet = ob.OBElementTable()
    b = ob.OBBuilder()

    csv_file = 'metatlas_inchi_inchikey.csv'
    mols = read_molecules_from_csv(csv_file)
    
    lpad = LaunchPad(host='mongodb03.nersc.gov', 
                     port=27017,
                     name='metatlas',
                     username='metatlas_admin',
                     password='2sssj2sssj2a')

    print lpad.get_fw_dict_by_id(4)
    launch_rocket(lpad, fw_id=4)
#    for _, mol_string in mols.iteritems():
#        setup_task = Firework(CreateOrcaInputTask(molecule_string=mol_string))
#        run_calculation = Firework(ComputeEnergyTask())
#        add_info_to_db = Firework(AddCalculationtoDBTask())
#
#        fw = Workflow([setup_task])
##        fw = Workflow([setup_task, run_calculation, add_info_to_db], 
##                {setup_task: [run_calculation], 
##                run_calculation: [add_info_to_db])
#        lpad.add_wf(fw)
#        raise

#    with Parallel(n_jobs=8, verbose=5) as p:
#        p(delayed(perform_work)(mol_string) \
#                for _, mol_string in mols.iteritems())

    #perform_work(mols['UCMIRNVEIXFBKS-UHFFFAOYSA-N'])
