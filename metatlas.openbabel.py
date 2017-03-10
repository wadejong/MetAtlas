#
# Base code for calculating most likely protonation sites
#
# Using OpenBabel as the toolkit to generate protonated and deprotonated structures
#

import openbabel as ob
import re
import csv
import subprocess 
import pymongo
import bson
from joblib import Parallel, delayed

def get_n_electrons(m):
  electron_count = 0
  for a in ob.OBMolAtomIter(m):
    electron_count += a.GetAtomicNum()
  return(electron_count)

def create_orca_input(m, db):
    charge = m.GetTotalCharge()
    mult = (1 if (get_n_electrons(m)+charge) % 2 == 0 else 2)
    input_name = 'scr/'+m.GetFormula()+'.inp'
    inp = open(input_name,'w')
    inp.write('!PM3 Opt \n%coords \n  CTyp xyz\n')
    inp.write(' Charge '+str(charge)+'\n Mult '+str(mult)+' \n coords\n')
    for atom in ob.OBMolAtomIter(m):
        inp.write('  '+obet.GetSymbol(atom.GetAtomicNum())+' '
                +str(atom.GetX())+' '+str(atom.GetY())+' '+str(atom.GetZ())+' \n')
    inp.write(' end\nend\n')
    inp.write('%geom\n MaxIter 200\n end\n')
    inp.write('%scf\n MaxIter 1500\n end\n')
    inp.close()
    db['calculationSetup'] = { "molecularSpinMultiplicity" : mult, 
            "charge" : charge, 
            "numberOfElectrons" : get_n_electrons(m), 
            "waveFunctionTheory" : "PM3" } 
    subprocess.call(['cp',input_name,'inp'+input_name.strip('scr')])
    return(input_name) 

def calculate_energy(fname, np):
    outfile = fname.strip('inp')+'out'
    out = open(outfile,'w')
    if np > 1:
        subprocess.call(['mpirun.openmpi', '-np', str(np), \
         '../../orca_3_0_3/orca',fname],stdout=myOutput)
    else:
        subprocess.call(['../../orca_3_0_3/orca',fname],stdout=out)
    return(outfile) 

def get_energy(m, f, db, mongoDB):
    db['_id'] = f.strip('.inp')
    c = mongoDB.find({db['_id']: {"$exists":"true"}}).limit(1)
    if c.count() > 0: 
        print 'database entry ', db['_id'], ' found. moving on.'
        egy = 0.0
    else:
        output_fname = calculate_energy(f, np)
        mol = {}
        atoms = []
        egy = 0.0
        converged = False
        start_atoms = False
        legy = False
        with open(output_fname,'r') as output:
            if 'THE OPTIMIZATION HAS CONVERGED' in output.read(): 
                print "optimization converged"
                converged = True
            else:
                egy = get_energy(m, f, db, mongoDB)

        with open(output_fname,'r') as output:
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
                        match = re.search(r'\s*(?P<atom>[A-Z][a-z]*)'+
                                        '\s*(?P<x>\-*[0-9]+\.[0-9]+)'+
                                        '\s*(?P<y>\-*[0-9]+\.[0-9]+)'+
                                        '\s*(?P<z>\-*[0-9]+\.[0-9]+)', line)

                        atom['elementSymbol'] = match.group('atom')
                        coords = [ match.group(i) for i in ['x', 'y', 'z'] ]
                        atom['cartesianCoordinates'] = \
                                { 'value' : coords, 'units' : 'Angstrom' }
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
            db['totalEnergy'] = {"value":egy, "units":"Hartree"}

            try:
                insert = mongoDB.insert_one(db)
            except pymongo.errors.OperationError:
                print "failed to add energy to database"
            else:
                print "successfully added energy to database"
        else:
            print "no energy was found for ", m.GetFormula()
            egy = get_energy(m, f, db, mongoDB)

    print 'get_energy is returning properly'
    return(egy)

def read_molecules_from_csv(fname):
    mols = {}
    with open(fname) as csvFile:
        csvReader = csv.reader(csvFile)
        for row in csvReader: 
            _,inchiString,inchiKey = row[0],row[1],row[2]
            mols[inchiKey] = inchiString
    return mols

def protonate(m, atom, db):
    mm = ob.OBMol(m)
    atom.IncrementImplicitValence()
    mm.AddHydrogens(atom)
    print 'protonating atom', mm.GetFormula()
    totalCharge = m.GetTotalCharge()
    m.SetTotalCharge(totalCharge+1)
    egy = getEnergy(m, db, mongoDB)
    if egy < protonationEnergy:
        protonatedEnergy = egy
        protonatedAtom = atom.GetIdx()
    atomToDelete = m.getAtom(m.NumAtoms())
    m.DeleteAtom(atomToDelete)
    m.SetTotalCharge(totalCharge)
    print('protonation',egy)

def deprotonate(m, atom, db):
    mm = ob.OBMol(m)
    mm.DeleteAtom(atom)
    print 'deprotonating atom', mm.GetFormula()
    mm.SetTotalCharge(m.GetTotalCharge()-1)
    try:
        egy = get_energy(mm, db, mongoDB)
    except:
        print "failed to assign deprot energy"
        
    if egy < deprotonatedEnergy:
        deprotonatedEnergy = egy
    for connectedAtom in ob.OBAtomAtomIter(atom):
        deprotonatedAtom = connectedAtom.GetIdx()
        print('deprotonation',egy)

def perform_work(mol_string):
    db = {} 
    mongoClient = pymongo.MongoClient()
    mongoDB = mongoClient.test.metatlas

    m = ob.OBMol()
    m.SetTitle(mol_string)
    obc.ReadString(m,mol_string)
    m.AddHydrogens()
    b.Build(m)
    obc.WriteFile(m, 'xyz/'+m.GetFormula()+'.xyz')
    print m.GetFormula(), str(m.NumAtoms())+' atoms'

    fname = create_orca_input(m, db)
    egy_neutral = get_energy(m, fname, db, mongoDB)

    print 'energy =', egy_neutral

#    Parallel(n_jobs=8)(delayed(deprotonate)(m, atom) \
#            for atom in ob.OBMolAtomIter(m) if atom.IsHydrogen)

#    Parallel(n_jobs=8)(delayed(protonate)(m, atom) \
#            for atom in ob.OBMolAtomIter(m) if not atom.IsHydrogen)

if __name__ == "__main__":
    obc = ob.OBConversion()
    obc.SetInAndOutFormats('inchi','xyz')
    obet = ob.OBElementTable()
    b = ob.OBBuilder()

    csv_file = 'metatlas_inchi_inchikey.csv' 
    mols = read_molecules_from_csv(csv_file)

    np =1
    with Parallel(n_jobs=8, verbose=5) as p:
        p(delayed(perform_work)(mol_string) \
                for _, mol_string in mols.iteritems())

    #perform_work(mols['UCMIRNVEIXFBKS-UHFFFAOYSA-N'])
