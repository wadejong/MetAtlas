#
# Base code for calculating most likely protonation sites
#
# Using OpenBabel as the toolkit to generate protonated and deprotonated structures
#

import openbabel
import uuid
import csv
import subprocess 
import pymongo
import bson
from joblib import Parallel, delayed
  
def mySplit(s, nvars):
  return (s.split() + [None] * nvars)[:nvars] 

def get_n_electrons(m):
  electron_count = 0
  for a in openbabel.OBMolAtomIter(m):
    electron_count += a.GetAtomicNum()
  return(electron_count)

def create_orca_input(m):
    charge = m.GetTotalCharge()
    mult = (1 if (get_n_electrons(m)+charge) % 2 == 0 else 2)
    input_name = 'scr/'+str(uuid.uuid4())+'.inp'
    inp = open(input_name,'w')
    inp.write('!PM3 Opt \n%coords \n  CTyp xyz\n')
    inp.write(' Charge '+str(charge)+'\n Mult '+str(mult)+' \n coords\n')
    for atom in openbabel.OBMolAtomIter(m):
        inp.write('  '+obet.GetSymbol(atom.GetAtomicNum())+' '+str(atom.GetX())+' '+str(atom.GetY())+' '+str(atom.GetZ())+' \n')
    inp.write(' end\nend\n')
    inp.write('%geom\n MaxIter 200\n end\n')
    inp.close()
    db['calculationSetup'] = { "molecularSpinMultiplicity" : mult, "charge" : charge, "numberOfElectrons" : get_n_electrons(m), "waveFunctionTheory" : "PM3" } 
    return(input_name) 
#
# Run Orca with generated input
#

def calculate_energy(fname, np):
    outfile = fname.strip('inp')+'out'
    out = open(outfile,'w')
    if np > 1:
        subprocess.call(['mpirun.openmpi', '-np', str(np), \
         '../../orca_3_0_3/orca',myInputFileName],stdout=myOutput)
    else:
        subprocess.call(['../../orca_3_0_3/orca',fname],stdout=out)
    return(outfile) 
#
# Get energy for given structure and charge
#

def get_energy(m):
    fname = create_orca_input(m)
    output_fname = calculate_energy(fname, np)
    output = open(output_fname,'r')
    o = output.readline()
    mol = {}
    atoms = []
    egy = 0.0
    while o:
        if o.find('THE OPTIMIZATION HAS CONVERGED')>=0:
            while o: 
                if o.find('CARTESIAN COORDINATES')>=0:
                    o = output.readline()
                    o = output.readline()
                    while o:
                        atom = {}
                        atomName,xcoord,ycoord,zcoord = mySplit(o,4)
                        if atomName == None: 
                            break
                        else:
                            atom['elementSymbol'] = atomName
                            atom['cartesianCoordinates'] = \
                                    { 'value' : [xcoord, ycoord, zcoord], 'units' : 'Angstrom' }
                            atoms.append(atom)
                            o = output.readline()
                            while o:
                                if o.find('Total Energy')>=0:
                                    egy = o.split()[3]
                            break
                        o = output.readline()
                    break
                o = output.readline()
            break
        o = output.readline()
#
#  Store data in MongoDb
#
    db['_id'] = fname.strip('.inp')
    db['parentInchi'] = m.GetTitle(mol_string)
    db['childInchi'] = obc.WriteString(m)
    obc.SetOutFormat('smi')
    db['childSmiles'] = obc.WriteString(m)
    obc.SetOutFormat('inchi')
    db['molecularFormula'] = m.GetFormula()
    mol['atoms'] = atoms
    db['molecule'] = mol
    db['totalEnergy'] = {"value":egy, "units":"Hartree"}
    insertResult = mongoDB.metatlas.insert_one(db)
    return(egy)

def protonate_molecule():
    pass

def read_molecules_from_csv(fname):
    mols = {}
    with open(fname) as csvFile:
        csvReader = csv.reader(csvFile)
        for row in csvReader: 
            _,inchiString,inchiKey = row[0],row[1],row[2]
            mols[inchiKey] = inchiString
    return mols

def perform_work(mol_string):
    obc.ReadString(m,mol_string)
    # I need this line to have the correct number of hydrogens
    formula = m.GetFormula()
    #   print(m.GetFormula())
    #   print(m.NumAtoms())
    #   print(m.NumConformers())
    #   print(obc.WriteString(m))
    m.AddHydrogens()
    b.Build(m)
    print(m.GetFormula())
    print(m.NumAtoms())
    print(m.NumConformers())
    print(m.Has2D())
    print(m.Has3D())
    obc.WriteFile(m, 'xyz/'+formula+'.xyz')
    #   obf.Setup(m)
    #   obf.ConjugateGradients(1000,1.0e-6)
    #   obf.UpdateCoordinates(m)
    m.SetTitle(mol_string)
    neutralEnergy = get_energy(m)
    protonatedEnergy, deprotonatedEnergy = 0.0, 0.0
    protonatedAtom, deprotonatedAtom = 0, 0
    
    for atom in openbabel.OBMolAtomIter(m):
        if atom.IsHydrogen:# Deprotonate
            mm = openbabel.OBMol(m)
            print 'deprotonating atom', mm.GetFormula()
            mm.DeleteAtom(atom)
            mm.SetTotalCharge(m.GetTotalCharge()-1)
            try:
                egy = get_energy(mm)
            except:
                print "failed to assign deprot energy"
                
            if egy < deprotonatedEnergy:
                deprotonatedEnergy = egy
            for connectedAtom in openbabel.OBAtomAtomIter(atom):
                deprotonatedAtom = connectedAtom.GetIdx()
                print('deprotonation',egy)
        else:
            protonate_molecule()
    #    #
    #    # Protonate, still needs work to make this go
    #    #
    #           atom.IncrementImplicitValence()
    #           m.AddHydrogens(atom)
    #           totalCharge = m.GetTotalCharge()
    #           m.SetTotalCharge(totalCharge+1)
    #           egy = getEnergy(m)
    #           if egy < protonationEnergy:
    #             protonatedEnergy = egy
    #             protonatedAtom = atom.GetIdx()
    #           atomToDelete = m.getAtom(m.NumAtoms())
    #           m.DeleteAtom(atomToDelete)
    #           m.SetTotalCharge(totalCharge)
    #           print('protonation',egy)
    print('lowest (de)protonation',deprotonatedEnergy,protonatedEnergy)

#
# Main routine reads metatlas inchi and converts to xyz
# Optimize structure
# Find deprotonation sites
# Find protonation sites
#

if __name__ == "__main__":
    obc = openbabel.OBConversion()
    obc.SetInAndOutFormats('inchi','xyz')
    obet = openbabel.OBElementTable()
    obf = openbabel.OBForceField.FindForceField('UFF')
    b = openbabel.OBBuilder()
    m = openbabel.OBMol()

    mongoClient = pymongo.MongoClient()
    mongoDB = mongoClient.test    
    db = {} 
    csv_file = 'metatlas_inchi_inchikey.csv'

    mols = read_molecules_from_csv(csv_file)

    np =1
    for key, string in mols.iteritems():
        perform_work(string)
