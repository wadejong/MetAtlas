#
# Base code for calculating most likely protonation sites
#
# Using OpenBabel as the toolkit to generate protonated and deprotonated structures
#

import openbabel, uuid, csv, subprocess, pymongo, bson
  
obc = openbabel.OBConversion()
obc.SetInAndOutFormats('inchi','xyz')
b = openbabel.OBBuilder()
m = openbabel.OBMol()
obElementTable = openbabel.OBElementTable()
obf = openbabel.OBForceField.FindForceField('UFF')

#
# Goal is to store the data in database, here trying MongoDb
#

mongoClient = pymongo.MongoClient()
mongoDB = mongoClient.test    
mongoDocument = {}


def mySplit(myString, numberOfVariables):
  return (myString.split() + [None] * numberOfVariables)[:numberOfVariables]


def numberOfElectrons(myMolecule):
  electronCount = 0
  for myAtom in openbabel.OBMolAtomIter(myMolecule):
    electronCount += myAtom.GetAtomicNum()
  return(electronCount)

#
# From structure and charge, generate an Orca input
#

def createInput(myMolecule):
  myCharge = myMolecule.GetTotalCharge()
  myMult = (1 if (numberOfElectrons(myMolecule)+myCharge) % 2 == 0 else 2)
  inputFileName = 'scr/'+str(uuid.uuid4())+'.inp'
  inputFile = open(inputFileName,'w')
  inputFile.write('!PM3 Opt \n%coords \n  CTyp xyz\n')
  inputFile.write(' Charge '+str(myCharge)+'\n Mult '+str(myMult)+' \n coords\n')
  for myAtom in openbabel.OBMolAtomIter(myMolecule):
    inputFile.write('  '+obElementTable.GetSymbol(myAtom.GetAtomicNum())+' '+str(myAtom.GetX())+' '+str(myAtom.GetY())+' '+str(myAtom.GetZ())+' \n')
  inputFile.write(' end\nend\n')
  inputFile.write('%geom\n MaxIter 200\n end\n')
  inputFile.close()
  mongoDocument['calculationSetup'] = { "molecularSpinMultiplicity" : myMult, "charge" : myCharge, "numberOfElectrons" : numberOfElectrons(myMolecule), "waveFunctionTheory" : "PM3" } 
  return(inputFileName)

#
# Run Orca with generated input
#

def calculateEnergy(myInputFileName):
  myOutputFileName = myInputFileName.strip('inp')+'out'
  myOutput = open(myOutputFileName,'wb')
  subprocess.call(['../../orca_3_0_3/orca',myInputFileName],stdout=myOutput)
  #subprocess.call(['mpirun.openmpi', '-np', '4', '../../orca_3_0_3/orca',myInputFileName],stdout=myOutput)
  return(myOutputFileName)

#
# Get energy for given structure and charge
#

def getEnergy(myMolecule):
  mongoDocument = {}
  myInputFileName = createInput(myMolecule)
  myOutputFileName = calculateEnergy(myInputFileName)
  myOutput = open(myOutputFileName,'r')
  outputLine = myOutput.readline()
  moleculeObj = {}
  atomsObj = []
  myEnergy = 0.0
  while outputLine:
    if outputLine.find('THE OPTIMIZATION HAS CONVERGED')>=0:
      while outputLine: 
        if outputLine.find('CARTESIAN COORDINATES')>=0:
          outputLine = myOutput.readline()
          outputLine = myOutput.readline()
          while outputLine:
            atomObj = {}
            atomName,xCoord,yCoord,zCoord = mySplit(outputLine,4)
            if atomName == None: 
              break
            else:
              atomObj['elementSymbol'] = atomName
              atomObj['cartesianCoordinates'] = { 'value' : [xCoord, yCoord, zCoord], 'units' : 'Angstrom' }
              atomsObj.append(atomObj)
            outputLine = myOutput.readline()
          while outputLine:
            if outputLine.find('Total Energy')>=0:
               myEnergy = outputLine.split()[3]
               break
            outputLine = myOutput.readline()
          break
        outputLine = myOutput.readline()
      break
    outputLine = myOutput.readline()
#
#  Store data in MongoDb
#
  mongoDocument['_id'] = myInputFileName.strip('.inp')
  mongoDocument['parentInchi'] = m.GetTitle(inchiString)
  mongoDocument['childInchi'] = obc.WriteString(myMolecule)
  obc.SetOutFormat('smi')
  mongoDocument['childSmiles'] = obc.WriteString(myMolecule)
  obc.SetOutFormat('inchi')
  mongoDocument['molecularFormula'] = myMolecule.GetFormula()
  moleculeObj['atoms'] = atomsObj
  mongoDocument['molecule'] = moleculeObj
  mongoDocument['totalEnergy'] = {"value":myEnergy, "units":"Hartree"}
  insertResult = mongoDB.metatlas.insert_one(mongoDocument)
  return(myEnergy)

#
# Main routine reads metatlas inchi and converts to xyz
# Optimize structure
# Find deprotonation sites
# Find protonation sites
#

with open('metatlas_inchi_inchikey.csv') as csvFile:
 csvReader = csv.reader(csvFile)
 for csvRow in csvReader:
   speciesNumber,inchiString,inchiKey = csvRow[0],csvRow[1],csvRow[2]
   print(inchiString)
   obc.SetInAndOutFormats('inchi','xyz')
   obc.ReadString(m,inchiString)
# I need this line to have the correct number of hydrogens
   formula = m.GetFormula()
   print(m.GetFormula())
   print(m.NumAtoms())
   print(m.NumConformers())
   print(obc.WriteString(m))
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
   m.SetTitle(inchiString)
   neutralEnergy = getEnergy(m)
   protonatedEnergy, deprotonatedEnergy = 0.0, 0.0
   protonatedAtom, deprotonatedAtom = 0, 0
#
# Deprotonate
#
   for atom in openbabel.OBMolAtomIter(m):
     if atom.IsHydrogen:
       mm = openbabel.OBMol(m)
       mm.DeleteAtom(atom)
       mm.SetTotalCharge(m.GetTotalCharge()-1)
       egy = getEnergy(mm)
       if egy < deprotonationEnergy:
         deprotonationEnergy = egy
         for connectedAtom in openbabel.OBAtomAtomIter(atom):
           deprotonatedAtom = connectedAtom.GetIdx()
       print('deprotonation',egy)
     else:
#
# Protonate, still needs work to make this go
#
       atom.IncrementImplicitValence()
       m.AddHydrogens(atom)
       totalCharge = m.GetTotalCharge()
       m.SetTotalCharge(totalCharge+1)
       egy = getEnergy(m)
       if egy < protonationEnergy:
         protonationEnergy = egy
         protonatedAtom = atom.GetIdx()
       atomToDelete = m.getAtom(m.NumAtoms())
       m.DeleteAtom(atomToDelete)
       m.SetTotalCharge(totalCharge)
       print('protonation',egy)
   print('lowest (de)protonation',deprotonationEnergy,protonationEnergy)
